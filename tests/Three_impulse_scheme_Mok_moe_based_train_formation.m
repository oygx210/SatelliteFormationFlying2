clear all;
consts = startup_formation_control();

% Testing 3impulse scheme that utilizes mean oe difference described in 
% Reference: Mok et al, Impulsive Control of Satellite Formation Flying using Orbital Period Difference

%% Setting system parameters
% Launch vehicle (LV) orbit [sma[m], ecc[-], inc[rad], RAAN[rad], AOP[rad], MA[rad]]
LV.coe(1) = 700e3 + consts.rEarth;                              % [m] sma
LV.coe(2) = 4e-4;                                               % [-] ecc
LV.coe(3) = get_SSO_inclination(LV.coe(1), LV.coe(2), consts);  %[rad] inc
LV.coe(4) = 0;                                                  % [rad] RAAN
LV.coe(5) = 0;                                                  % [rad] AOP
LV.coe(6) = 0;                                                  % [rad] Mean anomaly
LV.rv = oe2rv(LV.coe, consts);

% Spacecraft parameters
spacecraft.DragArea = 0.0481;           % [m^2] cross-section area perpendicular to the incoming airflow
spacecraft.Cdrag = 2.2;                 % atmospheric drag coefficient
spacecraft.mass = 4;                    % [kg]

% formation parameters
rho = 1000; % radius of relative orbit
alpha = 0; % phase of the satellite at the relative orbit

formation.N_sats = 2;
formation.deployment_interval = 30; % seconds
formation.drifting_time = consts.day2sec/24; % time when satellites are orbiting without control
formation.ISD = 1000;
formation.ISD_acceptable_error = 200;
formation.geometry(:,1) = [0 0 0 0];
formation.geometry(:,2) = [rho sqrt(3)/2*rho 0 0];

%CubeSat deploler
deployer.dep2orb = [1 0 0;
                    0 1 0;
                    0 0 1];

deployer.dV_dir = [-1.18; 0; 0];      % [m/s]
deployer.dV_mean = 0;                 % [m/s]
deployer.dV_sigma = 0.07/3;           % [m/s]
deployer.dV_sigma_wide = 0.035/3;     % [m/s]

disp('Undocking');
[t_undocking, rv_ECI_undocking] = undocking(LV.rv, deployer, formation, spacecraft, consts);
disp([num2str(formation.N_sats),' satellites were successfully released']);
% look at the ISD

disp('Uncontrolled flying');
options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec_drifting, rv_ECI_drifting] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:formation.drifting_time, rv_ECI_undocking(:,end), options_precision);
rv_ECI_drifting = rv_ECI_drifting';

rv_ECI = [rv_ECI_undocking, rv_ECI_drifting];
t = [t_undocking; t_vec_drifting];


% ISD = vecnorm(rv_ECI_drifting(1:3,end) - rv_ECI_drifting(7:9,end));
plot_ISD(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, formation, 'h');
plot_mean_oe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'm');


rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,end), formation_geometry, consts);

for i = 1:formation.N_sats
    rv_ECI_required(:,i) = orb2ECI(rv_ECI(1:6,end), rv_orb_required(:,i));
end

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% 1. Determines difference in mean orbital elements(moe)
leader.coe = vecRV2OE(rv_ECI(1:6), consts); leader.moe = rv2moe(rv_ECI(1:6), consts);
follower.coe = vecRV2OE(rv_ECI(7:12), consts); follower.moe = rv2moe(rv_ECI(7:12), consts);

dM_moe = 2*atan(formation.ISD/2/leader.moe(1));

follower.moe_required = leader.moe;
follower.moe_required(7) = leader.moe(7) - dM_moe;

follower.delta_moe = follower.moe_required - follower.moe;

if (follower.moe(5) > pi && follower.moe_required(5) < pi)
    follower.delta_moe(5) = 2*pi - follower.moe(5) + follower.moe_required(5);
end
if (follower.moe(7) > pi && follower.moe_required(7) < pi)
    follower.delta_moe(7) = 2*pi - follower.moe(7) + follower.moe_required(7);
end

%2. Calculates dV for three impuslive maneuvers and argument of latitude for first correction
% dV1 corrects di,dRAAN at argument of latitude theta_h
% dV2 & dV3 corrects da, de, dAOP, dM at perigee and apogee respectively

maneuvers = calculate_3impulses(follower.moe, follower.delta_moe, consts);
maneuvers_dV = [maneuvers.dV1; maneuvers.dV2; maneuvers.dV3];

if maneuvers.theta_h < follower.coe(8) 
   t_span = (2*pi - follower.coe(8) + maneuvers.theta_h)*sqrt(follower.moe(1)^3/consts.muEarth);
else
   t_span = (maneuvers.theta_h - follower.coe(8))*sqrt(follower.moe(1)^3/consts.muEarth);
end

% waiting to apply the first implulse at theta
rv_init = rv_ECI;
[t_vec1, rv_ECI1] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
t = t_vec1;
rv_ECI1 = rv_ECI1';
rv_ECI = rv_ECI1;

maneuvers_t(1,1) = t(end);

follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV1], consts);
follower.coe = vecRV2OE(follower.rv_ECI, consts); follower.moe = rv2moe(follower.rv_ECI, consts);
rv_ECI(7:12,end) = follower.rv_ECI;

t_span = (2*pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);

rv_init = rv_ECI(:,end);
[t_vec2, rv_ECI2] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
rv_ECI2 = rv_ECI2';
rv_ECI = [rv_ECI, rv_ECI2];
t = [t; t_vec2 + t(end)];

maneuvers_t(2,1) = t(end);
follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV2], consts);
follower.coe = vecRV2OE(follower.rv_ECI,consts); follower.moe = rv2moe(follower.rv_ECI,consts);
rv_ECI(7:12,end) = follower.rv_ECI;

% time before apocenter considering that pericenter has slightly shifted after maneuvering from its previous location
if follower.coe(6) > pi
    t_span = (3*pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);
else
    t_span = (pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);
end

rv_init = rv_ECI(:,end);
[t_vec3, rv_ECI3] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
rv_ECI3 = rv_ECI3';
rv_ECI = [rv_ECI, rv_ECI3];
t = [t; t_vec3 + t(end)];

maneuvers_t(3,1) = t(end);
follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV3], consts);
rv_ECI(7:12,end) = follower.rv_ECI;

plot_mean_oe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'h');