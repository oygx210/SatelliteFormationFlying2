% Main script for formation control

clear all;
clf;
consts = startup_formation_control();

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

% Formation parameters
% formation.min_IPD = 500; % the required minimum inter pixel distance
% [~,~,c1,c2,alpha] = positions(formation.min_IPD, "Sk", consts.muEarth/LV.coe(1)^3);
% % figure;
% % plot(c1.*cos(alpha), c1.*sin(alpha), 'ok', 'MarkerSize', 12);
% % axis off
% 
% formation.N_sats = length(c1);
% c1 = c1';
% c2 = c2';
% alpha = alpha';
% c3 = zeros(1,formation.N_sats);
% 
% formation.geometry = [c1; c2; c3; alpha];
formation.geometry = [0 0 0 0
                      0 0 1000 0];
formation.geometry = formation.geometry';

formation.N_sats = size(formation.geometry, 2);
formation.tracking_error = 0.2; 
formation.deployment_interval = 30; % seconds
formation.drifting_time = consts.day2sec/24; % time when satellites are orbiting without control
formation.FF_time = consts.day2sec*1; % formation flying experiment time

% CubeSat deployer 
deployer.dep2orb = [1 0 0;
                    0 1 0;
                    0 0 1];

deployer.dV_dir = [-1.18; 0; 0];      % [m/s]
deployer.dV_mean = 0;                 % [m/s]
deployer.dV_sigma = 0.07/3;           % [m/s]
deployer.dV_sigma_wide = 0.035/3;     % [m/s]

%% Phase A - Undocking & uncontrolled flying

disp('Undocking');
[t_undocking, rv_ECI_undocking] = undocking(LV.rv, deployer, formation, spacecraft, consts);
disp([num2str(formation.N_sats),' satellites were successfully released']);
% look at the ISD

disp('Uncontrolled flying');
options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec_drifting, rv_ECI_drifting] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:formation.drifting_time, rv_ECI_undocking(:,end), options_precision);
rv_ECI_drifting = rv_ECI_drifting';


for i = 1:formation.N_sats

    rv_orb(:,i) = ECI2orb(rv_ECI_drifting(1:6,end), rv_ECI_drifting(i*6-5:i*6,end), consts);
    rv_orb_required(:,i) = get_rv_from_analytic_HCW_solution(rv_ECI_drifting(1:6,end), formation.geometry(:,i), consts);
    rv_ECI_required(:,i) = orb2ECI(rv_ECI_drifting(1:6,end), rv_orb_required(:,i), consts);
    
end

% sat1.moe = rv2moe(rv_ECI_drifting(1:6,end), consts);
% sat1.coe = vecRV2OE(rv_ECI_drifting(1:6,end), consts);
% sat1.moe_required = rv2moe(rv_ECI_required(:,1), consts);
% sat1.coe_required = vecRV2OE(rv_ECI_required(:,1), consts);
% sat1.delta_coe = sat1.coe_required - sat1.coe;
% sat1.delta_moe = sat1.moe_required - sat1.moe;
% 
% 
% sat2.moe = rv2moe(rv_ECI_drifting(7:12,end), consts);
% sat2.coe = vecRV2OE(rv_ECI_drifting(7:12,end), consts);
% sat2.moe_required = rv2moe(rv_ECI_required(:,2), consts);
% sat2.coe_required = vecRV2OE(rv_ECI_required(:,2), consts);
% sat2.delta_coe = sat2.coe_required - sat2.coe;
% sat2.delta_moe = sat2.moe_required - sat2.moe;
% 
% figure(1);
% plot3(rv_orb(1,:), rv_orb(2,:),rv_orb(3,:), 'or');
% hold on;
% plot3(rv_orb_required(1,:), rv_orb_required(2,:),rv_orb_required(3,:), 'ok');
% xlabel('tangential,m');
% ylabel('normal,m');
% zlabel('radial,m');
% legend('Current', 'Required');
% axis equal;

%% Phase B - Formation flying (deployment and maintenance)
disp('FF experiment start');
[t_vec_FF, rv_FF] = formation_flying(rv_ECI_drifting(:,end), formation.FF_time, consts, spacecraft, formation);

rv_ECI = [rv_ECI, rv_FF];
t_vec = [t_vec; t_vec_FF+t_vec(end)];



