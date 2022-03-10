clear all;

% This is to try deployment and reconf with the aid of continuous control
% 1. обновить орбиту и моменты для перестроений
% 2. решить задачу для одного спутника и перестроений с разных точек на орбите и на
% относительные орбиту с разными положениями на относительной орбите

%% Initial conditions and simulation parameters
consts = startup_formation_control();

global environment A B K T
environment = 'point mass';

T = 0; % global time

spacecraft.dry_mass = 17;                                       % [kg]
spacecraft.propellant_mass = 1; % kg - propellant mass 
spacecraft.thruster_Isp = 214; % seconds, Busek BGT-X5
spacecraft.thrust = 180e-3; % thrust of satellite propulsion [N], Busek BGT-X5
% spacecraft.thrust = 500e-3; % thrust of satellite propulsion [N], Busek BGT-X5
spacecraft.u_max = spacecraft.thrust / (spacecraft.dry_mass + spacecraft.propellant_mass); % max contol vector magnitude - unit thrust 

k_rev2rep = 14;
k_day2rep = 1;
% orbital elements
orbit_epoch_GD = datetime(2023, 1, 1, 0, 0, 0);
orbit_epoch_JD = juliandate(orbit_epoch_GD);

oe = zeros(6,1);
oe(1) = 900e3 + consts.rEarth;
oe(2) = 0; % ecc, 
oe(3) = 0;
orbit_epoch_GD = datetime(orbit_epoch_JD, 'convertfrom','juliandate');
oe(4) = 0;
oe(5) = 0; % AOP, deg
oe(6) = 0; % M, deg - actually argument of latitude
% the actual orbital elements should correspond to the optimal coverage

formation.coe = oe;
mean_motion = sqrt(consts.muEarth/formation.coe(1)^3);
formation.orbit_epoch = orbit_epoch_GD;
formation.rv = oe2rv(formation.coe, consts);

theta_min = deg2rad(10);
b = consts.rEarth / oe(1);
gamma_max = asin(b * cos(theta_min));
beta_max = pi/2 - gamma_max - theta_min;

d_max = oe(1) / cos(theta_min) * sin(pi/2 - gamma_max - theta_min);
ISD_min = deg2rad(1/60) * d_max;

rho = ISD_min * 4.5;
rho = 0.4e4;

formation.IPD_min = ceil(ISD_min); 
formation.ISD_safe = 30; % meters

demonstration{1,1}.HCW_constants = [[0;0;0;0], [rho; rho; 0; (pi/2 + pi/4)]];
rv_orb = get_rv_from_analytic_HCW_solution(formation.rv, demonstration{1,1}.HCW_constants(:,2), consts);

rhs_mag = sqrt((-2 * mean_motion * rv_orb(6))^2 + (-mean_motion^2 * rv_orb(2))^2 + (2 * mean_motion * rv_orb(4) + 3 * mean_motion^2 * rv_orb(3))^2);
u2rhs_max = spacecraft.u_max / rhs_mag;

demonstration{1,1}.deployment_time = orbit_epoch_GD;
demonstration{1,1}.reconfiguration_time = orbit_epoch_GD + seconds(days(1));
demonstration{1,1}.demo_time = [42000 : 43000];
demonstration{1,1}.demo_midpoint = orbit_epoch_GD + seconds(42500);

formation.N_sats = size(demonstration{1,1}.HCW_constants,2);
formation.N_active_sats = formation.N_sats-1;
formation.fuel_level = ones(formation.N_active_sats,1)*spacecraft.propellant_mass;
% formation.tracking_error = round(formation.IPD_min/10); % tracking error in meters
formation.tracking_error = 1;
formation.tracking_error_rho = 1; % tracking error in meters
formation.tracking_error_v = 0.01; % tracking error in meters

% LQR definition
switch environment
    case 'point mass'
        n = sqrt(consts.muEarth/formation.coe(1)^3);
        C = [0 0 0
             0 -n^2 0 
             0 0 3*n^2];
        D = [0 0 -2*n
             0 0 0
             2*n 0 0];
        A = [[zeros(3), eye(3)]; [C,D]];
    case 'J2'
        n = sqrt(consts.muEarth/formation.coe(1)^3);
        s = 3*consts.J2*consts.rEarth_equatorial^2/2/formation.coe(1)^2*(1+3*cos(2*formation.coe(3)));
        c = sqrt(1 + s);

        C = [0 0 0
             0 -(3*c^2-2)*n^2 0 
             0 0 (5*c^2-2)*n^2];
        D = [0 0 -2*n*c
             0 0 0
             2*n*c 0 0];
        A = [[zeros(3), eye(3)]; [C,D]];
end

B = [zeros(3); eye(3)];

Q = diag([1e-6, 1e-6, 1e-6, 1e-2, 1e-2, 1e-2]);

R = diag([1; 1; 1]);

[K,~,~] = lqr(A,B,Q,R);

% formation satellites are at the same orbit
for i = 1:formation.N_sats
    rv_ECI(i*6-5:i*6,1) = formation.rv;
end

orbit_period = 2*pi/n;
T_maintenance = orbit_period * 2;
T_rude_control = T_maintenance;

formation.geometry = demonstration{1,1}.HCW_constants;
[t_vec, rv_ECI, maneuvers, formation_fuel_level] = continuous_control(rv_ECI, T_maintenance, T_rude_control, consts, spacecraft, formation, 2);


