%% The script is used to compute cost matrix C for reconfiguration from one
% orbital configuration to another and its further maintenance for a specific period of time 

% Input:
% 1. Formation's target orbit state vector
% 2. Current orbital configuration
% 3. Target orbital configuration
% 6. Time for maintenance of the new orbital configuration
 
% Output: 
% 1. Cost matrix for reconfiguration
% 2. Maintenance cost for all relative trajectories 

% Sctipt outline: 
% 1. Calculate target orbit at reconf time
% 2. Based on the HCW constants find ECI state vectors of formation
% satellites
% 3. Calculate reconfiguration matrix taking intro account impulsive and
% continuous post control 
% 4. Compute maintenance cost

clear all;

%% Inputs for cost matrix calculation data
consts = startup_formation_control();
load('MD_output_Acta_Paper.mat');
n = sqrt(consts.muEarth/MD.target_orbit(1)^3);

spacecraft.dry_mass = 17;                                       % [kg]
spacecraft.propellant_mass = 1; % kg - propellant mass 
spacecraft.thruster_Isp = 214; % seconds, Busek BGT-X5
spacecraft.thrust = 100e-3; % thrust of satellite propulsion [N], Busek BGT-X5
spacecraft.u_max = spacecraft.thrust / (spacecraft.dry_mass + spacecraft.propellant_mass); % max contol vector magnitude - unit thrust 

formation.IPD_min = MD.IPD_min; 
formation.ISD_safe = 30; % meters
formation.final_orbit_epoch = MD.final_orbit_epoch;
formation.N_sats = size(MD.HCW_constants_demo1,2);
formation.N_active_sats = formation.N_sats-1;
formation.fuel_level = ones(formation.N_active_sats,1)*spacecraft.propellant_mass;
formation.tracking_error = round(formation.IPD_min/10); % tracking erro in meters

global environment A B K T
environment = 'J2';
T = 0;
% LQR definition
switch environment
    case 'point mass'
        n = sqrt(consts.muEarth/MD.target_orbit(1)^3);
        C = [0 0 0
             0 -n^2 0 
             0 0 3*n^2];
        D = [0 0 -2*n
             0 0 0
             2*n 0 0];
        A = [[zeros(3), eye(3)]; [C,D]];
    case 'J2'
        n = sqrt(consts.muEarth/MD.target_orbit(1)^3);
        s = 3*consts.J2*consts.rEarth_equatorial^2/2/MD.target_orbit(1)^2*(1+3*cos(2*MD.target_orbit(3)));
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

Q = diag([1e-7, 1e-7, 1e-7, 1e-9, 1e-9, 1e-9]);
R = diag([1; 1; 1]);

[K,~,~] = lqr(A,B,Q,R);

reconfiguration_time_12 = MD.reconfiguration_time(2);
reconfiguration_time_21 = MD.orbit_epoch + days(1);

orbit_period = 2*pi/n;

u1 = mod(n*seconds(reconfiguration_time_12 - MD.orbit_epoch), 2*pi);
u2 = mod(n*seconds(reconfiguration_time_21 - MD.orbit_epoch), 2*pi);

target_orbit_12 = MD.target_orbit; target_orbit_21 = MD.target_orbit;
target_orbit_12(6) = u1; target_orbit_21(6) = u2; 
target_orbit_12 = oe2rv(target_orbit_12, consts);
target_orbit_21 = oe2rv(target_orbit_21, consts);

Cost_matrix_12 = get_cost_matrix(target_orbit_12, MD.HCW_constants_demo1, MD.HCW_constants_demo2, seconds(days(1)/2 - seconds(orbit_period*3)), formation, spacecraft, consts);
T = 0;
Cost_matrix_21 = get_cost_matrix(target_orbit_21, MD.HCW_constants_demo2, MD.HCW_constants_demo1, seconds(days(1)/2 - seconds(orbit_period*3)), formation, spacecraft, consts);

Cost_matrices.HCW_configuration1 = MD.HCW_constants_demo1;
Cost_matrices.HCW_configuration1 = MD.HCW_constants_demo2;
Cost_matrices.target_orbit = MD.target_orbit;
Cost_matrices.orbit_epoch = MD.orbit_epoch;
Cost_matrices.Cost_matrix_12 = Cost_matrix_12;
Cost_matrices.Cost_matrix_21 = Cost_matrix_21;


%% saving data
save('C:\GoogleDrive\SatelliteFormationFlying\data\Cost_matrices_with_inputs','Cost_matrices');

function cost_matrix = get_cost_matrix(target_orbit, HCW_constants_current, HCW_constants_required, maintenance_time, formation, spacecraft, consts)
tic;
rv_orb_current = get_rv_from_analytic_HCW_solution(target_orbit, HCW_constants_current, consts);

for i = 1:formation.N_sats
    rv_ECI_current(i*6-5:i*6,1) = orb2ECI(target_orbit, rv_orb_current(:,i), consts);
end

reconfiguration_matrix_dV = get_reconfiguration_matrix(rv_ECI_current, HCW_constants_required, formation, spacecraft, consts);

% maintenance cost
rv_orb_required = get_rv_from_analytic_HCW_solution(target_orbit, HCW_constants_required, consts);
formation.geometry = HCW_constants_required;

for i = 1:formation.N_sats
    rv_ECI_required(i*6-5:i*6,1) = orb2ECI(target_orbit, rv_orb_required(:,i), consts);
end
[~, ~, maneuvers_maintenance] = maintenance(rv_ECI_required, maintenance_time, 0, formation, spacecraft, consts);

cost_matrix.reconfiguration = reconfiguration_matrix_dV;
cost_matrix.maintenance = maneuvers_maintenance;
toc;
end