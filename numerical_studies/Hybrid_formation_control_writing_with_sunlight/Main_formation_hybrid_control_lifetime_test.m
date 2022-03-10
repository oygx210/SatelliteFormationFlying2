clear all;

%% Initial conditions and simulation parameters
consts = startup_formation_control();
load('MD_output_JASR.mat');

global environment A B K T
environment = 'J2';

T = 0; % global time

spacecraft.dry_mass = 17;                                       % [kg]
spacecraft.propellant_mass = 1; % kg - propellant mass 
spacecraft.thruster_Isp = 214; % seconds, Busek BGT-X5
spacecraft.thrust = 180e-3; % thrust of satellite propulsion [N], Busek BGT-X5
spacecraft.u_max = spacecraft.thrust / (spacecraft.dry_mass + spacecraft.propellant_mass); % max contol vector magnitude - unit thrust 

formation.coe = MD.target_orbit;
mean_motion = sqrt(consts.muEarth/formation.coe(1)^3);
formation.orbit_epoch = MD.orbit_epoch;
formation.rv = oe2rv(formation.coe, consts);
formation.IPD_min = MD.IPD_min; 
formation.ISD_safe = 30; % meters
formation.final_orbit_epoch = MD.final_orbit_epoch;
formation.N_sats = size(MD.HCW_constants_demo1,2);
formation.N_active_sats = formation.N_sats-1;
formation.fuel_level = ones(formation.N_active_sats,1)*spacecraft.propellant_mass;
formation.tracking_error = 1; % tracking error in meters

T = 0;

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

Q = diag([1e-7, 1e-7, 1e-7, 1e-9, 1e-9, 1e-9]);
R = diag([1; 1; 1]);

[K,~,~] = lqr(A,B,Q,R);

mode = 'recon_test';
if mode == 'maint_test'
    image_index = 2;
    maintenance_duration = days(30);
    demonstration{1,1}.HCW_constants = MD.HCW_constants_demo2;
    demonstration{1,1}.deployment_time = MD.orbit_epoch;
    demonstration{1,1}.reconfiguration_time = demonstration{1,1}.deployment_time + maintenance_duration;
    demonstration{1,1}.demo_time = [seconds(MD.demonstration_time(image_index,1) - MD.orbit_epoch),seconds(MD.demonstration_time(image_index,2) - MD.orbit_epoch)];
    demonstration{1,1}.demo_midpoint = MD.demo2_midpoint;
    demonstration{1,1}.Cost_matrix_dV = zeros(50);
elseif mode == 'recon_test'    
    mode2 = 2;
    if mode2 == 1
        demonstration{1,1}.HCW_constants = MD.HCW_constants_demo1;
        demonstration{1,1}.deployment_time = MD.orbit_epoch;
        demonstration{1,1}.reconfiguration_time = demonstration{1,1}.deployment_time + hours(12);
        demonstration{1,1}.demo_time = [seconds(MD.demonstration_time(1,1) - MD.orbit_epoch),seconds(MD.demonstration_time(1,2) - MD.orbit_epoch)];        
        demonstration{1,1}.demo_midpoint = MD.demo1_midpoint;
        demonstration{1,1}.Cost_matrix_dV = zeros(50);

        demonstration{2,1}.HCW_constants = MD.HCW_constants_demo2;
        demonstration{2,1}.deployment_time = demonstration{1,1}.reconfiguration_time;
        demonstration{2,1}.reconfiguration_time = demonstration{2,1}.deployment_time + hours(12);
        demonstration{2,1}.demo_time = [seconds(MD.demonstration_time(2,1) - MD.orbit_epoch),seconds(MD.demonstration_time(2,2) - MD.orbit_epoch)];        
        demonstration{2,1}.demo_midpoint = MD.demo2_midpoint;
        load('Cost_matrix_JASR.mat');
        Cost_matrix_maintenance = ones(size(MD.HCW_constants_demo1,2)-1).*cost_matrix_dV.maintenance_cost';
        Cost_matrix_reconfiguration = cost_matrix_dV.reconfiguration_cost; 
        Cost_matrix_dV = Cost_matrix_reconfiguration + Cost_matrix_maintenance;
        demonstration{2,1}.Cost_matrix_dV = Cost_matrix_dV;

        for i = 3:20
            demonstration{i,1}.deployment_time = demonstration{i-1,1}.reconfiguration_time;
            demonstration{i,1}.reconfiguration_time = demonstration{i,1}.deployment_time + hours(12);
            demonstration{i,1}.demo_time = seconds(demonstration{i,1}.deployment_time - formation.orbit_epoch) + [seconds(hours(6)), seconds(hours(6))+100];

            if mod(i,2) == 1
                demonstration{i,1}.HCW_constants = MD.HCW_constants_demo1;
            else
                demonstration{i,1}.HCW_constants = MD.HCW_constants_demo2;
            end   
        end
        
        for i = 3:length(demonstration)
            cost_matrix_dV = get_cost_matrix(demonstration{i-1,1}, demonstration{i,1}, formation, spacecraft, consts);
            Cost_matrix_maintenance = ones(formation.N_active_sats).*cost_matrix_dV.maintenance_cost';
            Cost_matrix_reconfiguration = cost_matrix_dV.reconfiguration_cost; 
            Cost_matrix_dV = Cost_matrix_reconfiguration + Cost_matrix_maintenance;
            demonstration{i,1}.Cost_matrix_dV = Cost_matrix_dV;
        end
    save('C:\GoogleDrive\SatelliteFormationFlying\data\lifetime_ASR_20demos', 'demonstration');

    elseif mode2 == 2
        load('lifetime_ASR_20demos.mat');
    end
   
end
for i = 1:7
    demonstration_new{i,1} = demonstration{i,1};
end
demonstration = demonstration_new;   
T = 0;

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

Q = diag([1e-7, 1e-7, 1e-7, 1e-9, 1e-9, 1e-9]);
R = diag([1; 1; 1]);

[K,~,~] = lqr(A,B,Q,R);

for i = 1:formation.N_sats
    rv_ECI(i*6-5:i*6,1) = formation.rv;
end
orbit_period = 2*pi/n;

%% Formation flying dynamics and control according to the mission scenario
disp('Image demonstration mission simulation started');
tic;
[t_vec, rv_ECI, maneuvers, post_coorection_maneuvers, fuel_consumption, t_events, formation_state, HCW_constants_assigned, match_matrix] = formation_hybrid_control_structured(rv_ECI(:,end), demonstration, formation, spacecraft, consts);
toc;

for i = 1:size(rv_ECI,2)
    oe(:,i) = rv2oe(rv_ECI(1:6,i), consts);
end

for i = 1:size(t_events,1)
    
    [~,index_events(i,1)] = min(abs(t_vec - t_events(i,1)));
    [~,index_events(i,2)] = min(abs(t_vec - t_events(i,2)));
    
end

[~,demo1_index] = min(abs(t_vec - seconds(demonstration{1, 1}.demo_midpoint - formation.orbit_epoch)));
[~,demo2_index] = min(abs(t_vec - seconds(demonstration{2, 1}.demo_midpoint - formation.orbit_epoch)));

u1_calculated = rad2deg(oe(8, demo1_index));
u2_calculated = rad2deg(oe(8, demo2_index));

[~,t_step_reconf] = min(abs(t_vec - t_events(4,2)));

disp('Postprocessing: Converting formation state vector from ECI to the orbital reference frame');
tic;
for j = 1:t_step_reconf

    rv_orb_required(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), HCW_constants_assigned(:,:,1), consts);

    for i = 1:formation.N_sats
        rv_orb(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
        deltaR(i,j) = vecnorm(rv_orb_required(1:3,i,j) - rv_orb(1:3,i,j));
        deltaV(i,j) = vecnorm(rv_orb_required(4:6,i,j) - rv_orb(4:6,i,j));
        c1(i,j) = rv_orb(4,i,j)/mean_motion + 2*rv_orb(3,i,j);

    end
end

for j = t_step_reconf+1:size(rv_ECI,2)

    rv_orb_required(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), HCW_constants_assigned(:,:,2), consts);

    for i = 1:formation.N_sats
        rv_orb(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
        deltaR(i,j) = vecnorm(rv_orb_required(1:3,i,j) - rv_orb(1:3,i,j));
        deltaV(i,j) = vecnorm(rv_orb_required(4:6,i,j) - rv_orb(4:6,i,j));
        c1(i,j) = rv_orb(4,i,j)/mean_motion + 2*rv_orb(3,i,j);        
    end
end
toc;