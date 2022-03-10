% Research sript for testing maintenance regime
clear all;
consts = startup_formation_control();

%% Initial conditions and simulation parameters

global environment
environment = 'J2';

spacecraft.dry_mass = 15;                                       % [kg]
spacecraft.propellant_mass = 1; % kg - propellant mass 
spacecraft.thruster_Isp = 65; % seconds

load('MD_output_Paris.mat');

formation.coe = MD.target_orbit;
formation.orbit_epoch = MD.orbit_epoch;
formation.rv = oe2rv(formation.coe, consts);

ISD_projected = 5000;
formation.geometry(:,1) = zeros(1,4);
formation.geometry(:,2) = [ISD_projected ISD_projected 0 0]';

formation.reconfiguration_flag = 0;
formation.N_sats = size(formation.geometry,2);
formation.N_active_sats = formation.N_sats-1;
formation.fuel_level = ones(formation.N_active_sats,1)*spacecraft.propellant_mass;
formation.tracking_error = 0.2; 
formation.FF_time = consts.day2sec*3;

rv_ECI(1:6,1) = formation.rv;
rv_orb_sat2 = get_rv_from_analytic_HCW_solution(rv_ECI(1:6), formation.geometry(:,2), consts);

rv_ECI(7:12,1) = orb2ECI(rv_ECI(1:6), rv_orb_sat2, consts);

%% Formation flying (relative orbit maintenance)

% options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
% [t_vec, rv_ECI] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[0 formation.FF_time], rv_ECI, options_precision);
% rv_ECI = rv_ECI';

disp('FF experiment start');
tic;
[t_vec, rv_ECI, maneuvers, T_corrections_start, T_corrections_end] = formation_flying(rv_ECI(:,end), formation.FF_time, consts, spacecraft, formation);
toc;

for j = 1:size(rv_ECI,2)

    for i = 1:formation.N_active_sats
        rv_orb(:,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*(i+1)-5:6*(i+1),j), consts);
    end

    rv_orb_required(:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), formation.geometry(:,2), consts);    
end

error = vecnorm(rv_orb_required(1:3,:) - rv_orb(1:3,:))./vecnorm(rv_orb_required(1:3,:));

figure('Name', '\epsilon(t)');
plot(t_vec/60, error(1:end-1));
xlabel('simulation time, min');
ylabel('\epsilon');
hold on;

figure('Name', 'Relative motion');
plot3(rv_orb(1,:), rv_orb(2,:), rv_orb(3,:));
hold on;
plot3(rv_orb_required(1,:), rv_orb_required(2,:), rv_orb_required(3,:));
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

% ISD = vecnorm(rv_ECI(1:3,:) - rv_ECI(7:9,:));
% figure;
% plot(t_vec/60,ISD);
% 
% plot_mean_oe_diff(t_vec, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'm');
% 
% for i = 1:size(rv_ECI,2)
%     rv_ECI_req(:,i) = orb2ECI(rv_ECI(1:6,i), rv_orb_required(:,i), consts);
% end

% 1. Plot error in oe
% 2. Plot oe(t) with flag at maneuvers(specific burns)
% resurrect 4 impulse maneuver;
