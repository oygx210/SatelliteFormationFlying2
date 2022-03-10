% Main script for formation control

clear all;
consts = startup_formation_control();

%% Initial conditions and simulation parameters

% Spacecraft parameters - 12U CubeSat biggest side area is considered
spacecraft.dry_mass = 15;                                       % [kg]
spacecraft.propellant_mass = 1; % kg - propellant mass 
spacecraft.thruster_Isp = 65; % seconds

global environment
environment = 'J2';

load('MD_output_Paris.mat');

formation.coe = MD.target_orbit;
formation.orbit_epoch = MD.orbit_epoch;
formation.rv = oe2rv(formation.coe, consts);
formation.ISD_min = MD.ISD_min; 
formation.ISD_safe = 30; % meters

formation.image_attitude_at_reconf1 = MD.image_attitude_at_reconf(1);
formation.geometry1 = orbital_configuration('AA', 'PCO', formation.ISD_min, formation.image_attitude_at_reconf1);

formation.geometry = formation.geometry1;

formation.image_attitude_at_reconf2 = MD.image_attitude_at_reconf(2);
formation.geometry2 = orbital_configuration("IAA_logo", "PCO", formation.ISD_min, formation.image_attitude_at_reconf2);
formation.reconfiguration_time = seconds(MD.reconfiguration_time(2) - formation.orbit_epoch); % at the time the reconfiguration from geometry1 to geometry2 should happen

formation.reconfiguration_flag = 2;

formation.N_sats = size(formation.geometry,2);
formation.N_active_sats = formation.N_sats-1;
formation.fuel_level = ones(formation.N_active_sats,1)*spacecraft.propellant_mass;
formation.tracking_error = 0.1; 
formation.deployment_interval = 30; % seconds
formation.drifting_time = 1000; % time when satellites are orbiting without control
formation.FF_time = consts.day2sec*1 - 2*3600; % formation flying experiment time

% CubeSat deployer 
deployer.dep2orb = [1 0 0;
                    0 1 0;
                    0 0 1];

deployer.dV_dir = [-1.18; 0; 0];      % [m/s]
deployer.dV_mean = 0;                 % [m/s]
deployer.dV_sigma = 0.07;           % [m/s]
deployer.dV_sigma_wide = 0.035;     % [m/s]

%% Undocking phase

% Experiment_option: 1 - undocking with random deployment errors and uncontrolled drifting
%                    2 - the satellites have the same sma (idealized case)

Experiment_option = 2;

switch Experiment_option
    case 1   
        % Undocking & uncontrolled flying

        disp('Undocking');
        [t_undocking, rv_ECI_undocking] = undocking(formation.rv, deployer, formation, spacecraft, consts);
        disp([num2str(formation.N_sats),' satellites were successfully released']);
        % look at the ISD
        
        for i = 1:formation.N_sats
            oe_undocking(:,i) = rv2oe(rv_ECI_undocking(i*6-5:i*6,end), consts);
        end
                        
        disp('Uncontrolled flying');
        options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [t_vec_drifting, rv_ECI_drifting] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),0:formation.drifting_time, rv_ECI_undocking(:,end), options_precision);
        
        rv_ECI_drifting = rv_ECI_drifting';

        for i = 1:formation.N_sats
            rv_orb(:,i) = ECI2orb(rv_ECI_drifting(1:6,end), rv_ECI_drifting(i*6-5:i*6,end), consts);
        end
        
        figure;
        plot3(rv_orb(1,:), rv_orb(2,:),rv_orb(3,:), 'ok', 0,0,0,'or');
        title('Orbital configuration after commissioning period');
        xlabel('tangential,m');
        ylabel('normal,m');
        zlabel('radial,m');
        axis equal;

    case 2        

        for i = 1:formation.N_sats
            rv_ECI_drifting(i*6-5:i*6,1) = formation.rv;
        end
end

%% Formation flying (deployment, maintenance, and reconfiguration)
disp('FF experiment start');
[t_vec_FF, rv_FF, maneuvers] = formation_flying(rv_ECI_drifting(:,end), formation.FF_time, consts, spacecraft, formation);

%% Visualization
visualization_flag = 1 ;

switch visualization_flag
    case 0
        
         rv_ECI = rv_FF;
         t_vec = t_vec_FF;
    
    case 1
    
         rv_ECI = rv_FF;
         t_vec = t_vec_FF;

    for j = 1:size(rv_ECI,2)
        for i = 1:formation.N_sats
            rv_orb(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
        end
    end

    for i = 1:6081
        rv_orb_required(:,:,i) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,i), formation.geometry, consts);
    end

    for i = 6081:size(rv_ECI,2)
        rv_orb_required(:,:,i) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,i), formation.geometry2, consts);
    end

    fig = figure;
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    xlim([-10000 10000]);
    ylim([-10000 10000]);
    zlim([-10000 10000]);
    pbaspect([1 1 1]);
    title('Orbital reference frame');
    view(-45,-75);
    hold on;

    for i = 1:size(rv_orb,3)/30

        a = plot3(rv_orb_required(1,:,i*20), rv_orb_required(2,:,i*20), rv_orb_required(3,:,i*20), 'or');
        hold on;
        b = plot3(rv_orb(1,:,i*20), rv_orb(2,:,i*20), rv_orb(3,:,i*20), 'ok');
        legend('Required configuration', 'Current configuration');

        drawnow;

    %   Take a Snapshot
        movieVector(i) = getframe(fig);   %manually specify getframe region    

        delete(b);
        delete(a)
    end
    plot3(rv_orb(1,:,end), rv_orb(2,:,end), rv_orb(3,:,end), 'ok');
    plot3(rv_orb_required(1,:,end), rv_orb_required(2,:,end), rv_orb_required(3,:,end), 'or');

    myWriter = VideoWriter('Test_optimized','MPEG-4');   %create an .mp4 file
    myWriter.FrameRate = 40;

    %   Open the VideoWriter object, write the movie, and close the file
    open(myWriter);
    writeVideo(myWriter, movieVector);
    close(myWriter);
end

% rv_ECI = [rv_ECI, rv_FF];
% t_vec = [t_vec; t_vec_FF+t_vec(end)];

