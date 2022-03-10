clear all;
consts = startup_formation_control_func();

%% Readme for Swarm_dispersion_BC (Boundary conditions)

% The script is used for experiments on the swarm dispersion during the mission
% The main goal is to define deployer parameters to ensure required intersatellite distance

% The experiment is made for two similar 3U CubeSats
% The force model takes into account J2 effect and atmospheric drag

% The CubeSat swarm is to be deployed from International Space Station 

% Experiment plan.
% 1. We generate set of deployment velocities for a pair of satellites to
% consider boundary conditions for formation dispersion;
% 2. We propagate formation motion (formation consists of 2 sats) for 2
% months which is considered as mission duration according to current ConOps
% 3. Based on the simulations we define requirement for deployment velocity
% which ensure correct ISL distance for the whole mission duration.

% All vectors are defined as columns vector 

%% Spacecraft parameters
spacecraft.DragArea = 0.0481;            % [m^2] cross-section area (worst case) perpendicular to the incoming airflow
spacecraft.Cdrag = 2.2;                  % atmospheric drag coefficient
spacecraft.mass = 4;                     % [kg] CubeSat 3U mass
spacecraft.max_ISL = 1000*consts.km2m;   % [m] maximum distance for Inter-satellite link

%% CubeSat deployer parameters
deployer.dep2orb = [1 0 0;
                    0 1 0;
                    0 0 1];

% Varibles in the experiment
deployer.dV_dir = [-1; 0; 0];         % [m/s]
deployer.dV_mean = 0;                 % [m/s]
deployer.dV_3sigma = 0.02;            % [m/s]
deployer.dV_3sigma_wide = 0.25;      % [m/s]

%% Launch vehicle (LV) orbit - ISS orbit
LV.oe = zeros(6,1);
% OE = [sma[km], ecc[-], inc[rad], RAAN[rad], AOP[rad], MA[rad]]
LV.oe(1) = 400e3 + consts.rEarth;       % [km] sma
LV.oe(2) = 0;                           % ecc
LV.oe(3) = 52*pi/180;                          % [rad] inc
LV.oe(4) = 0;                           % [rad] RAAN
LV.oe(5) = 0;                           % [rad] AOP
LV.oe(6) = 0;                           % [rad] Mean anomaly

LV.rv = oe2rv(LV.oe, consts);
LV.mean_motion = sqrt(consts.muEarth/LV.oe(1)^3);

%% Boundary conditions simulation
N_conditions = 8;

% Simulation parameters
simulation_time = 1*consts.day2sec; % s
stopevent = @(t, rv) connection_is_lost(t, rv, spacecraft);
options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',  stopevent);
tspan = 0:simulation_time;

% Set of deployment conditions
Target.vec_error = [-1*deployer.dV_3sigma deployer.dV_3sigma_wide 0    
                    -1*deployer.dV_3sigma 0 deployer.dV_3sigma_wide
                    -1*deployer.dV_3sigma -deployer.dV_3sigma_wide 0
                    -1*deployer.dV_3sigma 0 -deployer.dV_3sigma_wide
                     deployer.dV_3sigma deployer.dV_3sigma_wide 0
                     deployer.dV_3sigma 0 deployer.dV_3sigma_wide 
                     deployer.dV_3sigma -deployer.dV_3sigma_wide 0
                     deployer.dV_3sigma 0 -deployer.dV_3sigma_wide];
                 
Target.vec_error = Target.vec_error';
Target.rv_orb_deployment = [zeros(1,N_conditions); zeros(1,N_conditions); zeros(1,N_conditions); 
                            deployer.dV_dir + Target.vec_error];                        

Target.rv_ECI_deployment = orb2ECI(LV.rv, Target.rv_orb_deployment, consts);
Target.dV_deployment_mag = sqrt(dot(Target.rv_orb_deployment(4:6,:),Target.rv_orb_deployment(4:6,:)));
Target.dV_deployment_theta = atand(sqrt(Target.rv_orb_deployment(5,:).^2 + Target.rv_orb_deployment(6,:).^2) ./ Target.rv_orb_deployment(4,:));
Target.dV_deployment_phi = atan2(Target.rv_orb_deployment(6,:),Target.rv_orb_deployment(5,:));

Chaser.vec_error = -1*Target.vec_error;
Chaser.rv_orb_deployment = [zeros(1,N_conditions); zeros(1,N_conditions); zeros(1,N_conditions); 
                            deployer.dV_dir + Chaser.vec_error];                        

Chaser.rv_ECI_deployment = orb2ECI(LV.rv, Chaser.rv_orb_deployment, consts);
Chaser.dV_deployment_mag = sqrt(dot(Chaser.rv_orb_deployment(4:6,:),Chaser.rv_orb_deployment(4:6,:)));
Chaser.dV_deployment_theta = atand(sqrt(Chaser.rv_orb_deployment(5,:).^2 + Chaser.rv_orb_deployment(6,:).^2) ./ Chaser.rv_orb_deployment(4,:));
Chaser.dV_deployment_phi = atan2(Chaser.rv_orb_deployment(6,:),Chaser.rv_orb_deployment(5,:));

figure(1);
% Check boundary conditions
for i = 1:8
    subplot(2,4,i);
    polarplot(Target.dV_deployment_phi(i), Target.dV_deployment_theta(i),'ok','MarkerSize', 12, 'MarkerFaceColor', 'k');
    hold on;
    polarplot(Chaser.dV_deployment_phi(i), Chaser.dV_deployment_theta(i),'or','MarkerSize', 12, 'MarkerFaceColor', 'r');
end
legend('Targets deployment dV', 'Chasers deployment dV')

rv_ECI_deployment = [Target.rv_ECI_deployment; Chaser.rv_ECI_deployment];

tic;
f = waitbar(0,'Please wait...');
for i = 1:N_conditions

    sol{1,i} = ode45(@(t, rv) J2_atmo_2sats(t, rv, consts, spacecraft), tspan, rv_ECI_deployment(:,i), options);

    waitbar(i/N_conditions, f, 'Please wait...');
end
close(f);
toc;

for i = 1:N_conditions
    Mission_duration(i) = sol{1,i}.x(end);
    Intersatellite_dist{i,:} =  sqrt(dot(sol{1,i}.y(1:3,:) - sol{1,i}.y(7:9,:),sol{1,i}.y(1:3,:) - sol{1,i}.y(7:9,:)));
    t_vec{i,:} = sol{1,i}.x(:);
end

[t_min, index_t_min] = min(Mission_duration);
ISL_max = zeros(N_conditions,1);

if t_min < simulation_time
    deployer.quality = 0;
    Swarm_LT = t_min/consts.day2sec;
    tq = 1:floor(t_min);

    for i = 1:N_conditions
        Intersatellite_dist_interp{i,:} =  interp1(sol{1,i}.x, Intersatellite_dist{i,:}, tq, 'PCHIP');
        Intersatellite_dist_vec(i,:) = Intersatellite_dist_interp{i,:};
        ISL_max(i) = Intersatellite_dist_vec(i,end);
    end
    
    
%     deval;
    Intersatellite_dist_mean = mean(Intersatellite_dist_vec);
    Intersatellite_dist_std = std(Intersatellite_dist_vec);

    figure;
    for i = 1:N_conditions
        plot(tq/consts.day2sec, Intersatellite_dist_vec(i,:)/consts.km2m);
        hold on;
    end
    
    plot(tq/consts.day2sec, Intersatellite_dist_mean/1000, '--k');
    hold on
    ErrorBarStep = ceil(t_min/60);
    errorbar(tq(1:ErrorBarStep:end)/consts.day2sec, Intersatellite_dist_mean(1:ErrorBarStep:end)/consts.km2m, Intersatellite_dist_std(1:ErrorBarStep:end)/consts.km2m);
    xlabel('time, days');
    ylabel('ISL, km');
    legend(['dV = ',num2str(abs(deployer.dV_dir(1)),2),' m/s; 3\sigma_{||} = ',num2str(deployer.dV_3sigma,2), ' m/s; 3\sigma_{\perp} = ', num2str(deployer.dV_3sigma_wide,2), ' m/s']);
    title('Swarm dispersion simulation with boundary conditions');

else  
    deployer.quality = 1;
    Swarm_LT = t_min/consts.day2sec;
    tq = 1:floor(t_min);

    for i = 1:N_conditions
        Intersatellite_dist_interp{i,:} =  interp1(sol{1,i}.x, Intersatellite_dist{i,:}, tq, 'PCHIP');
        Intersatellite_dist_vec(i,:) = Intersatellite_dist_interp{i,:};
        ISL_max(i) = Intersatellite_dist_vec(i,end);
    end
    
    Intersatellite_dist_mean = mean(Intersatellite_dist_vec);
    Intersatellite_dist_std = std(Intersatellite_dist_vec);

    figure;
    for i = 1:N_conditions
        plot(tq/consts.day2sec, Intersatellite_dist_vec(i,:)/consts.km2m);
        hold on;
    end
    
    plot(tq/consts.day2sec, Intersatellite_dist_mean/1000, '--k');
    hold on
    ErrorBarStep = ceil(simulation_time/60);
    errorbar(tq(1:ErrorBarStep:end)/consts.day2sec, Intersatellite_dist_mean(1:ErrorBarStep:end)/consts.km2m, Intersatellite_dist_std(1:ErrorBarStep:end)/consts.km2m);
    xlabel('time, days');
    ylabel('ISL, km');
    legend(['dV = ',num2str(abs(deployer.dV_dir(1)),2),' m/s; 3\sigma_{||} = ',num2str(deployer.dV_3sigma,2), ' m/s; 3\sigma_{\perp} = ', num2str(deployer.dV_3sigma_wide,2), ' m/s']);
    title('Swarm dispersion simulation with boundary conditions');
end

ISL_max = max(ISL_max)/1000;

%% saving to file
save(strcat('simulation_results\swarm_dispersion\BoundaryCond_dV_',num2str(abs(deployer.dV_dir(1)),2),'_3sigma_par_',num2str(deployer.dV_3sigma,2), '_3sigma_perp_', num2str(deployer.dV_3sigma_wide,2), '.mat'),...
    'spacecraft',...
    'deployer',...
    'LV',...
    'simulation_time',...
    'N_conditions',...
    'Target',...
    'Chaser',...
    'rv_ECI_deployment',...
    'Intersatellite_dist',...
    't_vec',...
    'Swarm_LT');

disp('Simulation is done. Data Saved Successfully.');