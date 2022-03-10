%% Readme for Swarm_dispersion_MC

% The script is used for experiments on the swarm dispersion during the mission
% The main goal is to define deployer parameters to ensure required intersatellite distance

% The experiment is made for two similar satellites corresponding to 3U CubeSat standard
% The force model takes into account J2 effect and atmospheric drag

% The CubeSat swarm is to be deployed from International Space Station 

% Experiment plan.
% 1. We randomly generate deployment velocity vector according to normal
% distribution function for various variances of deployment velocity error;
% 2. We propagate formation motion (formation consists of 2 sats) for 2
% months which is considered as mission duration according to current ConOps
% 3. During the modeling we keep tracking intersatellite distance and stop
% integrating when it exceeds the maximum ISL distance for the CubeSats Com system
% 4. Based on the simulations we define requirement for deployment velocity
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

deployer.dV_dir = [-1; 0; 0];         % [m/s]
deployer.dV_mean = 0;                   % [m/s]
deployer.dV_sigma = 0.025/3;             % [m/s]
deployer.dV_sigma_wide = 0.02/3;% [m/s]

%% Launch vehicle (LV) orbit - ISS orbit
LV.oe = zeros(6,1);
% OE = [sma[km], ecc[-], inc[rad], RAAN[rad], AOP[rad], MA[rad]]
LV.oe(1) = 400e3 + consts.rEarth;       % [km] sma
LV.oe(2) = 0;                           % ecc
LV.oe(3) = 52;                          % [rad] inc
LV.oe(4) = 0;                           % [rad] RAAN
LV.oe(5) = 0;                           % [rad] AOP
LV.oe(6) = 0;                           % [rad] Mean anomaly

LV.rv = oe2rv(consts, LV.oe);
LV.mean_motion = sqrt(consts.muEarth/LV.oe(1)^3);

%% Monte Carlo simulation

% Simulation parameters
simulation_time = 10*consts.day2sec; % s
stopevent = @(t, rv) connection_is_lost(t, rv, spacecraft);
options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',  stopevent);
tspan = 0:simulation_time;

N_samples = 10000; % number of samples

% Set of deployment conditions
Target.rv_orb_deployment = [zeros(1,N_samples); zeros(1,N_samples); zeros(1,N_samples); 
                            deployer.dep2orb * [deployer.dV_dir(1)*ones(1,N_samples) + normrnd(deployer.dV_mean*ones(1,N_samples), deployer.dV_sigma*ones(1,N_samples));
                            normrnd(zeros(1,N_samples), deployer.dV_sigma_wide*ones(1,N_samples));
                            normrnd(zeros(1,N_samples), deployer.dV_sigma_wide*ones(1,N_samples))]];                        

Target.rv_ECI_deployment = orb2ECI(LV, Target.rv_orb_deployment);
Target.dV_deployment_mag = sqrt(dot(Target.rv_orb_deployment(4:6,:),Target.rv_orb_deployment(4:6,:)));
Target.dV_deployment_theta = atand(sqrt(Target.rv_orb_deployment(5,:).^2 + Target.rv_orb_deployment(6,:).^2) ./ Target.rv_orb_deployment(4,:));
Target.dV_deployment_phi = atan2(Target.rv_orb_deployment(6,:),Target.rv_orb_deployment(5,:))*180/pi;

figure(1);
subplot(1,3,1);
histogram(Target.dV_deployment_mag);
xlabel('deployment velocity, [m/s]');
ylabel('number of samples');
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
        '3\sigma (dV_{||}) = ', num2str(deployer.dV_sigma, 2), ' [m/s]']);

subplot(1,3,2);
histogram(Target.dV_deployment_theta);
xlabel('deployment direction error, degrees');
ylabel('number of samples');
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
       '3\sigma (dV_{\perp}) = ', num2str(deployer.dV_sigma_wide, 2), ' [m/s]']);

title('Target satellite deployment conditions');
subplot(1,3,3);
polarplot(Target.dV_deployment_phi, Target.dV_deployment_theta, '+m');
ax = gca;
rruler = ax.RAxis;
rruler.Label.String = 'elevation';

Chaser.rv_orb_deployment = [zeros(1,N_samples); zeros(1,N_samples); zeros(1,N_samples); 
                            deployer.dep2orb * [deployer.dV_dir(1)*ones(1,N_samples) + normrnd(deployer.dV_mean*ones(1,N_samples), deployer.dV_sigma*ones(1,N_samples));
                            normrnd(zeros(1,N_samples), deployer.dV_sigma_wide*ones(1,N_samples));
                            normrnd(zeros(1,N_samples), deployer.dV_sigma_wide*ones(1,N_samples))]];                        

Chaser.rv_ECI_deployment = orb2ECI(LV, Chaser.rv_orb_deployment);
Chaser.dV_deployment_mag = sqrt(dot(Chaser.rv_orb_deployment(4:6,:),Chaser.rv_orb_deployment(4:6,:)));
Chaser.dV_deployment_theta = atand(sqrt(Chaser.rv_orb_deployment(5,:).^2 + Chaser.rv_orb_deployment(6,:).^2) ./ Chaser.rv_orb_deployment(4,:));
Chaser.dV_deployment_phi = atan2(Chaser.rv_orb_deployment(6,:),Chaser.rv_orb_deployment(5,:))*180/pi;

figure(2);
subplot(1,3,1);
histogram(Chaser.dV_deployment_mag);
xlabel('deployment velocity, [m/s]');
ylabel('number of samples');
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
        '3\sigma (dV_{||}) = ', num2str(deployer.dV_sigma, 2), ' [m/s]']);

subplot(1,3,2);
histogram(Chaser.dV_deployment_theta);
xlabel('deployment direction error, degrees');
ylabel('number of samples');
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
       '3\sigma (dV_{\perp}) = ', num2str(deployer.dV_sigma_wide, 2), ' [m/s]']);

title('Chaser satellite deployment conditions');
subplot(1,3,3);
polarplot(Chaser.dV_deployment_phi, Chaser.dV_deployment_theta, '+m');
ax = gca;
rruler = ax.RAxis;
rruler.Label.String = 'elevation';

rv_ECI_deployment = [Target.rv_ECI_deployment; Chaser.rv_ECI_deployment];

tic;
f = waitbar(0,'Please wait...');
for i = 1:N_samples

    sol{1,i} = ode45(@(t, rv) J2_atmo_2sats(t, rv, consts, spacecraft), tspan, rv_ECI_deployment(:,i), options);

    waitbar(i/N_samples,f,'Please wait...');
end
close(f);
toc;

%% saving to file
save(strcat('MonteCarlo_dV_',num2str(abs(deployer.dV_dir(1)),2),'_3sigma_par_',num2str(deployer.dV_sigma*3,2), '_3sigma_perp_', num2str(deployer.dV_sigma_wide*3,2), '.mat'),...
    'spacecraft',...
    'deployer',...
    'LV',...
    'simulation_time',...
    'N_samples',...
    'Target',...
    'Chaser',...
    'rv_ECI_deployment',...
    'sol');
disp('Simulation is done. Data is stored.');
