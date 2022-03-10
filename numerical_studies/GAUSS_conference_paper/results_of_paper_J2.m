Error = [0.1, 0.2];
L = [1e3; 5e3; 10e3];
Time_between_corrections = [2.1*3600; 3.5*3600];

Mission_duration = 60*86400;
mSat = 4; %kg
Isp = 280;

g0 = 9.8;

% ! Run firstly startup_formation_control -> Orbital_dynamics_J2.m -> the
% script

% modeling parameters
% Modeling was made using Orbital_dynamics_J2.m
% The parameters are the following r = 3000 m, alpha = pi/2;

for i = 1:length(Time_between_corrections)

    N_correction_maneuvers(i) = round(Mission_duration / Time_between_corrections(i));
    
end

%% Chaser relative orbit correction
dV.ten_percent = simple_control(consts, chaser.rv_HCW_ECI(:,Time_between_corrections(1) + 1), chaser.rv_ECI(:,Time_between_corrections(1) + 1));
dV.ten_percent_total = dV.ten_percent*N_correction_maneuvers(1);
mprop.ten_percent = mSat*(1-exp(-dV.ten_percent_total/Isp/g0));

dV.twenty_percent = simple_control(consts, chaser.rv_HCW_ECI(:,Time_between_corrections(2) + 1), chaser.rv_ECI(:,Time_between_corrections(2) + 1));
dV.twenty_total = dV.ten_percent*N_correction_maneuvers(2);
mprop.ten_percent = mSat*(1-exp(-dV.ten_percent_total/Isp/g0));


%% Deployment
dV.deployment = simple_control(consts, target.rv_ECI(:,1), chaser.rv_ECI(:,1));
mprop.ten_percent = mSat*(1-exp(-dV.deployment/Isp/g0));

