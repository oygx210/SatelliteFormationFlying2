% The script aims to find attitude of an image for demonstration & save
% orbital parameters for formation dynamics and control study

target_orbit = table2array(Trade_off_table(targ_orb,2:7)); % 
target_orbit(1) = consts.rEarth + target_orbit(1)*1000;    % [m] sma
target_orbit(2) = target_orbit(2);                         % [-] ecc
target_orbit(3) = target_orbit(3)*pi/180;                  % [rad] inc
target_orbit(4) = target_orbit(4)*pi/180;                  % [rad] RAAN
target_orbit(5) = target_orbit(5)*pi/180;                  % [rad] AOP
target_orbit(6) = target_orbit(6)*pi/180;                  % [rad] Mean anomaly

orbit_epoch = table2array(Trade_off_table(targ_orb,8));
mean_motion = sqrt(consts.muEarth/target_orbit(1)^3);
ISD_min = table2array(Trade_off_table(targ_orb,10));
alpa_demo1_demonstration = 0;
alpa_demo2_demonstration = 0;

reconf_1_time = POI.t(1);
alpha_reconf1 = mod(alpa_demo1_demonstration - mean_motion*(seconds(morning_demo(targ_orb) - reconf_1_time)), 2*pi);

reconf2_time = POI.t(morning_demo_stop_step(targ_orb));
alpha_reconf2 = mod(alpa_demo2_demonstration - mean_motion*(seconds(evening_demo(targ_orb) - reconf2_time)), 2*pi);

MD.target_orbit = target_orbit;
MD.orbit_epoch = orbit_epoch;
MD.reconfiguration_time = [reconf_1_time;reconf2_time];
MD.image_attitude_at_reconf = [alpha_reconf1; alpha_reconf2];
MD.ISD_min = round(ISD_min);
save('E:\Spacecraft-formation-control\numerical_studies\Letters_impuslive_control\Mission_design\MD_output_Paris.mat','MD');