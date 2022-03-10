% The script aims to find the target orbit for demonstration based on Mission_design_step1

%% Input parameters
POI_latitude = deg2rad(48.8566); % degrees
min_sat_elevation_for_demonstration = 10;

h_min = 350e3; % min admissible orbit altitude
h_max = 800e3; % max admissible orbit altitude
N = 1:40; % range of possible number of revolutions

I0 = 1360;                  % The average intensity of solar energy at the earth distance (W/m^2)
Iref = (2.54e-6)/683;       % Intensity of the light from the brightness star in W/m^2 
alpha = 9.308422677e-3;     % angular size of the Sun from the Earth
rho = 0.92;                 % The reflectivity coefficient for the thin Mylar film coated with aluminum 
m_required = -8; % required magnitude for demonstration
I_req = Iref*10^(-m_required/2.5);

ISD_angular = deg2rad(1/60); % min angular distance between formation satellite defined by human eye resolution

%% Mission design procedure, Step 1

% The case when POI is located at northern hemisphere
R_Sun_z_sign = sign(round(Sun.r(3,1)));

if R_Sun_z_sign > 0 
    
    TA1 = POI_latitude;
    TA2 = pi - POI_latitude;
    dTA = TA2 - TA1;

    h_orbitX = -Sun.r(1,1);
    h_orbitY = -Sun.r(2,1);
    RAAN = atan2(h_orbitX, -h_orbitY);
    
    if RAAN < 0
        RAAN = 2*pi + RAAN;
    end
else
    
    TA1 = pi - POI_latitude;
    TA2 = POI_latitude;
    dTA = 2*pi - TA1 + TA2;
    
    if dTA < 0 
        dTA = 2*pi + dTA;
    end

    h_orbitX = Sun.r(1,1);
    h_orbitY = Sun.r(2,1); 
    RAAN = atan2(h_orbitX, -h_orbitY);
    
    if RAAN < 0
        RAAN = 2*pi + RAAN;
    end

end

%% Mission design procedure, Step 2
dT = time_between_demonstrations;
Time_initial = POI.t(1);

T_orb = dT'./(N + dTA/(2*pi));
h = (T_orb.^2*consts.muEarth/4/pi^2).^(1/3) - consts.rEarth;

logical = ~(h > h_min & h < h_max);
h(logical) = NaN;
logical = ~(h <= max(h,[],2));
h(logical) = NaN;
logical = isnan(h);
h(logical) = [];
h = h';

mean_motion = sqrt(consts.muEarth./(h + consts.rEarth).^3);
incl = get_SSO_inclination(h + consts.rEarth, zeros(length(h),1), consts);

%% Final data processing
TA_initial = mod(TA1*ones(length(h),1) - mean_motion.*seconds(morning_demo - Time_initial)',2*pi);

for i = 1:length(sun_elevation_at_demo)
   Orbit_epoch(i) =  POI.t(1);
end

Trade_off_table = table;

Trade_off_matrix = [sun_elevation_at_demo', h, zeros(length(sun_elevation_at_demo),1),...
                   rad2deg(incl), ones(length(sun_elevation_at_demo),1)*rad2deg(RAAN), zeros(length(sun_elevation_at_demo),1), rad2deg(TA_initial)];

Trade_off_table.Variables = Trade_off_matrix;
Trade_off_table.Properties.VariableNames = {'Sun elevation at demo','orbit altitude','ecc','incl', 'RAAN', 'AOP', 'TA'};
Trade_off_table.epoch = Orbit_epoch';

oe_initial = Trade_off_matrix(:,2:7);
oe_initial(:,1) = oe_initial(:,1) + consts.rEarth;
oe_initial(:,3) = deg2rad(oe_initial(:,3));
oe_initial(:,4) = deg2rad(oe_initial(:,4));
oe_initial(:,6) = deg2rad(oe_initial(:,6));

oe_initial = oe_initial';

%% orbital dynamics simulation
for i = 1:size(oe_initial,2)
    rv_initial(:,i) = oe2rv(oe_initial(:,i), consts);
end

spacecraft.DragArea = 340e-3*200e-3;                            % [m^2] cross-section area perpendicular to the incoming airflow
spacecraft.Cdrag = 2.2;                                         % atmospheric drag coefficient
spacecraft.mass = 16;                                           % [kg]

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec, rv_ECI] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:consts.day2sec, rv_initial, options_precision);
rv_ECI = rv_ECI';

% figure(1);
% for j = 1:size(Trade_off_table,1)
%         plot3(rv_ECI(6*j-5,1:8000), rv_ECI(6*j-4,1:8000), rv_ECI(6*j-3,1:8000));
%         hold on;
% end
% earth_sphere('m');

%% Raw demonstration parameters calculation

for i = 1:size(Trade_off_table,1)
    r_morning(:,:,i) = rv_ECI(6*i-5:6*i-3,morning_demo_step(i) - 600:morning_demo_step(i) + 600);
    rPOI_morning(:,:,i) = POI.r(:, morning_demo_step(i) - 600:morning_demo_step(i) + 600)*1000;
    rPOI_norm_morning(i,:) = vecnorm(rPOI_morning(:,:,i));

    rPOI2Sat_morning(:,:,i) = r_morning(:,:,i) - rPOI_morning(:,:,i);   % vector from POI to SAT
    rPOI2Sat_norm_morning(i,:) = vecnorm(rPOI2Sat_morning(:,:,i));
    rSun_morning(:,:,i) = Sun.r(:, morning_demo_step(i) - 600:morning_demo_step(i) + 600)*1000;
    rSun_norm_morning(i,:) = vecnorm(rSun_morning(:,:,i));
    
    rSat2POI_morning(:,:,i) = rPOI_morning(:,:,i) - r_morning(:,:,i); % vectors from Sat to POI
    rSat2POI_norm_morning(i,:) = vecnorm(rSat2POI_morning(:,:,i));
    
    rSat2Sun_morning(:,:,i)  = rSun_morning(:,:,i) - r_morning(:,:,i);
    rSat2Sun_norm_morning(i,:) = vecnorm(rSat2Sun_morning(:,:,i));
    
    Sat_elevation_morning(i,:) = rad2deg(pi/2 - acos(dot(rPOI2Sat_morning(1:3,:,i),rPOI_morning(1:3,:,i))./rPOI2Sat_norm_morning(i,:)./rPOI_norm_morning(i,:)));
    min_Sat_elevation_morning(i,:) = min(Sat_elevation_morning(i,:));
    max_Sat_elevation_morning(i,:) = max(Sat_elevation_morning(i,:));
    
    incident_angle_morning(i,:) = rad2deg(acos(dot(rSat2POI_morning(1:3,:,i),rSat2Sun_morning(1:3,:,i))./rSat2POI_norm_morning(i,:)./rSat2Sun_norm_morning(i,:))/2);
    min_incident_angle_morning(i,:) = min(incident_angle_morning(i,:));
    max_incident_angle_morning(i,:) = max(incident_angle_morning(i,:));
    
    r_evening(:,:,i) = rv_ECI(6*i-5:6*i-3,evening_demo_step(i) - 600:evening_demo_step(i) + 600);
    rPOI_evening(:,:,i) = POI.r(:, evening_demo_step(i) - 600:evening_demo_step(i) + 600)*1000;
    rPOI_norm_evening(i,:) = vecnorm(rPOI_evening(:,:,i));
    rPOI2Sat_evening(:,:,i) = r_evening(:,:,i) - rPOI_evening(:,:,i);   % vector from POI to SAT
    rPOI2Sat_norm_evening(i,:) = vecnorm(rPOI2Sat_evening(:,:,i));
    rSun_evening(:,:,i) = Sun.r(:, evening_demo_step(i) - 600:evening_demo_step(i) + 600)*1000;
    rSun_norm_evening(i,:) = vecnorm(rSun_evening(:,:,i));

    rSat2POI_evening(:,:,i) = rPOI_evening(:,:,i) - r_evening(:,:,i); % vectors from Sat to POI
    rSat2POI_norm_evening(i,:) = vecnorm(rSat2POI_evening(:,:,i));

    rSat2Sun_evening(:,:,i)  = rSun_evening(:,:,i) - r_evening(:,:,i);
    rSat2Sun_norm_evening(i,:) = vecnorm(rSat2Sun_evening(:,:,i));

    incident_angle_evening(i,:) = rad2deg(acos(dot(rSat2POI_evening(1:3,:,i),rSat2Sun_evening(1:3,:,i))./rSat2POI_norm_evening(i,:)./rSat2Sun_norm_evening(i,:))/2);
    min_incident_angle_evening(i,:) = min(incident_angle_evening(i,:));
    max_incident_angle_evening(i,:) = max(incident_angle_evening(i,:));

    Sat_elevation_evening(i,:) = rad2deg(pi/2 - acos(dot(rPOI2Sat_evening(1:3,:,i),rPOI_evening(1:3,:,i))./rPOI2Sat_norm_evening(i,:)./rPOI_norm_evening(i,:)));
    min_Sat_elevation_evening(i,:) = min(Sat_elevation_evening(i,:));
    max_Sat_elevation_evening(i,:) = max(Sat_elevation_evening(i,:));
    
end

for i = 1:size(Trade_off_table,1)
    
   for j = 1:(size(Sat_elevation_morning,2)-1)  
       if Sat_elevation_morning(i,j+1) > min_sat_elevation_for_demonstration && Sat_elevation_morning(i,j) < min_sat_elevation_for_demonstration
            morning_demo_start_step(i)  = morning_demo_step(i) - 600 + j;
            morning_demo_start_step_internal(i) = j;
       end
       if Sat_elevation_morning(i,j+1) < min_sat_elevation_for_demonstration && Sat_elevation_morning(i,j) > min_sat_elevation_for_demonstration
            morning_demo_stop_step(i)  = morning_demo_step(i) - 600 + j;
            morning_demo_stop_step_internal(i) = j;
       end
    
       if Sat_elevation_evening(i,j+1) > min_sat_elevation_for_demonstration && Sat_elevation_evening(i,j) < min_sat_elevation_for_demonstration
            evening_demo_start_step(i)  = evening_demo_step(i) - 600 + j;
            evening_demo_start_step_internal(i) = j;
       end
       if Sat_elevation_evening(i,j+1) < min_sat_elevation_for_demonstration && Sat_elevation_evening(i,j) > min_sat_elevation_for_demonstration
            evening_demo_stop_step(i)  = evening_demo_step(i) - 600 + j;
            evening_demo_stop_step_internal(i) = j;
       end
   end   
   
   if i > length(morning_demo_start_step) || i > length(evening_demo_start_step) ...
      || isnan(morning_demo_start_step(i)) || isnan(evening_demo_start_step(i))
        morning_demo_start_step(i) = NaN;
        morning_demo_stop_step(i) = NaN;
        evening_demo_start_step(i) = NaN;
        evening_demo_stop_step(i) = NaN;
        morning_demo_duration(i) = NaN;
        evening_demo_duration(i) = NaN;
        mean_demo_duration(i) = NaN;
        max_d(i) = NaN;
        min_ISD(i) = NaN;
        Ar(i) = NaN;
        side(i) = NaN;
   else
        morning_demo_duration(i) = morning_demo_stop_step(i) - morning_demo_start_step(i);
        evening_demo_duration(i) = evening_demo_stop_step(i) - evening_demo_start_step(i);
        mean_demo_duration(i) = (morning_demo_duration(i) + evening_demo_duration(i))/2;

        d_morning = vecnorm(rv_ECI(6*i-5:6*i-3,morning_demo_start_step(i):morning_demo_stop_step(i)) - POI.r(:,morning_demo_start_step(i):morning_demo_stop_step(i))*1000);
        d_evening = vecnorm(rv_ECI(6*i-5:6*i-3,evening_demo_start_step(i):evening_demo_stop_step(i)) - POI.r(:,evening_demo_start_step(i):evening_demo_stop_step(i))*1000);
        max_d = max([d_morning d_evening]);
        min_ISD(i) = 2*tan(ISD_angular/2)*max_d;

        gamma_morning = incident_angle_morning(i,morning_demo_start_step_internal(i):morning_demo_stop_step_internal(i));
        gamma_evening = incident_angle_morning(i,evening_demo_start_step_internal(i):evening_demo_stop_step_internal(i));

        theta_sat_morning = Sat_elevation_morning(i,morning_demo_start_step_internal(i):morning_demo_stop_step_internal(i));
        theta_sat_evening = Sat_elevation_evening(i,evening_demo_start_step_internal(i):evening_demo_stop_step_internal(i));

        tau_morning = atmospheric_transmissivity(deg2rad(theta_sat_morning));
        tau_evening = atmospheric_transmissivity(deg2rad(theta_sat_evening));
        
        I_sp_morning = tau_morning.*cosd(gamma_morning).*sind(theta_sat_morning)./d_morning.^2;
        I_sp_evening = tau_evening.*cosd(gamma_evening).*sind(theta_sat_evening)./d_evening.^2;
       
        % check this for different cases
        figure(2)
        plot(I_sp_morning);
        hold on;
        plot(I_sp_evening, '--');
        title('specific pixel intenstity during morning and evening demonstrations');
        legend('morning demo', 'evening demo');
        ylabel('I_{specific}, W/m^2');
        
        [I_sp_morning_min, index_morning] = min(I_sp_morning);
        [I_sp_evening_min, index_evening] = min(I_sp_evening);
        
        if I_sp_morning_min < I_sp_evening_min            
            I_sp_min = I_sp_morning_min;
        else
            I_sp_min = I_sp_evening_min;
        end
        
        Ar(i) = I_req*4*tan(alpha/2)^2/(I0*rho*I_sp_min);
        side(i) = sqrt(Ar(i));
        
   end
   
   if i == targ_orb
        I_sp_morning_target = I_sp_morning;
        I_sp_evening_target = I_sp_evening;
   end
    
end

Trade_off_table.demo_duration = (mean_demo_duration/60)'; 
Trade_off_table.min_ISD = min_ISD';
Trade_off_table.SS_side_6 = side';
    
Trade_off_table.("orbit altitude") = Trade_off_table.("orbit altitude")./1000;

% writetable(Trade_off_table,'E:\Spacecraft-formation-control\numerical_studies\Letters_impuslive_control\Mission_designs','Delimiter',',');

%% Calculating and plotting demonstration parameters
% theta_sat, incident angle, magnitude, d

targ_orb = 5; % index of the target orbti in the table

theta_sat = Sat_elevation_morning(targ_orb, morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb));
tau = atmospheric_transmissivity(deg2rad(theta_sat));
gamma = incident_angle_morning(targ_orb,morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb));
theta_sat = Sat_elevation_morning(targ_orb,morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb));
I_reflector_morning = I0*Ar(targ_orb)*rho*tau.*cosd(gamma).*sind(theta_sat)./(4*rSat2POI_norm_morning(targ_orb,morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb)).^2*tan(alpha/2)^2);
m_reflector_morning = -2.5*log10(I_reflector_morning/Iref);

theta_sat = Sat_elevation_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb));
tau = atmospheric_transmissivity(deg2rad(theta_sat));
gamma = incident_angle_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb));
theta_sat = Sat_elevation_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb));
I_reflector_evening = I0*Ar(targ_orb)*rho*tau.*cosd(gamma).*sind(theta_sat)./(4*rSat2POI_norm_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb)).^2*tan(alpha/2)^2);
m_reflector_evening = -2.5*log10(I_reflector_evening/Iref);

for i = 1:(morning_demo_stop_step(targ_orb) - morning_demo_start_step(targ_orb))+1
    t_vec_morning_demo(i) = POI.t(1) + seconds(morning_demo_start_step(targ_orb)) + seconds(i) - seconds(1);
end

for i = 1:(evening_demo_stop_step(targ_orb) - evening_demo_start_step(targ_orb))+1
    t_vec_evening_demo(i) = POI.t(1) + seconds(evening_demo_start_step(targ_orb)) + seconds(i) - seconds(1);
end

% % check that demonstration parameters are met
% figure(3);
% subplot(1,3,1);
% plot(t_vec_morning_demo, Sat_elevation_morning(targ_orb, morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time');
% ylabel('\theta_{sat}, deg');
% 
% subplot(1,3,2);
% plot(t_vec_morning_demo, incident_angle_morning(targ_orb, morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time');
% ylabel('\gamma, deg');
% 
% subplot(1,3,3);
% plot(t_vec_morning_demo, m_reflector_morning);
% grid on;
% hold on;
% xlabel('time');
% ylabel('m');
% legend('m_{max} = -8', 'm_{max} = -6');
% 
% figure(4);
% subplot(1,3,1);
% plot(t_vec_evening_demo, Sat_elevation_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time');
% ylabel('\theta_{sat}, deg');
% 
% subplot(1,3,2);
% plot(t_vec_evening_demo, incident_angle_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time');
% ylabel('\gamma, deg');
% 
% subplot(1,3,3);
% plot(t_vec_evening_demo, m_reflector_evening);
% grid on;
% hold on;
% xlabel('time');
% ylabel('m');
% legend('m_{max} = -8', 'm_{max} = -6');
% 
% figure(5);
% subplot(2,2,1);
% plot(t_vec_morning_demo, Sun_elevation(morning_demo_start_step(targ_orb):morning_demo_stop_step(targ_orb)));
% grid on;
% hold on;
% xlabel('time (UTC+2)');
% ylabel('\theta_{Sun}, deg');
% legend('morning demonstration');
% subplot(2,2,2);
% plot(t_vec_morning_demo, Sat_elevation_morning(targ_orb, morning_demo_start_step_internal(targ_orb):morning_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time (UTC+2)');
% ylabel('\theta_{sat}, deg');
% legend('morning demonstration');
% 
% subplot(2,2,3);
% plot(t_vec_evening_demo, Sun_elevation(evening_demo_start_step(targ_orb):evening_demo_stop_step(targ_orb)));
% grid on;
% hold on;
% xlabel('time (UTC+2)');
% ylabel('\theta_{Sun}, deg');
% legend('evening demonstration');
% subplot(2,2,4);
% plot(t_vec_evening_demo, Sat_elevation_evening(targ_orb, evening_demo_start_step_internal(targ_orb):evening_demo_stop_step_internal(targ_orb)));
% grid on;
% hold on;
% xlabel('time (UTC+2)');
% ylabel('\theta_{sat}, deg');
% legend('evening demonstration');

figure(6);
subplot(1,3,1);
plot(I_sp_morning_target);
hold on;
plot(I_sp_evening_target);
hold off;
xlabel('time, seconds');
ylabel('I_{sp}');
legend('morning demonstration','evening demonstration');
grid on;

subplot(1,3,2);
plot(t_vec_morning_demo, m_reflector_morning);
grid on;
hold on;
legend('A_r = 21.5 m^2', 'A_r = 135.5 m^2');
xlabel('time (UTC+2)');
ylabel('m');

subplot(1,3,3);
plot(t_vec_evening_demo, m_reflector_evening);
grid on;
hold on;
xlabel('time');
ylabel('m');
legend('A_r = 21.5 m^2', 'A_r = 135.5 m^2');

