clear all;

%% Initializing
global environment
environment = 'J2';
consts = startup_formation_control();

load('simulation_u_optimized.mat');

data1 = load('schedule_54.mat');
data2 = load('schedule_100.mat');
data3 = load('schedule_200.mat');

cum_cost = 0;
for i = 1:size(data1.schedule,1)
    cum_cost = cum_cost + data1.schedule(i,5);
    cum_cost_54(i) = cum_cost;
end

cum_cost = 0;
for i = 1:size(data2.schedule,1)
    cum_cost = cum_cost + data2.schedule(i,5);
    cum_cost_100(i) = cum_cost;
end
cum_cost = 0;
for i = 1:size(data3.schedule,1)
    cum_cost = cum_cost + data3.schedule(i,5);
    cum_cost_200(i) = cum_cost;
end

figure();
plot(data1.schedule(:,3), cum_cost_54/1000);
hold on;
plot(data2.schedule(:,3), cum_cost_100/1000);
hold on;
plot(data3.schedule(:,3), cum_cost_200/1000);
hold on;
xlabel('time, seconds');
ylabel('Cumulative demonstration cost, kUSD');

% orbit epoch
orbit_epoch_GD = datetime(2022, 12, 22, 0, 0, 0);
orbit_epoch_JD = juliandate(orbit_epoch_GD);
k_rev2rep = 14;
k_day2rep = 1;
% orbital elements
oe = zeros(6,1);
oe(1) = get_SSO_RGT_orbit_sma(k_rev2rep, k_day2rep, consts);
oe(2) = 0; % ecc, 
oe(3) = get_SSO_inclination(oe(1), oe(2), consts);
oe(4) = get_RAAN_for_terminator_orbit(orbit_epoch_GD);
oe(5) = 0; % AOP, deg
oe(6) = u_optimal; % M, deg - actually argument of latitude
% oe(6) = deg2rad(140); % M, deg - actually argument of latitude

spacecraft = [];

% Importing database with cities parameters
cities_database = readtable('Table_cities');

% Map resolution parameters
earth = imread('earthmap1k.jpg');
scale1 = size(earth,1);
scale2 = size(earth,2);

% T_simulation = round(2*pi*sqrt(oe(1)^3/consts.muEarth))*2;                % seconds
T_simulation = consts.day2sec;                % seconds
t_step = 30;

% demonstration parameters
theta_sat_min = deg2rad(10); % deg
theta_sun_max = deg2rad(-5); % deg
demo_duration = 60; % s
foot_print_area = 54;

% % prepating images for demonstration
% ISD_angular = deg2rad(1/60); % min angular distance between formation satellite defined by human eye resolution
% d_max = roots([1, -2*consts.rEarth*cos(pi/2 + theta_sat_min),consts.rEarth^2 - oe(1)^2]);
% d_max = max(d_max);
% ISD_min = 2*tan(ISD_angular/2)*d_max;
% 
% % XY_HCW_constants = Image2HCW_constants("XY", "PCO", ISD_min, 0);
% rho_phi_HCW_constants = Image2HCW_constants("rho_phi", "PCO", ISD_min, 0);
% 

%% Preparing cost_matrix and area_matrix
cities_database_matrix = table2array(cities_database(:,2:end));
cities_database_matrix(:,1:2) = deg2rad(cities_database_matrix(:,1:2));
cost_matrix = zeros(scale1,scale2);
area_matrix = zeros(scale1,scale2);
[M, N] = size(cost_matrix);
phi_range = linspace(-pi, pi, N+1);
phi_range = phi_range(1:end-1);
th_range = linspace(-pi/2, pi/2, M+1);
th_range = th_range(1:end-1);

for i = 1:size(cities_database_matrix,1)
    [~,ind1] = min(abs(cities_database_matrix(i,1) - th_range));
    [~,ind2] = min(abs(cities_database_matrix(i,2) - phi_range));
    cost_matrix(ind1, ind2) = cities_database_matrix(i,5);
    area_matrix(ind1, ind2) = cities_database_matrix(i,3);
    
    if area_matrix(ind1, ind2) > foot_print_area
        cost_matrix(ind1, ind2) = cost_matrix(ind1, ind2) * foot_print_area / area_matrix(ind1, ind2);
%         cost_matrix(ind1, ind2) = cost_matrix(ind1, ind2);
    end
        
end

%% Orbit propagation and SV convertion to ECEF coordinate system
rv_init = oe2rv(oe, consts);
dt = [0:t_step:T_simulation];
options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec, rv_ECI] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), dt, rv_init, options_precision);
rv_ECI = rv_ECI';
t_vec_JD = orbit_epoch_JD + t_vec/consts.day2sec;
r_ECEF = ECItoECEF(t_vec_JD,rv_ECI(1:3,:));


plot_ground_track_with_cost_map(r_ECEF, t_vec_JD, cost_matrix, theta_sat_min, theta_sun_max, consts)

% Cost analysis
tic;
schedule = calculate_cities_coverage_fixed_demo_time(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts);
toc;
Total_cost = sum(schedule(:,end));

% load("schedule_1day_optimal_GT.mat");
% coverage_animation(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, cost_matrix, schedule, consts);

