clear all;

%% Initializing
global environment
environment = 'J2';
consts = startup_formation_control();

spacecraft = [];

% Importing database with cities parameters
cities_database = readtable('Table_cities');

% Map resolution parameters
earth = imread('earthmap1k.jpg');
scale1 = size(earth,1);
scale2 = size(earth,2);

% demonstration parameters
theta_sat_min = deg2rad(10); % deg
theta_sun_max = deg2rad(-5); % deg
demo_duration = 60; % s
foot_print_area = 200;

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

options = optimset('PlotFcns',@optimplotfval);
options.Tolx = 1e-5;
options.TolFun = 1e3;
% [t_JD_optimal, C_optimal] = fminbnd(@(orbit_epoch_JD)find_total_demonstration_cost(orbit_epoch_JD, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts),...
%                         orbit_epoch_JD_min,orbit_epoch_JD_max, options);

% u_range = linspace(0, 2*pi - 5*pi/180, round((2*pi - 5*pi/180) / (5*pi/180)));
% 
% tic;
% for i = 1:length(u_range)
%     C(i) = find_total_demonstration_cost(u_range(i), theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts);
% end
% toc;
% 
% figure();
% plot(u_range, C);
% 
u0 = 0.619592;
u_range = [0, (2*pi -1e-3)];

tic;
[u_optimal, C_optimal, exitflag, output] = fminsearch(@(u)find_total_demonstration_cost(u, u_range, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts),...
                        u0, options);
toc;
function C = find_total_demonstration_cost(u, u_range, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts)

    if u < u_range(1) || u > u_range(2)
        C = 0;
    else
        % SSO RGT orbit parameters
        k_rev2rep = 14;
        k_day2rep = 1;
        % orbital elements
        orbit_epoch_GD = datetime(2022, 12, 22, 0, 0, 0);
        orbit_epoch_JD = juliandate(orbit_epoch_GD);

        oe = zeros(6,1);
        oe(1) = get_SSO_RGT_orbit_sma(k_rev2rep, k_day2rep, consts);
        oe(2) = 0; % ecc, 
        oe(3) = get_SSO_inclination(oe(1), oe(2), consts);
        orbit_epoch_GD = datetime(orbit_epoch_JD, 'convertfrom','juliandate');
        oe(4) = get_RAAN_for_terminator_orbit(orbit_epoch_GD);
        oe(5) = 0; % AOP, deg
        oe(6) = u; % M, deg - actually argument of latitude
        disp(u);
%         T_simulation = round(2*pi*sqrt(oe(1)^3/consts.muEarth));                % seconds
        T_simulation = consts.day2sec;                % seconds
        t_step = 30;
        spacecraft = [];
        %% Orbit propagation and SV convertion to ECEF coordinate system
        rv_init = oe2rv(oe, consts);
        dt = [0:t_step:T_simulation];
        options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [t_vec, rv_ECI] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), dt, rv_init, options_precision);
        rv_ECI = rv_ECI';
        t_vec_JD = orbit_epoch_JD + t_vec/consts.day2sec;
        r_ECEF = ECItoECEF(t_vec_JD,rv_ECI(1:3,:));
                
        % Demonstration cost calculation
        % tic;
        schedule = calculate_cities_coverage_fixed_demo_time(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts);
        % toc;
        C = -sum(schedule(:,end));
    end
    disp(C);
end