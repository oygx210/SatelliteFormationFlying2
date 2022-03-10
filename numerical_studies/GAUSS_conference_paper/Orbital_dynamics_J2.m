%% All vectors are defined as vector columns

%% target satellite initial orbit

% Orbital elements
target.oe_in = zeros(6,1);
target.oe_in(1) = 350e3 + consts.rEarth;                                          % sma, km
target.oe_in(2) = 0;                                                              % ecc
target.oe_in(3) = get_SSO_inclination(consts, target.oe_in(1), target.oe_in(2));  % inc, rad
target.oe_in(4) = 0;                                                              % RAAN, rad
target.oe_in(5) = 0;                                                              % AOP, rad
target.oe_in(6) = 0;                                                              % Mean anomaly, rad

% State vector
target.rv_in = oe2rv(consts, target.oe_in);

mean_motion = sqrt(consts.muEarth/target.oe_in(1)^3);
orbit_period = 2 * pi / mean_motion;

target.rv_in = oe2rv(consts, target.oe_in);

%% Chaser satellites initial conditions

% Setting constants to define relative trajectory of chaser SC
alpha = pi/2;
r = 0;
c1 = r;
c2 = sqrt(3)/2*r;
c3 = 30000;

chaser.rv_orb_in = get_rv_from_analytic_HCW_solution(c1, c2, c3, alpha, target.oe_in(5) + target.oe_in(6), mean_motion);
chaser.rv_ECI_in = orb2ECI(target.rv_in, chaser.rv_orb_in, mean_motion);
   
%% Integration

% Conditions for modelling
dt = 1; % s
% simulation_time = orbit_period*30; % s
simulation_time = 3*86400; % s


options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
tspan = [0:dt:simulation_time];

% Target satellite, ECI
target_init = target.rv_in;
[target.dt_vec, target.rv_ECI] = ode45(@(t, rv)J2(t, rv, consts), tspan, target_init, options);
target.rv_ECI = target.rv_ECI';

% Chaser satellite, ECI
chaser_init = chaser.rv_ECI_in;
[chaser.dt_vec, chaser.rv_ECI] = ode45(@(t, rv)J2(t, rv, consts), tspan, chaser_init, options);
chaser.rv_ECI = chaser.rv_ECI';

% Chaser satellite, Orb, linearized dynamics
chaser_init_orb = chaser.rv_orb_in;
[chaser.dt_vec, chaser.rv_orb_HCW] = ode45(@(t, rv_orb)hcw(t, rv_orb, mean_motion), tspan, chaser_init_orb, options);
chaser.rv_orb_HCW = chaser.rv_orb_HCW';
chaser.r_norm_HCW_ODE = sqrt(chaser.rv_orb_HCW(1,:).^2 + chaser.rv_orb_HCW(2,:).^2 + chaser.rv_orb_HCW(3,:).^2);

for i = 1:length(chaser.dt_vec)
    
    chaser.rv_HCW_ECI(:,i) = orb2ECI(target.rv_ECI(:,i) ,chaser.rv_orb_HCW(:,i), mean_motion);
    
end


% Converts Chaser satellite position from inertial to orbital reference
% frame
for i = 1:length(chaser.dt_vec)
    chaser.rv_orb(:,i) = ECI2orb(target.rv_ECI(:,i), chaser.rv_ECI(:,i), mean_motion);
    
    chaser.r_orb_norm(i) = norm(chaser.rv_orb(:,i));
end

%% Results demonstration

figure(1);
plot_Earth_meters();
hold on;
plot3(chaser.rv_ECI(1,:), chaser.rv_ECI(2,:), chaser.rv_ECI(3,:));
hold on;
plot3(target.rv_ECI(1,:), target.rv_ECI(2,:), target.rv_ECI(3,:));
% title('Target and Chaser satellite trajectories in ECI');
xlabel('x, m');
ylabel('y, m');
zlabel('z, m');
grid on;
axis equal;

figure(2);
% subplot(1,2,1);
% plot3(chaser.rv_ECI(1,:) - target.rv_ECI(1,:), chaser.rv_ECI(2,:) - target.rv_ECI(2,:),...
%         chaser.rv_ECI(3,:) - target.rv_ECI(3,:));
% hold on;
% xlabel('x, m');
% ylabel('y, m');
% zlabel('z, m');
% title('(R_{chaser} - R_{target}) inertial rf');

subplot(1,2,1);
plot3(chaser.rv_orb(1,:),chaser.rv_orb(2,:),chaser.rv_orb(3,:), 'k', 0,0,0, 'ok');
hold on;
xlabel('x, m');
ylabel('y, m');
zlabel('z, m');
% title('Chaser position in orbital rf');
% axis equal;
% 
subplot(1,2,2);
plot(chaser.dt_vec./3600, chaser.r_orb_norm, 'k');
hold on;
xlabel('t, hours');
ylabel('dist, m');
% title('norm(R_{chaser} - R_{target})');