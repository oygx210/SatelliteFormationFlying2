% Near-circular orbit;
% Launch from ISS via P-POD with same parameters as GAUSS P-POD has

clear all;
consts = startup_formation_control();

% Spacecraft parameters
spacecraft.DragArea = 0.0481;           % [m^2] cross-section area perpendicular to the incoming airflow
spacecraft.Cdrag = 2.2;                 % atmospheric drag coefficient
spacecraft.mass = 4;                    % [kg]

% Launch vehicle (LV) orbit [sma[m], ecc[-], inc[rad], RAAN[rad], AOP[rad], MA[rad]]
LV.coe(1) = 400e3 + consts.rEarth;      % [m] sma
LV.coe(2) = 4e-4;                       % [-] ecc
LV.coe(3) = 52*consts.deg2rad;          %[rad] inc
LV.coe(4) = 0;                          % [rad] RAAN
LV.coe(5) = 0;                          % [rad] AOP
LV.coe(6) = 0;                          % [rad] Mean anomaly
LV.rv = oe2rv(LV.coe, consts);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec, rv_ECI] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:2*pi*sqrt(LV.coe(1)^3/consts.muEarth)*10, LV.rv, options);
rv_ECI = rv_ECI';

coe = zeros(8, size(rv_ECI, 2));

for i = 1:size(rv_ECI, 2)
    coe(:,i) = rv2coe(rv_ECI(:,i), consts); 
    mean_coe(:,i) = osc2mean(coe(:,i), consts);
end

time_scale = 'h';
if time_scale == 's'
    time_res = 1;
elseif time_scale == 'm'
    time_res = 60;
elseif time_scale == 'h'
    time_res = 60*60;    
elseif time_scale == 'd'
    time_res = 60*60*24;
end


figure;
subplot(2,4,1);
plot(t_vec/time_res, coe(1,:), t_vec/time_res, mean_coe(1,:));
legend('Osculating oe', 'Mean oe');
ylabel('sma, m');
xlabel(['time, ', time_scale]);

subplot(2,4,2);
plot(t_vec/time_res, coe(2,:), t_vec/time_res, mean_coe(2,:));
ylabel('ecc, -');
xlabel(['time, ', time_scale]);

title('Orbital elements');

subplot(2,4,3);
plot(t_vec/time_res, coe(3,:), t_vec/time_res, mean_coe(3,:));
ylabel('inc, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,4);
plot(t_vec/time_res, coe(4,:), t_vec/time_res, mean_coe(4,:));
ylabel('RAAN, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,5);
plot(t_vec/time_res, coe(5,:), t_vec/time_res, mean_coe(5,:));
ylabel('AOP, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,6);
plot(t_vec/time_res, coe(6,:), t_vec/time_res, mean_coe(6,:));
ylabel('\nu, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,7);
plot(t_vec/time_res, coe(7,:), t_vec/time_res, mean_coe(7,:));
ylabel('M, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,8);
plot(t_vec/time_res, coe(8,:), t_vec/time_res, mean_coe(8,:));
ylabel('u, rad');
xlabel(['time, ', time_scale]);

