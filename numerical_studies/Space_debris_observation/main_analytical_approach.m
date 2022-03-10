clear all;

%% Initialization
global environment 
environment = 'J2';
syms t

consts = startup_formation_control();

spacecraft.dr_observation_max = 500e3; 

orbit_epoch = datetime(2022, 1, 1, 0, 0, 0);
misson_duration = seconds(days(1));

% Observer orbit - circular SSO, 700km, terminator
oe_o(1) = consts.rEarth + 700e3;
oe_o(2) = 0;
oe_o(3) = get_SSO_inclination(oe_o(1), oe_o(2), consts);
oe_o(4) = get_RAAN_for_terminator_orbit(orbit_epoch);
oe_o(5:6) = 0;

rv_obs_fun = analytical_rv_function(oe_o, 'm, s', consts);

% Debris orbit - circular SSO, 900km, perp to terminator

oe_d(1) = consts.rEarth + 750e3;
oe_d(2) = 0;
% oe_d(3) = get_SSO_inclination(oe_o(1), oe_o(2), consts);
oe_d(3) = 0;
oe_d(4) = get_RAAN_for_terminator_orbit(orbit_epoch) + pi/2;
oe_d(5:6) = 0;

rv_debr_fun = analytical_rv_function(oe_d, 'm, s', consts);
dr = norm(rv_obs_fun(1:3) - rv_debr_fun(1:3));
equation = dr - sym(spacecraft.dr_observation_max) == sym(0);
tau_analytical = custom_solver(equation, t, 1e4);

f = dr - spacecraft.dr_observation_max;
fplot(f, [1 10*consts.day2sec]);

% a_d_range = consts.rEarth + [600e3 800e3];
% i_d_range = [96*pi/180 98*pi/180];
% RAAN_d_range = [oe_o(4) - 0.005 oe_o(4) + 0.005];
% M_d_range = [0 2*pi];
% variables_range = [a_d_range; i_d_range; RAAN_d_range; M_d_range];
% variables0 = [consts.rEarth + 700e3; 97*pi/180; oe_o(4); 0];
% 
% % tau = time_to_observation(variables, variables_range, rv_obs_fun, dr_observation_max, consts);
% fun = @(variables)time_to_observation(variables, variables_range, rv_target_fun, spacecraft.dr_observation_max, consts);
% tau_max = fminsearch(fun,variables0);
