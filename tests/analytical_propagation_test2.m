clear all;

% 1. Two satellites in circular orbits in equatorial plane moving in central gravity field - 1 day
%   The first one - analytical
%   The second one - numerical in central 
% 2. Two satellites in eccentric orbits 0.3 under central gravity field
%   For to situations solve analytic equation showing tau max
% 3. Derive a formula for tau max and check that if put two satellites at
%   the location it shows the right results
% 3. Apply NelderMead for the search of the mean anomaly corresponding to
% tau max and it coincides with analytic formula results

% make pre check on the Minimum Orbit Intersection Distance

consts = startup_formation_control();
global environment
environment = 'point mass';
syms t
assume(t, 'real');

dM = pi;
oe_o = [consts.rEarth + 700e3; 0; 0; 0; 0; 0];
oe_d = [consts.rEarth + 800e3; 0; 0; 0; 0; 0];
r_obs = 200e3;

d_phi = acos((oe_o(1)^2 + oe_d(1)^2 - r_obs^2) / (2 * oe_o(1) * oe_d(1)));
oe_o = [consts.rEarth + 700e3; 0; 0; 0; 0; 0];

n_o = sqrt(consts.muEarth / oe_o(1)^3);
n_d = sqrt(consts.muEarth / oe_d(1)^3);
tau_estimated = (2*pi - 2*d_phi) / (n_o - n_d);

rv_analytical_eq_circle_o = analytical_rv_function(oe_o, 'm, s', consts);
rv_analytical_eq_circle_d = analytical_rv_function(oe_d, 'm, s', consts);
dr = norm(rv_analytical_eq_circle_o(1:3) - rv_analytical_eq_circle_d(1:3));
f = dr - r_obs;
fplot(f, [1 10*consts.day2sec]);

% eq = dr - sym(r_obs) == sym(0);
% 
% tau_analytical = cutsom_solver(eq, t, 5e4);
% disp([tau_estimated, tau_analytical]);

variables_range = [consts.rEarth + 800e3, consts.rEarth + 900e3
                   0, 2*pi ];
               
variables0 = [consts.rEarth + 900e3; pi];

fun = @(variables)time_to_observation(variables, variables_range, rv_analytical_eq_circle_o, r_obs, consts);
M_d_optimal = fminsearch(fun,variables0);   

oe_d = [consts.rEarth + 800e3; 0; 0; 0; 0; M_d_optimal];
rv_analytical_eq_circle_d = analytical_rv_function(oe_d, 'm, s', consts);
dr = norm(rv_analytical_eq_circle_o(1:3) - rv_analytical_eq_circle_d(1:3));
eq = dr - sym(r_obs) == sym(0);
tau_neldermead = custom_solver(eq, t, 5e4);
disp([tau_estimated, tau_neldermead]);

% f = dr - r_obs;
% fplot(f, [1 10*consts.day2sec]);



