clear all;

%% Initialization
global environment 
syms t
spacecraft = [];
consts = startup_formation_control();

rv_target_orbit = [-3670e3; -3870e3; 4400e3; 4.7e3; -7.4e3; 1e3]; % Book example
target_orbit = rv2oe(rv_target_orbit, consts);
oe = [target_orbit(1:5);target_orbit(7)];

r_target_an2 = analytical_rv_function(oe, consts);
r_test = double(subs(r_target_an2, t, seconds(days(4))));

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
environment = 'J2';

[t_vec, rv_ECI_numerical_J2] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [1:4*consts.day2sec], rv_target_orbit, options_precision);
rv_ECI_numerical_J2 = rv_ECI_numerical_J2'; 

r_book = [9672e3; 4320e3; -8691e3]; % Result of analytical propagation from the book
dr_book_analytical = vecnorm(r_book - r_test(1:3));
dr_book_num = vecnorm(r_book - rv_ECI_numerical_J2(1:3,end));