

Torb = 2*pi*sqrt(formation.coe(1)^3/consts.muEarth);

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, rv_ECI] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[0 Torb], formation.rv, options_precision);        
rv_ECI = rv_ECI';

r = rv_ECI(1:3,:);
v = rv_ECI(4:6,:);

t_desired = 1:1:Torb;

s = spline(t,rv_ECI,t_desired);

figure;
plot3(r(1,:), r(2,:), r(3,:), 'ok');
hold on;
plot3(s(1,:), s(2,:), s(3,:), '*r');

