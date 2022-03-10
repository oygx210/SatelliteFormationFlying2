function [t, rv] = formation_configuration_maintenance(rv_ECI, consts, spacecraft, formation)

option_tolerance = odeset('RelTol',1e-12,'AbsTol',1e-12);

leader.rv = rv_ECI(1:6);
leader.oe = rv2oe(leader.rv, consts);
leader.mean_oe = osc2mean(leader.oe, consts);

follower.rv = rv_ECI(7:12);
follower.oe = rv2oe(follower.rv, consts);
follower.mean_oe = osc2mean(follower.oe, consts);

follower.mean_oe_required = leader.mean_oe;
dM_required = 2*atan(formation.ISD/2/follower.mean_oe_required(1));
follower.mean_oe_required(7) = follower.mean_oe_required(7) - dM_required;

delta_mean_oe = follower.mean_oe_required - follower.mean_oe;

if (follower.mean_oe(5) > pi && follower.mean_oe_required(5) < pi)
    delta_mean_oe(5) = 2*pi - follower.mean_oe(5) + follower.mean_oe_required(5);
end

if (follower.mean_oe(7) > pi && follower.mean_oe_required(7) < pi)
    delta_mean_oe(7) = 2*pi - follower.mean_oe(7) + follower.mean_oe_required(7);
end

% Define mean_oe differences
da = delta_mean_oe(1);
de = delta_mean_oe(2);
di = delta_mean_oe(3);
dRAAN = delta_mean_oe(4);
dAOP = delta_mean_oe(5);
dM = delta_mean_oe(7);

% The maneuvers compensate difference in a,e,inc, RAAN, AOP, M

% Precalculated dVs
r = follower.rv(1:3);
v = follower.rv(4:6);
h_norm = vecnorm(cross(r,v));
r_norm = vecnorm(r);

dVh = [0;h_norm/r_norm*sqrt(di^2 + dRAAN^2*(sin(follower.mean_oe(3)))^2); 0];

dVp = [0; 0; 0];
n = sqrt(consts.muEarth/follower.mean_oe(1)^3);
etta = sqrt(1 - follower.mean_oe(2)^2);
dVp(1) = n*follower.mean_oe(1)*etta/4 * (da/follower.mean_oe(1) + de/(1 + follower.mean_oe(2)));
dVp(3) = -n*follower.mean_oe(1)/4 * ((1 + follower.mean_oe(2))^2/etta * (dAOP + dRAAN*cos(follower.mean_oe(3))) + dM);

dVa = [0; 0; 0];
etta = sqrt(1 - follower.mean_oe(2)^2);
dVa(1) = n*follower.mean_oe(1)*etta/4 * (da/follower.mean_oe(1) - de/(1 - follower.mean_oe(2)));
dVa(3) = n*follower.mean_oe(1)/4 * ((1 - follower.mean_oe(2))^2/etta * (dAOP + dRAAN*cos(follower.mean_oe(3))) + dM);

theta_h = atan2(dRAAN*sin(follower.mean_oe(3)),di);

if theta_h < 0
    theta_h = 2*pi + theta_h;
end

if follower.oe(8) < theta_h
    t_span = (theta_h - follower.oe(8))*sqrt(follower.mean_oe(1)^3/consts.muEarth);
else 
    t_span = (2*pi - follower.oe(8)  + theta_h)*sqrt(follower.mean_oe(1)^3/consts.muEarth);
end

rv_init = rv_ECI;
[t_vec1, rv_ECI1] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, option_tolerance);
t = t_vec1;
rv_ECI1 = rv_ECI1';
rv_ECI = rv_ECI1;

follower.rv = rv_ECI(7:12, end);
follower.oe = rv2oe(follower.rv, consts);
follower.mean_oe = osc2mean(follower.oe, consts);


rv_ECI(7:12,end) = orb2ECI(rv_ECI(7:12,end), [0;0;0;dVh] ,consts);

% plot_coe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'd');
% plot_mean_oe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'd');

follower.oe = rv2oe(rv_ECI(7:12,end), consts);
follower.mean_oe = osc2mean(follower.oe, consts);

% Then the second impulse is planned to be at periapsis
    
t_span = (2*pi - follower.mean_oe(6))*sqrt(follower.mean_oe(1)^3/consts.muEarth);

rv_init = rv_ECI(:,end);
[t_vec2, rv_ECI2] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, option_tolerance);
t = [t; t_vec2 + t(end)];
rv_ECI2 = rv_ECI2';
rv_ECI = [rv_ECI, rv_ECI2];

follower.rv = rv_ECI(7:12, end);
follower.oe = rv2oe(follower.rv, consts);
follower.mean_oe = osc2mean(follower.oe, consts);

rv_ECI(7:12, end) = orb2ECI(rv_ECI(7:12, end), [0;0;0;dVp], consts);
follower.oe = rv2oe(rv_ECI(7:12,end), consts);
follower.mean_oe = osc2mean(follower.oe, consts);

t_span = pi*sqrt(follower.mean_oe(1)^3/consts.muEarth);
% if follower.mean_oe(8) < pi
%     t_span = (pi - follower.mean_oe(8))*sqrt(follower.mean_oe(1)^3/consts.muEarth);
% else
%     t_span = (2*pi - follower.mean_oe(8) + pi)*sqrt(follower.mean_oe(1)^3/consts.muEarth);
% end

rv_init = rv_ECI(:,end);
[t_vec3, rv_ECI3] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, option_tolerance);
t = [t; t_vec3 + t(end)];
rv_ECI3 = rv_ECI3';
rv_ECI = [rv_ECI, rv_ECI3];

follower.rv = rv_ECI(7:12, end);
follower.oe = rv2oe(follower.rv, consts);
follower.mean_oe = osc2mean(follower.oe, consts);

rv_ECI(7:12, end) = orb2ECI(rv_ECI(7:12, end), [0; 0; 0; dVa], consts);

leader.mean_oe = osc2mean(rv2oe(rv_ECI(1:6, end), consts),consts);
follower.mean_oe = osc2mean(rv2oe(rv_ECI(7:12, end), consts),consts);
delta_mean_oe2 = leader.mean_oe - follower.mean_oe;
deltas = [delta_mean_oe2, delta_mean_oe];

plot(t, vecnorm(rv_ECI(1:3,:) - rv_ECI(7:9,:)));
plot_coe_diff(t,  rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'd');
plot_mean_oe_diff(t,  rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'd');



rv = rv_ECI;