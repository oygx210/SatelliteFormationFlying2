function [t, rv] = formation_maintenance_3impulses_Mok_coe_based(rv_ECI, consts, spacecraft, formation)

% Mok et al.,  Impulsive Control of Satellite FormationFlying using Orbital Period Difference

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% 1. We determine difference in classical orbital elements
leader.coe = vecRV2OE(rv_ECI(1:6), consts);
follower.coe = vecRV2OE(rv_ECI(7:12), consts);

dM = 2*atan(formation.ISD/2/leader.coe(1));

follower.coe_required = leader.coe;
follower.coe_required(7) = leader.coe(7) - dM;

follower.delta_coe = follower.coe_required - follower.coe;
if (follower.coe(5) > pi && follower.coe_required(5) < pi)
    follower.delta_coe(5) = 2*pi - follower.coe(5) + follower.coe_required(5);
end
if (follower.coe(7) > pi && follower.coe_required(7) < pi)
    follower.delta_coe(7) = 2*pi - follower.coe(7) + follower.coe_required(7);
end

% We calculate three impulses and location for di&dRAAN correction
maneuvers = calculate_3impulses(rv_ECI(7:12), follower.coe, follower.delta_coe, consts);

% To correct the difference in classical orbital elements we use the scheme:
% 1. based on current argument of latitude and the required to perform di&dRAAN correction we calculate time to apply the impulse
% 2. We propagate both satellites untill the follower reaches the required argument of latitude to correct di&dRAAN
% 3. We apply the impulse and recalculate classical oe -> based on new coe we  calculate time to reach perigee of follower's orbit
% 4. We propagate to perigee of followers orbit and appy the second impulse there
% 5. We update coe of the follower and fly to apogee where the third impulse should be applied

if maneuvers.theta_h < follower.coe(8) 
   t_span = (2*pi - follower.coe(8) + maneuvers.theta_h)*sqrt(follower.coe(1)^3/consts.muEarth);
else
   t_span = (maneuvers.theta_h - follower.coe(8))*sqrt(follower.coe(1)^3/consts.muEarth);
end

% waiting to apply the first implulse at theta
rv_init = rv_ECI;
[t_vec1, rv_ECI1] = ode45(@(t, rv) central_gravity_Formation(t, rv, consts), 0:t_span, rv_init, options);
t = t_vec1;
rv_ECI1 = rv_ECI1';
rv_ECI = rv_ECI1;

follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV1], consts);
follower.coe = vecRV2OE(follower.rv_ECI, consts);
rv_ECI(7:12,end) = follower.rv_ECI;

t_span = (2*pi - follower.coe(6)) * sqrt(follower.coe(1)^3/consts.muEarth);

rv_init = rv_ECI(:,end);
[t_vec2, rv_ECI2] = ode45(@(t, rv) central_gravity_Formation(t, rv, consts), 0:t_span, rv_init, options);
rv_ECI2 = rv_ECI2';
rv_ECI = [rv_ECI, rv_ECI2];
t = [t; t_vec2 + t(end)];

follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV2], consts);
follower.coe = vecRV2OE(follower.rv_ECI,consts);
rv_ECI(7:12,end) = follower.rv_ECI;

% time before apocenter considering that pericenter has slightly shifted after maneuvering from its previous location
if follower.coe(6) > pi
    t_span = (3*pi - follower.coe(6)) * sqrt(follower.coe(1)^3/consts.muEarth);
else
    t_span = (pi - follower.coe(6)) * sqrt(follower.coe(1)^3/consts.muEarth);
end

rv_init = rv_ECI(:,end);
[t_vec3, rv_ECI3] = ode45(@(t, rv) central_gravity_Formation(t, rv, consts), 0:t_span, rv_init, options);
rv_ECI3 = rv_ECI3';
rv_ECI = [rv_ECI, rv_ECI3];
t = [t; t_vec3 + t(end)];

follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV3], consts);
rv_ECI(7:12,end) = follower.rv_ECI;

plot_coe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'h');
plot_ISD(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'h');

rv = rv_ECI;

end