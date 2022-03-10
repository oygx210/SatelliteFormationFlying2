function [t, rv_ECI, maneuvers_dV, maneuvers_t] = formation_maintenance_3impulses_Mok_moe_based(rv_ECI, consts, spacecraft, formation)

% The maneuvers are adopted from Mok et al., Impulsive Control of Satellite Formation Flying using Orbital Period Difference
% The idea to utilize mean orbital elements is adopted from Schaub et al., Impulsive Feedback Control to Establish Specific Mean Orbit Elements of Spacecraft Formations

% To correct the difference in mean orbital elements we utilize the following scheme:
% 1. based on current argument of latitude and the required to perform di&dRAAN correction we calculate time to apply the impulse
% 2. We propagate both satellites untill the follower reaches the required argument of latitude to correct di&dRAAN
% 3. We apply the impulse and recalculate classical oe -> based on new coe we  calculate time to reach perigee of follower's orbit
% 4. We propagate to perigee of followers orbit and appy the second impulse there
% 5. We update coe of the follower and fly to apogee where the third impulse should be applied

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% 1. Determines difference in mean orbital elements(moe)
leader.coe = vecRV2OE(rv_ECI(1:6), consts); leader.moe = rv2moe(rv_ECI(1:6), consts);
follower.coe = vecRV2OE(rv_ECI(7:12), consts); follower.moe = rv2moe(rv_ECI(7:12), consts);

dM_moe = 2*atan(formation.ISD/2/leader.moe(1));

follower.moe_required = leader.moe;
follower.moe_required(7) = leader.moe(7) - dM_moe;

follower.delta_moe = follower.moe_required - follower.moe;

if (follower.moe(5) > pi && follower.moe_required(5) < pi)
    follower.delta_moe(5) = 2*pi - follower.moe(5) + follower.moe_required(5);
end
if (follower.moe(7) > pi && follower.moe_required(7) < pi)
    follower.delta_moe(7) = 2*pi - follower.moe(7) + follower.moe_required(7);
end

%2. Calculates dV for three impuslive maneuvers and argument of latitude for first correction
% dV1 corrects di,dRAAN at argument of latitude theta_h
% dV2 & dV3 corrects da, de, dAOP, dM at perigee and apogee respectively

maneuvers = calculate_3impulses(follower.moe, follower.delta_moe, consts);
maneuvers_dV = [maneuvers.dV1; maneuvers.dV2; maneuvers.dV3];

if maneuvers.theta_h < follower.coe(8) 
   t_span = (2*pi - follower.coe(8) + maneuvers.theta_h)*sqrt(follower.moe(1)^3/consts.muEarth);
else
   t_span = (maneuvers.theta_h - follower.coe(8))*sqrt(follower.moe(1)^3/consts.muEarth);
end

% waiting to apply the first implulse at theta
rv_init = rv_ECI;
[t_vec1, rv_ECI1] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
t = t_vec1;
rv_ECI1 = rv_ECI1';
rv_ECI = rv_ECI1;

maneuvers_t(1,1) = t(end);

follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV1], consts);
follower.coe = vecRV2OE(follower.rv_ECI, consts); follower.moe = rv2moe(follower.rv_ECI, consts);
rv_ECI(7:12,end) = follower.rv_ECI;

t_span = (2*pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);

rv_init = rv_ECI(:,end);
[t_vec2, rv_ECI2] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
rv_ECI2 = rv_ECI2';
rv_ECI = [rv_ECI, rv_ECI2];
t = [t; t_vec2 + t(end)];

maneuvers_t(2,1) = t(end);
follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV2], consts);
follower.coe = vecRV2OE(follower.rv_ECI,consts); follower.moe = rv2moe(follower.rv_ECI,consts);
rv_ECI(7:12,end) = follower.rv_ECI;

% time before apocenter considering that pericenter has slightly shifted after maneuvering from its previous location
if follower.coe(6) > pi
    t_span = (3*pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);
else
    t_span = (pi - follower.coe(6)) * sqrt(follower.moe(1)^3/consts.muEarth);
end

rv_init = rv_ECI(:,end);
[t_vec3, rv_ECI3] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:t_span, rv_init, options);
rv_ECI3 = rv_ECI3';
rv_ECI = [rv_ECI, rv_ECI3];
t = [t; t_vec3 + t(end)];

maneuvers_t(3,1) = t(end);
follower.rv_ECI = orb2ECI(rv_ECI(7:12,end), [0; 0; 0; maneuvers.dV3], consts);
rv_ECI(7:12,end) = follower.rv_ECI;

plot_mean_oe_diff(t, rv_ECI(1:6,:), rv_ECI(7:12,:), consts, 'h');
end