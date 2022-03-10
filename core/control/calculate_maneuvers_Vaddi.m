function [t_span, dV1, dV2] = calculate_maneuvers_Vaddi(rv_ECI_current, rv_ECI_required, consts)

% The maneuvers are adopted from Vaddi et al., Formation Establishment and Reconfiguration Using Impulsive Control

coe_current = rv2oe(rv_ECI_current, consts);
coe_required = rv2oe(rv_ECI_required, consts);

eqoe_current = coe2oe_equinoctial(coe_current);
eqoe_required = coe2oe_equinoctial(coe_required);

% coe [ a; e; i; RAAN, AOP, nu, M, u]
% oe_equinoctial [a,q1,q2,i,RAAN,lambda], where 
% q1 = e*cos(AOP), q2 = e*sin(AOP), lambda = AOP + M;

delta_eqoe = eqoe_required - eqoe_current;
A = [eqoe_current, eqoe_required, delta_eqoe];

if (eqoe_current(5) > pi && eqoe_required(5) < pi && abs(delta_eqoe(5)) > pi)
    delta_eqoe(5) = 2*pi - eqoe_current(5) + eqoe_required(5);
elseif (eqoe_current(5) < pi && eqoe_required(5) > pi && abs(delta_eqoe(5)) > pi)
    delta_eqoe(5) = 2*pi - eqoe_required(5) + eqoe_current(5);
end

dq1 = delta_eqoe(2);
dq2 = delta_eqoe(3);
di = delta_eqoe(4);
dRAAN = delta_eqoe(5);

gamma = sqrt(coe_current(1)/consts.muEarth);

theta_h = atan3(dRAAN*sin(coe_current(3)),di);

if theta_h < coe_current(8) 
   t_span = (2*pi - coe_current(8) + theta_h)*sqrt(coe_current(1)^3/consts.muEarth);
else
   t_span = (theta_h - coe_current(8))*sqrt(coe_current(1)^3/consts.muEarth);
end

dV1 = [0;
       1/gamma * sqrt(di^2 + dRAAN^2*sin(coe_current(3))^2);
       sqrt(dq1^2 + dq2^2)/2/gamma];

dV2 = [0;
       0;
       -sqrt(dq1^2 + dq2^2)/2/gamma];
% disp([dV1,dV2]);
% [t_span; dV1; dV2];

end