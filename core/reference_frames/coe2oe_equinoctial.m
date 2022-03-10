function oe_equinoctial = coe2oe_equinoctial(coe)

% coe [ a; e; i; RAAN, AOP, nu, M, u]
% oe_equinoctial [a,q1,q2,i,RAAN,lambda], where 
% q1 = e*cos(AOP), q2 = e*sin(AOP), lambda = AOP + M;

oe_equinoctial = [coe(1); coe(2)*cos(coe(5)); coe(2)*sin(coe(5));...
                  coe(3); coe(4); coe(5) + coe(7)];

end