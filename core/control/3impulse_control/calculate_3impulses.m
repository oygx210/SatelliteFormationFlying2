function maneuvers = calculate_3impulses(oe, delta_oe, consts)

% Reference: Mok et al, Impulsive Control of Satellite Formation Flying using Orbital Period Di?erence

da = delta_oe(1);
de = delta_oe(2);
di = delta_oe(3);
dRAAN = delta_oe(4);
dAOP = delta_oe(5);
dM = delta_oe(7);

theta_h = atan3(dRAAN*sin(oe(3)),di);

gamma = sqrt(oe(1)/consts.muEarth);

dV1 = [0;
       1/gamma * sqrt(di^2 + dRAAN^2*sin(oe(3))^2);
       0];

n = sqrt(consts.muEarth/oe(1)^3);
etta = sqrt(1-oe(2)^2);

dV2 = [n*oe(1)*etta/4 * (da/oe(1) + de/(1 + oe(2)));
       0;
       -n*oe(1)/4 * ((1+oe(2))^2/etta * (dAOP + dRAAN*cos(oe(3))) + dM)];

dV3 = [n*oe(1)*etta/4 * (da/oe(1) - de/(1 - oe(2))); 
       0;
       -n*oe(1)/4 * ((1-oe(2))^2/etta * (dAOP + dRAAN*cos(oe(3))) + dM)];

maneuvers.theta_h = theta_h;
maneuvers.dV1 = dV1;
maneuvers.dV2 = dV2;
maneuvers.dV3 = dV3;



