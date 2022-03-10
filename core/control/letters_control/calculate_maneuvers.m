function maneuvers = calculate_maneuvers(rv_ECI_current, rv_ECI_required, consts)


moe = rv2moe(rv_ECI_current, consts);
coe = vecRV2OE(rv_ECI_current, consts);

moe_required = rv2moe(rv_ECI_required, consts);

delta_moe = moe_required - moe;

% here might be some tricky things
if (moe(5) > pi && moe_required(5) < pi && delta_moe(5) > pi)
    delta_moe(5) = 2*pi - moe(5) + moe_required(5);
end
if (moe(7) > pi && moe_required(7) < pi && delta_moe(7) > pi)
    delta_moe(7) = 2*pi - moe(7) + moe_required(7);
end

da = delta_moe(1);
de = delta_moe(2);
di = delta_moe(3);
dRAAN = delta_moe(4);
dAOP = delta_moe(5);
dM = delta_moe(7);

theta_h = atan3(dRAAN*sin(moe(3)),di);

if theta_h < coe(8) 
   t_span = (2*pi - coe(8) + theta_h)*sqrt(moe(1)^3/consts.muEarth);
else
   t_span = (theta_h - coe(8))*sqrt(moe(1)^3/consts.muEarth);
end

gamma = sqrt(moe(1)/consts.muEarth);

dV1 = [0;
       1/gamma * sqrt(di^2 + dRAAN^2*sin(moe(3))^2);
       0];

n = sqrt(consts.muEarth/moe(1)^3);
etta = sqrt(1-moe(2)^2);

dV2 = [n*moe(1)*etta/4 * (da/moe(1) + de/(1 + moe(2)));
       0;
       -n*moe(1)/4 * ((1+moe(2))^2/etta * (dAOP + dRAAN*cos(moe(3))) + dM)];

dV3 = [n*moe(1)*etta/4 * (da/moe(1) - de/(1 - moe(2))); 
       0;
       -n*moe(1)/4 * ((1-moe(2))^2/etta * (dAOP + dRAAN*cos(moe(3))) + dM)];
   
maneuvers = [t_span; dV1; dV2; dV3];

end