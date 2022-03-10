function [t_span, dV_perigee, dV_apogee] = calculate_sma_ecc_correction(rv_ECI_current, rv_ECI_required, consts)

% the function calculated dV to adjust sma and ecc
moe = rv2moe(rv_ECI_current, consts);
coe = rv2oe(rv_ECI_current, consts);

moe_required = rv2moe(rv_ECI_required, consts);

delta_moe = moe_required - moe;

n = sqrt(consts.muEarth/moe(1)^3);
etta = sqrt(1-moe(2)^2);

a = moe(1);
ecc = moe(2);
da = delta_moe(1);
decc = delta_moe(2);

if abs(da) < 1.5
    da = 0;
    decc = 0;
end

dV_perigee = [n*a*etta/4*(da/a + decc/(1+ecc));
              0;
              0];

dV_apogee = [n*a*etta/4*(da/a - decc/(1-ecc));
              0;
              0];

t_span = (2*pi - coe(6))*sqrt(moe(1)^3/consts.muEarth);

end