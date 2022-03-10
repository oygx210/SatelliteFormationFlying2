function [RAAN_dot, AOP_dot] = get_J2_precession(oe, consts)


C = - 3 / 2 * (sqrt(consts.muEarth)*consts.J2*consts.rEarth_equatorial^2);

RAAN_dot = C / ((1 - oe(2)^2)^2*oe(1)^(7/2)) * cos(oe(3));

AOP_dot = C / ((1 - oe(2)^2)^2 * oe(1)^(7/2)) * (5/2 * sin(oe(3))^2 - 2);


end