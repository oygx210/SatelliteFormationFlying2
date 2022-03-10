clear all;

% ISD and SR are defined for a give orbit
consts = startup_formation_control;

k_rev2rep = 14;
k_day2rep = 1;
% orbital elements
orbit_epoch_GD = datetime(2023, 1, 1, 0, 0, 0);
orbit_epoch_JD = juliandate(orbit_epoch_GD);

oe = zeros(6,1);
oe(1) = get_SSO_RGT_orbit_sma(k_rev2rep, k_day2rep, consts);
oe(2) = 0; % ecc, 
oe(3) = get_SSO_inclination(oe(1), oe(2), consts);
orbit_epoch_GD = datetime(orbit_epoch_JD, 'convertfrom','juliandate');
oe(4) = get_RAAN_for_terminator_orbit(orbit_epoch_GD);
oe(5) = 0; % AOP, deg
oe(6) = 0; % M, deg - actually argument of latitude
% the actual orbital elements should correspond to the optimal coverage

rv_ECI = oe2rv(oe, consts);

theta_min = deg2rad(10);
b = consts.rEarth / oe(1);
gamma_max = asin(b * cos(theta_min));
beta_max = pi/2 - gamma_max - theta_min;

d_max = oe(1) / cos(theta_min) * sin(pi/2 - gamma_max - theta_min);
ISD_min = deg2rad(1/60) * d_max;

rho = ISD_min * 10;

formation_geometry = [rho; sqrt(3) / 2 * rho; 0; 0];
rv_HCW = get_rv_from_analytic_HCW_solution(rv_ECI, formation_geometry, consts);

rv_ECI(7:12) = orb2ECI(rv_ECI, rv_HCW, consts);


