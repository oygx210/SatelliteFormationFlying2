function inclination = get_SSO_inclination(a, ecc, consts)

% The formula for the RAAN drift is taken from D.A. Vallado Fundamentals of
% Astrodynamics, page 649, eq 9-37

% Inputs: constants, a - semi-major axis [m] and inclination
% Outputs: SSO orbit inclination [rad]

p = a.*(1-ecc.^2);
n = sqrt(consts.muEarth./a.^3);

inclination = acos( -(consts.EarthMeanMotion*2*p.^2) ./ (3*n*consts.rEarth^2*consts.J2)); 

end