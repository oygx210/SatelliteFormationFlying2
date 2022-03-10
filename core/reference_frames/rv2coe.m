function coe = rv2coe(rv, consts)

mu = consts.muEarth;
tolerance = 1e-7;

% Author: Shamil Biktimirov
% Reference: Vallado D., Fundamentals of Astrodynamics and Applications,
% fourth edition, pp. 112-116

% Input: rv - satellite state vector [m, m/s]
%        consts - constants for formation constrol model
% Output: a [m] - semi-major axis
%         ecc [-] - eccentricity
%         inc [rad] - inclination, 0:2*pi
%         RAAN [rad] - right ascencion of ascending node, 0:2*pi
%         AOP [rad] - argument of pericenter, 0:2*pi
%         nu [rad] - true anomaly, 0:2*pi
%         M [rad] - mean anomaly, 0:2*pi
%         u [rad] - argument of latitude, 0:2*pi

% Calculating h, angular momentumv vector 
r = rv(1:3);
v = rv(4:6);
r_norm = vecnorm(r);
v_norm = vecnorm(v);

h = cross(r,v);
h_norm = vecnorm(h);
h_unit = h./h_norm;

% Calculating n = cross(K, h) - vector to node, where K = [0; 0; 1]
n = [-h(2); h(1); 0];
n_norm = vecnorm(n);
n_unit = n./n_norm;

% Calculating eccentricity vector
e = ((v_norm^2 - mu/r_norm)*r - dot(r,v)*v)/mu;
ecc = vecnorm(e);

if ecc >= 1
    warning('The orbit is not elliptic');
end

% Calculating semi-major axis
a = (2/r_norm - v_norm^2/mu)^(-1);

% Calculating inc, RAAN, AOP, nu
inc = acos(h_unit(3));

RAAN = acos(n_unit(1));
if n(2) < 0
    RAAN = 2*pi - RAAN;
end

AOP = acos(dot(n, e)/ n_norm /ecc);
if e(3) < 0
    AOP = 2*pi - AOP;
end

nu = acos(dot(e, r)/ecc/ r_norm);
if dot(r,v) < 0 
    nu  = 2*pi - nu;
end

if ecc < tolerance
    AOP = 0;
    
    nu = acos(dot(n, r)/ n_norm/ r_norm);
    if r(3) < 0
        nu = 2*pi - nu;
    end

end

u = acos(dot(n, r)/ n_norm/ r_norm);
if r(3) < 0
    u = 2*pi - u;
end

[~, M] = newtonnu (ecc, nu);

coe = [a; ecc; inc; RAAN; AOP; nu; M; u];




