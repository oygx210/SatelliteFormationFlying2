function rv_prime = J2_atmo(t, rv, consts, spacecraft)

% Right-hand side for the equation of motion of a satellite taking into
% account:
% 1. Oblated Earth
% 2. Atmospheric drag (piecewise exponential atmospheric density function)

R = consts.rEarth_equatorial;
J2 = consts.J2;
mu = consts.muEarth;
delta = 3/2*J2*mu*R^2;

x = rv(1);
y = rv(2); 
z = rv(3);
r = [rv(1); rv(2); rv(3)];

% J2 effect
acceleration_J2 = delta*r/norm(r)^5*(5*z^2/norm(r)^2 - 1)- 2*delta/norm(r)^5*[0; 0; z];

% Atmospheric drag
vRelativeECI = rv(4:6) - cross(consts.wEarth, rv(1:3))';
rhoAtmo = CIRA72(consts, (norm(rv(1:3)) - consts.rEarth));
acceleration_AD = - 0.5 * spacecraft.Cdrag * spacecraft.DragArea / spacecraft.mass * rhoAtmo * vRelativeECI * norm(vRelativeECI); 


rv_prime = [rv(4); rv(5); rv(6);...
            -consts.muEarth*rv(1)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(1) + acceleration_AD(1);...
            -consts.muEarth*rv(2)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(2) + acceleration_AD(2);...
            -consts.muEarth*rv(3)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(3) + acceleration_AD(3)];
        
end