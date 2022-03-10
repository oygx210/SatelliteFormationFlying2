function rv_prime = J2(t, rv, consts)

% Right-hand side for the equation of motion of a satellite orbiting
% oblated Earth

R = consts.rEarth_equatorial;
J2 = consts.J2;
mu = consts.muEarth;

delta = 3/2*J2*mu*R^2;

x = rv(1);
y = rv(2); 
z = rv(3);
r=[rv(1);rv(2);rv(3)];

acceleration_J2 = delta*r/norm(r)^5*(5*z^2/norm(r)^2 - 1)- 2*delta/norm(r)^5*[0; 0; z];


rv_prime = [rv(4); rv(5); rv(6);...
            -consts.muEarth*rv(1)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(1);...
            -consts.muEarth*rv(2)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(2);...
            -consts.muEarth*rv(3)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2) + acceleration_J2(3)];
end