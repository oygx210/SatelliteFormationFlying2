function rv_prime = J2_atmo_Formation(t, rv, consts, spacecraft)

%% The function builds right-hand side for the equation of motion of N satellites. 

% Right-hand side for the equation of motion takes into account:
% 1. Oblated Earth
% 2. Atmospheric drag (piecewise exponential atmospheric density function)

N_sats = length(rv)/6;
rv = reshape(rv, [6,N_sats]);

R = consts.rEarth_equatorial;
J2 = consts.J2;
mu = consts.muEarth;
delta = 3/2*J2*mu*R^2;

r_prime = [rv(4:6,:); zeros(3,N_sats)];
r_norm = vecnorm(rv(1:3,:));

% Central gravity
A_cg = [zeros(3,N_sats); -consts.muEarth./(r_norm.^3).*rv(1:3,:)];

% J2 perturbation
% check the formula Vallado or Chazov
A_j2 = delta./(r_norm.^5).*rv(1:3,:).*(5*rv(3,:).^2./r_norm.^2 - 1) - 2*delta./r_norm.^5.*[zeros(1,N_sats);zeros(1,N_sats);rv(3,:)];
A_j2 = [zeros(3,N_sats); A_j2];

% Atmospheric drag
vRelativeECI = rv(4:6,:) - cross(ones(1,N_sats).*consts.wEarth, rv(1:3,:));
rhoAtmo = CIRA72(consts, (vecnorm(rv(1:3,:)) - consts.rEarth));
A_ad = - 0.5 * spacecraft.Cdrag * spacecraft.DragArea / spacecraft.mass * rhoAtmo .* vRelativeECI .* vecnorm(vRelativeECI); 
A_ad = [zeros(3,N_sats); A_ad];

rv_prime = r_prime + A_cg + A_j2 + A_ad;
rv_prime = reshape(rv_prime, [6*N_sats,1]);

end
