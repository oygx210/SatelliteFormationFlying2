function rv_prime = central_gravity_Formation(t, rv, consts)

N_sats = length(rv)/6;
rv = reshape(rv, [6,N_sats]);

r_prime = [rv(4:6,:); zeros(3,N_sats)];
r_norm = vecnorm(rv(1:3,:));

% Central gravity
A_cg = [zeros(3,N_sats); -consts.muEarth./(r_norm.^3).*rv(1:3,:)];

rv_prime = r_prime + A_cg;
rv_prime = reshape(rv_prime, [6*N_sats,1]);

end
