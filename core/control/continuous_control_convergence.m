function [value_out,isterminal,direction] = continuous_control_convergence(t, rv, consts, spacecraft, formation)

% Reference orbits
rv_orb_predicted = get_rv_from_analytic_HCW_solution(rv(1:6), formation.geometry, consts);
rv = reshape(rv, [6,formation.N_sats]);

for i = 1:formation.N_sats
    rv_orb(:,i) = ECI2orb(rv(:,1), rv(:,i), consts);
end

isterminal = 1;
direction = 0;

convergence = 100; % meters
number_of_converged_pixels = sum(vecnorm(rv_orb(1:3,:) - rv_orb_predicted(1:3,:)) < convergence);

if number_of_converged_pixels == formation.N_sats
    value_out = 0;
else
    value_out = 1;
end

end