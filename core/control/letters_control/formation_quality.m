function [value_out,isterminal,direction] = formation_quality(t, rv, consts, spacecraft, formation)

% formation_quality is an event function that tracks error in position for
% all satellite in formation and tracks ISL distance between target and
% chaser(i)

% Reference orbits
rv_orb_predicted = get_rv_from_analytic_HCW_solution(rv(1:6), formation.geometry, consts);
rv = reshape(rv, [6,formation.N_sats]);

for i = 1:formation.N_sats
    rv_orb(:,i) = vecECI2orb(rv(:,1), rv(:,i), consts);
end

isterminal = 1;
direction = 0;

number_of_broken_pixels = sum(vecnorm(rv_orb(1:3,:) - rv_orb_predicted(1:3,:))./vecnorm(rv_orb_predicted(1:3,:)) > spacecraft.tracking_error);

if number_of_broken_pixels > 0
    value_out = 0;
else
    value_out = 1;
end


end