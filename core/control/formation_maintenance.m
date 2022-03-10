function [t, rv] = formation_maintenance(rv_ECI, consts, spacecraft, formation)

% here the Bi-Impulsive maneuvers are to be done in order to adjust
% reference orbit for satellites in formation

rv_orb_predicted = get_rv_from_analytic_HCW_solution(rv_ECI(1:6), formation.geometry, consts);

for i = 1:formation.N_sats
    rv_ECI_predicted(:,i) = orb2ECI(rv_ECI(1:6), rv_orb_predicted(:,i), consts);
end

rv_ECI = reshape(rv_ECI, [6,formation.N_sats]);

schedule = schedule_burns(rv_ECI(:,1), consts, formation);

dV = calculate_maneuver(rv_ECI(1:6,1), rv_ECI, rv_ECI_predicted, consts);


for i = 1:formation.N_sats
    [t, rv(i*6-5:i*6,:)] = singlesat_orbit_maintenance(rv_ECI(:,1), rv_ECI(:,i), schedule(:,i), dV(1:3,i), dV(4:6,i), consts, spacecraft);
end

end