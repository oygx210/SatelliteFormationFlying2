function cost_matrix = calculate_cost_matrix(rv_ECI_initial, demonstration, formation, spacecraft, consts)
    
    delta_u = mod(n*seconds(demonstration.deployment_time - formation.orbit_epoch), 2*pi);
    target_orbit_coe = formation.coe; 
    target_orbit_coe(6) = mod(formation.coe(6) + delta_u, 2*pi);
    target_orbit = oe2rv(target_orbit_coe, consts);

    rv_orb_current = get_rv_from_analytic_HCW_solution(target_orbit, HCW_constants_current, consts);

    for i = 1:formation.N_sats
        rv_ECI_current(i*6-5:i*6,1) = orb2ECI(target_orbit, rv_orb_current(:,i), consts);
    end


end