function [t_vec, rv_ECI, HCW_constants_assigned, matchMatrix, maneuvers_reconfiguration, continuous_maneuvers_dV, fuel_consumption_reconfiguration] = reconfiguration(rv_ECI, required_geometry, Cost_matrix_dV, formation, spacecraft, consts)

    global T    
    maneuvers_reconfiguration = [];

    for i = 1:formation.N_active_sats
        spacecraft_wet_mass_updated = (spacecraft.dry_mass + formation.fuel_level(i,end))*exp(-Cost_matrix_dV(i,:)/spacecraft.thruster_Isp/consts.g);
        reconfiguration_matrix_fuel(i,:) = (spacecraft.dry_mass + formation.fuel_level(i,end))*ones(1, formation.N_active_sats) - spacecraft_wet_mass_updated;
    end

   % Solving assignment problem and assigning satellites to new set of reference trajectories
   [matchMatrix, ~] = maneuverAssignment(reconfiguration_matrix_fuel, formation.fuel_level(:,end));    

    for i = 1:(formation.N_sats-1)
        formation.geometry(:,matchMatrix(i,1)+1) = required_geometry(:,matchMatrix(i,2)+1);
    end

    HCW_constants_assigned = formation.geometry;

    mode = 1;
    [t_vec_m, rv_ECI, impulsive_maneuvers_dV, ~] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);
    t_reconf_start = t_vec_m(1) + T;
    t_vec = t_vec_m + T;
    T = t_vec(end);

    mode = 1;
    [t_vec_out, rv_ECI_out, continuous_maneuvers_dV, ~] = continuous_post_correction(rv_ECI(:,end), consts, spacecraft, formation, mode);
    t_vec = [t_vec; t_vec_out(2:end)];
    rv_ECI = [rv_ECI, rv_ECI_out(:,2:end)];    
    T = t_vec(end);

    maneuvers_reconfiguration = impulsive_maneuvers_dV + continuous_maneuvers_dV;
    fuel_level_reconfiguration = write_off_fuel(formation.fuel_level,[1:formation.N_active_sats],maneuvers_reconfiguration, spacecraft, consts);
    fuel_consumption_reconfiguration = formation.fuel_level - fuel_level_reconfiguration(:,end);
    

    disp(['Reconfiguration took ' num2str((t_vec(end) - t_reconf_start)/60) ' minutes, impulsive corrections: ' num2str(t_vec_m(end)/60) ' min, continuous control: ' num2str((t_vec_out(end) - t_vec_out(1))/60) ' min']);
