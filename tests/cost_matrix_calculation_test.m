function cost_matrix = cost_matrix_calculation_test(rv_ECI_initial, demonstration, formation, spacecraft, consts)
    
    options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
    global T;

    formation.geometry = [];
    
    f = waitbar(0,'Calculating reconf cost matrix & reconf parameters');    
    state_after_correction = zeros(formation.N_active_sats,formation.N_active_sats,6);
    dt_impulsive_correction = zeros(formation.N_active_sats,formation.N_active_sats);

    for i = 1:(size(demonstration.HCW_constants,2)-1)
        f = waitbar(i/(size(demonstration.HCW_constants,2)-1));        
        for j = 1:(size(demonstration.HCW_constants,2)-1)
            rv_ECI = [];
            t_vec = [];
            formation.geometry = [];

            formation.geometry(:,1) = [0; 0; 0; 0];
            formation.geometry(:,2) = demonstration.HCW_constants(:,1+j);
            formation.N_sats = size(formation.geometry,2);
            formation.N_active_sats = size(formation.geometry,2) - 1;
            rv_ECI(1:6,1) = rv_ECI_initial(1:6,1);            
            rv_ECI(7:12,1) = rv_ECI_initial(6*(i+1)-5:6*(i+1),1);

            mode = 3;
            [t_vec, rv_ECI, impulsive_maneuvers_dV1(i,j), ~] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);
            state_after_correction1(i,j,:) = rv_ECI(:,end);
            dt_impulsive_correction1(i,j) = t_vec(end);
            geometry(i,j,:) = formation.geometry(:,2);
        end
    end
    
    for i = 1:(size(demonstration.HCW_constants,2)-1)
        f = waitbar(i/(size(demonstration.HCW_constants,2)-1));        
        for j = 1:(size(demonstration.HCW_constants,2)-1)
            rv_ECI = [];
            t_vec = [];
            formation.geometry = [];

            formation.geometry(:,1) = [0; 0; 0; 0];
            formation.geometry(:,2) = demonstration.HCW_constants(:,1+j);
            formation.N_sats = size(formation.geometry,2);
            formation.N_active_sats = size(formation.geometry,2) - 1;
            rv_ECI(1:6,1) = rv_ECI_initial(1:6,1);            
            rv_ECI(7:12,1) = rv_ECI_initial(6*(i+1)-5:6*(i+1),1);

            mode = 3;
            [t_vec, rv_ECI, impulsive_maneuvers_dV1(i,j), ~] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);
            state_after_correction1(i,j,:) = rv_ECI(:,end);
            dt_impulsive_correction1(i,j) = t_vec(end);
            geometry(i,j,:) = formation.geometry(:,2);
        end
    end

%     for i = 1:(size(demonstration.HCW_constants,2)-1)
%         for j = 1:(size(demonstration.HCW_constants,2)-1)
%             error_dV(i,j) = impulsive_maneuvers_dV2(i,j) - impulsive_maneuvers_dV1(i,j);
%             
%         end
%     end
    
    close(f);
end
