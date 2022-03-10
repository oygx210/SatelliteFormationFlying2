function Reconfiguration_matrix = get_reconfiguration_matrix(rv_ECI, HCW_constants_required, formation, spacecraft, consts)
    
    Reconfiguration_matrix = zeros(formation.N_sats-1);
    f = waitbar(0,'Calculating reconfiguration matrix');
    
    formation.geometry = [];
    
    for i = 1:(size(HCW_constants_required,2)-1)
        f = waitbar(i/(size(HCW_constants_required,2)-1));
        
        for j = 1:(size(HCW_constants_required,2)-1)
            formation.geometry(:,1) = [0; 0; 0; 0];
            formation.geometry(:,2) = HCW_constants_required(:,1+j);
            formation.N_sats = 2;
            rv(1:6,1) = rv_ECI(1:6,1);
            
            rv(7:12,1) = rv_ECI(6*(i+1)-5:6*(i+1),1);
            mode = 3;
            [~, ~, maneuvers_out,~] = multisatellite_orbit_correction_3_impulse(rv, consts, spacecraft, formation, mode);

            Reconfiguration_matrix(i,j) = maneuvers_out;
        end
    end 
    close(f);
end