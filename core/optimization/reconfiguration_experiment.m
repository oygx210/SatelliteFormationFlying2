function reconfiguration_experiment(rv_ECI, consts, spacecraft, formation, reconfiguration_matrix_dV)

    T_local = 0;
    rv = rv_ECI(:,end);
    formation.geometry = formation.geometry1;
    quality = formation_quality(T_local, rv, consts, spacecraft, formation);

    tic;
    reconfiguration_matrix_dV1 = get_reconfiguration_matrix(rv, consts, spacecraft, formation);
    toc;              

    reconfiguration_matrix_dV2 = reconfiguration_matrix_dV;

    for i = 1:10

        reconfiguration_matrix_fuel1 = zeros(formation.N_active_sats);            
        for j = 1:formation.N_active_sats               
            spacecraft_wet_mass_updated = (spacecraft.dry_mass + formation.fuel_level(j,end))*exp(-reconfiguration_matrix_dV1(j,:)/spacecraft.thruster_Isp/consts.g);
            reconfiguration_matrix_fuel1(j,:) = (spacecraft.dry_mass + formation.fuel_level(j,end))*ones(1, formation.N_active_sats) - spacecraft_wet_mass_updated;
        end

        [~, Fleft] = maneuverAssignment(reconfiguration_matrix_fuel1, formation.fuel_level(:,end));    
        formation.fuel_level = [formation.fuel_level, Fleft];

        reconfiguration_matrix_fuel2 = zeros(formation.N_active_sats);            
        for j = 1:formation.N_active_sats               
            spacecraft_wet_mass_updated = (spacecraft.dry_mass + formation.fuel_level(j,end))*exp(-reconfiguration_matrix_dV2(j,:)/spacecraft.thruster_Isp/consts.g);
            reconfiguration_matrix_fuel2(j,:) = (spacecraft.dry_mass + formation.fuel_level(j,end))*ones(1, formation.N_active_sats) - spacecraft_wet_mass_updated;
        end
        
        [~, Fleft] = maneuverAssignment(reconfiguration_matrix_fuel2, formation.fuel_level(:,end));    
        formation.fuel_level = [formation.fuel_level, Fleft];
    end

    for i = 1:size(formation.fuel_level,2)-1
        sigma_reconf(i) = std(formation.fuel_level(:,1+i));
        plot(formation.fuel_level(:,1+i)*1000);
        xlabel('satellite number');
        ylabel('Fuel level, gramms');
        hold on;
    end
    figure;
    plot(sigma_reconf);
    xlabel('reconfiguration number');
    ylabel('fuel distribution \sigma');
    
end