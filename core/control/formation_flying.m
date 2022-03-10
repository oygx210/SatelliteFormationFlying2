function [t_vec,rv_ECI, maneuvers, T_correction_start, T_correction_end] = formation_flying(rv_ECI, simulation_time, consts, spacecraft, formation)

    T_local = 0;
    t_vec =[];
    T_correction_start = [];
    T_correction_end = [];
    maneuvers = []; % dV magnitudes for corrections
    
    % Event functions
    tracking_formation_quality = @(t, rv) formation_quality(t, rv, consts, spacecraft, formation);
    options_quality = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', tracking_formation_quality);

    rv = rv_ECI(:,end);
    quality = formation_quality(T_local, rv, consts, spacecraft, formation);

    %% Deployment
    while quality == 0 
        disp('Deployment');
        mode = 1;
        T_correction_start = [T_correction_start, T_local];
        [t_vec_m, rv_ECI_m, maneuvers_dV, formation.fuel_level] = multisatellite_orbit_correction_4_impulse(rv, consts, spacecraft, formation, mode);
        maneuvers = [maneuvers, maneuvers_dV];
        t_vec = [t_vec; t_vec_m];
        rv_ECI = [rv_ECI, rv_ECI_m];
        T_local = t_vec(end);

        T_correction_end = [T_correction_end, T_local];        
        quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
        rv = rv_ECI(:,end);
    end
    
%% Normal operation depending on mode: 0 - no reconfiguration, 1 - reconfiguration, 2 - reconfiguration with maximin optimization and collision check

    switch formation.reconfiguration_flag
        case 0
            %% Formation has same geometry during the mission - no reconfiguration
            while T_local < simulation_time

                disp('Maintenance');
                rv = rv_ECI(:,end);
                [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [T_local simulation_time], rv, options_quality);
                rv_ECI_FF = rv_ECI_FF';
                rv_ECI = [rv_ECI, rv_ECI_FF];
                t_vec = [t_vec; t_vec_FF];

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                T_local = t_vec(end);

                if T_local < simulation_time && quality == 0

                    while quality == 0
                    disp(['Correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);
                    mode = 1;

                    T_correction_start = [T_correction_start, T_local];                    
%                     [t_vec_maneuvering, rv_maneuvering, maneuvers_dV, formation.fuel_level] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);
                    mode = 1;
                    [t_vec_maneuvering, rv_maneuvering, maneuvers_dV, formation.fuel_level] = multisatellite_orbit_correction_4_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);

                    rv_ECI = [rv_ECI, rv_maneuvering];
                    t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
                    rv = rv_ECI(:,end);
                    T_local = t_vec(end);
                    T_correction_end = [T_correction_end, T_local];        
                    
                    disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

                    maneuvers = [maneuvers, maneuvers_dV];

                    quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                    
                    if T_local > simulation_time
                        return;
                    end
                    
                    end
                    

                end

                if T_local >= simulation_time
                    disp('Mission with no reconfiguration complete!');
                end
            end
            
        case 1
            %% Control with non optimized reconfiguration
            while T_local < formation.reconfiguration_time

                disp('Maintenance');
                rv = rv_ECI(:,end);
%                 [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), [T_local formation.reconfiguration_time], rv, options_quality);
                [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) central_gravity_Formation(t, rv, consts), [T_local formation.reconfiguration_time], rv, options_quality);
                rv_ECI_FF = rv_ECI_FF';
                rv_ECI = [rv_ECI, rv_ECI_FF];
                t_vec = [t_vec; t_vec_FF];

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                T_local = t_vec(end);

                if T_local < formation.reconfiguration_time && quality == 0

                    while quality == 0

                    disp(['Correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);
                    % Disp sats whose positions are to be corrected 

                    [t_vec_maneuvering, rv_maneuvering] = multisatellite_orbit_correction_4_impulse(rv_ECI(:,end), consts, spacecraft, formation);

                    rv_ECI = [rv_ECI, rv_maneuvering];
                    t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
                    rv = rv_ECI(:,end);
                    T_local = t_vec(end);
                    disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

                    maneuvers = [maneuvers, maneuvers_dV];

                    quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                    end

                end                
            end
            
            disp('Reconfiguration time has come');
            formation.geometry = formation.geometry2;
            rv = rv_ECI(:,end);
            quality = formation_quality(T_local, rv, consts, spacecraft, formation);

            while quality == 0 
                disp('Orbit correction');
                [t_vec_m, rv_ECI_m, maneuvers_dV] = multisatellite_orbit_correction_4_impulse(rv, consts, spacecraft, formation);
                maneuvers = [maneuvers, maneuvers_dV];

                t_vec = [t_vec; t_vec_m];
                rv_ECI = [rv_ECI, rv_ECI_m];
                T_local = t_vec(end);

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                rv = rv_ECI(:,end);
            end
            
            while T_local < simulation_time

                disp('Maintenance');
                rv = rv_ECI(:,end);
                [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), [T_local simulation_time], rv, options_quality);
                rv_ECI_FF = rv_ECI_FF';
                rv_ECI = [rv_ECI, rv_ECI_FF];
                t_vec = [t_vec; t_vec_FF];

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                T_local = t_vec(end);

                if T_local < simulation_time && quality == 0

                    while quality == 0

                    disp(['Correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);
                    % Disp sats whose positions are to be corrected 

                    [t_vec_maneuvering, rv_maneuvering] = multisatellite_orbit_correction_4_impulse(rv_ECI(:,end), consts, spacecraft, formation);

                    rv_ECI = [rv_ECI, rv_maneuvering];
                    t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
                    rv = rv_ECI(:,end);
                    T_local = t_vec(end);
                    disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

                    maneuvers = [maneuvers, maneuvers_dV];

                    quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                    end

                end                
            end
            disp('Mission with non optimized reconfiguration complete!');  
            
        case 2
            %% Reconfiguration with fuel consumption optimization and collision avoidance
            while T_local < formation.reconfiguration_time

                disp('Maintenance');
                rv = rv_ECI(:,end);
                [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [T_local formation.reconfiguration_time], rv, options_quality);
                rv_ECI_FF = rv_ECI_FF';
                rv_ECI = [rv_ECI, rv_ECI_FF];
                t_vec = [t_vec; t_vec_FF];

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                T_local = t_vec(end);

                if T_local < formation.reconfiguration_time && quality == 0

                    while quality == 0

                    disp(['Correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);

                    mode = 1;
                    % add 4 impulse maneuvers
                    [t_vec_maneuvering, rv_maneuvering, maneuvers_out, fuel_level] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation, mode);

                    rv_ECI = [rv_ECI, rv_maneuvering];
                    t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
                    rv = rv_ECI(:,end);
                    T_local = t_vec(end);
                    disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

                    maneuvers = [maneuvers, maneuvers_dV];

                    quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                    end

                end                
            end
            %% Optimization 
            disp('Reconfiguration time has come');
            rv = rv_ECI(:,end);
            formation.geometry = formation.geometry2;
            quality = formation_quality(T_local, rv, consts, spacecraft, formation);
            
            tic;
            reconfiguration_matrix_dV = get_reconfiguration_matrix(rv, consts, spacecraft, formation);
            toc;              
            reconfiguration_matrix_fuel = zeros(formation.N_active_sats);
            
            for i = 1:formation.N_active_sats               
                spacecraft_wet_mass_updated = (spacecraft.dry_mass + formation.fuel_level(i,end))*exp(-reconfiguration_matrix_dV(i,:)/spacecraft.thruster_Isp/consts.g);
                reconfiguration_matrix_fuel(i,:) = (spacecraft.dry_mass + formation.fuel_level(i,end))*ones(1, formation.N_active_sats) - spacecraft_wet_mass_updated;
            end
            
            [matchMatrix, ~] = maneuverAssignment(reconfiguration_matrix_fuel, formation.fuel_level(:,end));    
             
            for i = 1:(formation.N_sats-1)
                formation.geometry(:,matchMatrix(i,1)+1) = formation.geometry2(:,matchMatrix(i,2)+1);
            end
            mode = 1;
            [t_vec_reconf, rv_ECI_reconf, maneuvers_dV] = multisatellite_orbit_correction_3_impulse(rv, consts, spacecraft, formation, mode);
            [collision_flag, ISD_min] = collision_check(t_vec_reconf, rv_ECI_reconf, formation);
                                               
            maneuvers_dV = [maneuvers, maneuvers_dV];
            rv_ECI = [rv_ECI, rv_ECI_reconf];
            t_vec = [t_vec; t_vec_reconf + t_vec(end)];
            T_local = t_vec(end);
            
            experimental_mode = 1;
            if experimental_mode == 1 
                reconfiguration_experiment(rv_ECI,consts, spacecraft, formation, reconfiguration_matrix_dV)
            end
                        
            while T_local < simulation_time

                disp('Maintenance');
                rv = rv_ECI(:,end);
%                 [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), [T_local simulation_time], rv, options_quality);
                [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) central_gravity_Formation(t, rv, consts), [T_local simulation_time], rv, options_quality);
                rv_ECI_FF = rv_ECI_FF';
                rv_ECI = [rv_ECI, rv_ECI_FF];
                t_vec = [t_vec; t_vec_FF];

                quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                T_local = t_vec(end);

                if T_local < simulation_time && quality == 0

                    while quality == 0

                    disp(['Correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);
                    % Disp sats whose positions are to be corrected 

                    [t_vec_maneuvering, rv_maneuvering] = multisatellite_orbit_correction_3_impulse(rv_ECI(:,end), consts, spacecraft, formation);

                    rv_ECI = [rv_ECI, rv_maneuvering];
                    t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
                    rv = rv_ECI(:,end);
                    T_local = t_vec(end);
                    disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

                    maneuvers = [maneuvers, maneuvers_dV];

                    quality = formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
                    end

                end                
            end
            disp('Mission with reconfiguration complete!');            
    end   
end