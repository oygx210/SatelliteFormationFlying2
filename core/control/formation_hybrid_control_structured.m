function [t_vec,rv_ECI, maneuvers, post_correction, fuel_consumption, t_events, formation_state, HCW_constants_assigned, match_matrix] = formation_hybrid_control_structured(rv_ECI, demonstration, formation, spacecraft, consts)
 
    global T;
    t_vec =[T];
    % Event functions

    formation_state = [];
    t_events = [];
    maneuvers = []; % dV magnitudes for corrections    
    fuel_consumption = [];
    post_correction = [];
    for k = 1:length(demonstration)
        %% Reconfiguration
        disp('Regime: Reconfiguration');
        disp(['Orbital configuration ' num2str(k)]);
        [t_vec_reconf, rv_ECI_reconf, HCW_constants_assigned(:,:,k),match_matrix(:,:,k), maneuvers_reconfiguration, post_correction, fuel_consumption_reconfiguration] = reconfiguration(rv_ECI(:,end), demonstration{k,1}.HCW_constants, demonstration{k,1}.Cost_matrix_dV, formation, spacecraft, consts);
        post_correction = [post_correction, post_correction];
        t_vec = [t_vec; t_vec_reconf(2:end)];
        rv_ECI = [rv_ECI, rv_ECI_reconf(:,2:end)];
        
        formation_state = [formation_state;1];
        if k == 1
            t_events = [t_vec_reconf(1), t_vec_reconf(end)];
        else
            t_events = [t_events; [t_vec_reconf(1), t_vec_reconf(end)]];
        end
        formation.geometry = HCW_constants_assigned(:,:,k);
        maneuvers = [maneuvers, maneuvers_reconfiguration];
        fuel_consumption = [fuel_consumption, fuel_consumption_reconfiguration];
        formation.fuel_level = formation.fuel_level - fuel_consumption_reconfiguration(:,end);

        %% 2. Maintenance (up to image demonstration)

        disp('Regime: Maintenance');
        T_maintenance = demonstration{k,1}.demo_time(1);
        dt_rude_control = seconds(minutes(20));

        [t_vec_maintenance, rv_ECI_maintenance, maneuvers_maintenance] = maintenance(rv_ECI(:,end),T_maintenance, dt_rude_control, formation, spacecraft, consts);
        t_vec = [t_vec; t_vec_maintenance(2:end)];
        rv_ECI = [rv_ECI, rv_ECI_maintenance(:,2:end)];
        T = t_vec(end);
        
        fuel_level_maintenance = write_off_fuel(formation.fuel_level,[1:formation.N_active_sats],maneuvers_maintenance, spacecraft, consts);
        fuel_consumption_maintenance = formation.fuel_level - fuel_level_maintenance(:,end);
        formation.fuel_level = fuel_level_maintenance(:,end);

        maneuvers = [maneuvers, maneuvers_maintenance];
        fuel_consumption = [fuel_consumption, fuel_consumption_maintenance];
        t_events = [t_events;[t_events(end,2), T]]; 
        formation_state = [formation_state;2];
    
        %% Stage 3. Image demonstration

        disp(['Demonstration ', num2str(k), ', Regime: StandBy']);
        options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [t_out, rv_ECI_out] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [T:demonstration{k,1}.demo_time(2)], rv_ECI(:,end), options_precision);
        rv_ECI_out = rv_ECI_out'; 
        rv_ECI = [rv_ECI, rv_ECI_out(:,2:end)];
        t_vec = [t_vec; t_out(2:end)];    
        T = t_vec(end);
        formation_state = [formation_state; 3];
        t_events = [t_events; [t_out(1), t_out(end)]];                
        maneuvers = [maneuvers, zeros(formation.N_active_sats,1)];
        fuel_consumption = [fuel_consumption, zeros(formation.N_active_sats,1)];

        %% 4. Maintenance (up to reconfiguration)

        disp('Regime: Maintenance');
        T_maintenance = seconds(demonstration{k,1}.reconfiguration_time - formation.orbit_epoch);
        dt_rude_control = seconds(minutes(15));
        [t_vec_maintenance, rv_ECI_maintenance, maneuvers_maintenance] = maintenance(rv_ECI(:,end),T_maintenance, dt_rude_control, formation, spacecraft, consts);
        t_vec = [t_vec; t_vec_maintenance(2:end)];
        rv_ECI = [rv_ECI, rv_ECI_maintenance(:,2:end)];
        T = t_vec(end);
        
        fuel_level_maintenance = write_off_fuel(formation.fuel_level,[1:formation.N_active_sats],maneuvers_maintenance, spacecraft, consts);
        fuel_consumption_maintenance = formation.fuel_level - fuel_level_maintenance(:,end);
        formation.fuel_level = fuel_level_maintenance(:,end);

        maneuvers = [maneuvers, maneuvers_maintenance];
        fuel_consumption = [fuel_consumption, fuel_consumption_maintenance];
        t_events = [t_events;[t_events(end,2), T]]; 
        formation_state = [formation_state;2];
    end
    
    disp(['Formation Flying Mission with ', num2str(size(demonstration,1)), ' Image Demonstations Complete!']);

end   
