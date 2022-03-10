function [t_vec, rv_ECI, maneuvers_maintenance] = maintenance(rv_ECI, T_maintenance, dt_rude_control, formation, spacecraft, consts)
    
    global T
    t_vec = [T];
    maneuvers_maintenance = zeros(formation.N_active_sats,1);
    tracking_formation_quality = @(t, rv) formation_quality(rv_ECI(:,end), consts, spacecraft, formation);
    options_quality = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', tracking_formation_quality);
    T_rude_control = T_maintenance - dt_rude_control;
    
    while T < T_maintenance

        quality = formation_quality(rv_ECI(:,end), consts, spacecraft, formation);

        while quality == 1 && T < T_rude_control
            [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [T + [0:10]], rv_ECI(:,end), options_quality);
            rv_ECI_FF = rv_ECI_FF';
            rv_ECI = [rv_ECI, rv_ECI_FF(:,2:end)];
            t_vec = [t_vec; t_vec_FF(2:end)];
            T = t_vec(end);
            quality = formation_quality(rv_ECI(:,end), consts, spacecraft, formation);        
        end

        if T < T_maintenance                   
            mode = 1;
            [t_vec_m, rv_ECI_m, continuous_maneuvers_dV, ~] = continuous_control(rv_ECI(:,end), T_maintenance, T_rude_control, consts, spacecraft, formation, mode);
            maneuvers_maintenance = maneuvers_maintenance + continuous_maneuvers_dV;                                
            rv_ECI = [rv_ECI, rv_ECI_m(:,2:end)];
            t_vec = [t_vec; t_vec_m(2:end)];        
            T = t_vec(end);
            quality = formation_quality(rv_ECI(:,end), consts, spacecraft, formation);                                     
        end                
    end

end