function [t_vec, rv_ECI, maneuvers, formation_fuel_level] = continuous_control(rv_ECI, T_maintenance, T_rude_control, consts, spacecraft, formation, mode)
    
    % mode = 1 - continuous control for relative orbits maintenance
    % mode = 2 - same as 1 + visualization on
    
    if mode == 2
        disp('    Continuous control');
    end
    global T;
    global K;
    maneuvers = zeros(formation.N_active_sats,1);
    options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
    tracking_formation_quality = @(t, rv) formation_quality(rv_ECI(:,end), consts, spacecraft, formation);
    options_quality = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', tracking_formation_quality);

    t_vec = T;
    
    convergence = 0;
    t_span_control = 0:5;
    e_r_max = 1;
    e_v_max = 0.01;
                
    rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,end), formation.geometry, consts);

    for i = 1:formation.N_sats
        rv_orb(:,i) = ECI2orb(rv_ECI(1:6,end), rv_ECI(6*i-5:6*i,end), consts);
    end
    if T < T_rude_control
        maneuvering_satellites = vecnorm(rv_orb(1:3,:) - rv_orb_required(1:3,:)) > formation.tracking_error;
    else
        maneuvering_satellites = ones(formation.N_sats,1);
    end
    
    if size(maneuvering_satellites,1) ~= formation.N_sats
        maneuvering_satellites = maneuvering_satellites';
    end
    
    while convergence ~= 1 || T < T_maintenance
        if T > T_rude_control
            maneuvering_satellites = [0;ones(formation.N_active_sats,1)];
        end
        e = rv_orb - rv_orb_required;
        e_r_mag = vecnorm(e(1:3,:));
        e_v_mag = vecnorm(e(4:6,:));
        u_switch = e_r_mag(:) > e_r_max | e_v_mag(:) > e_v_max;
        if size(u_switch,1) ~= formation.N_sats
            u_switch = u_switch';
        end
        try
            u_switch = u_switch.*maneuvering_satellites;
        catch
            disp(num2str(u_switch));
            disp(num2str(maneuvering_satellites));

        end
        z_orb = rv_ECI(1:3,end) / norm(rv_ECI(1:3,end));
        y_orb = cross(rv_ECI(1:3,end), rv_ECI(4:6,end));
        y_orb = y_orb(1:3) / norm(y_orb(1:3));
        x_orb = cross(y_orb, z_orb);
        orb2ECI_matrix = [x_orb y_orb z_orb];
             
        for i = 1:formation.N_sats
            u_orb = K*e(:,i)*u_switch(i);
            if norm(u_orb) > spacecraft.u_max
                u_orb = u_orb./norm(u_orb)*spacecraft.u_max;
            end   
            u_ECI(:,i) = orb2ECI_matrix*u_orb;
        end
    
        u_control = [zeros(3,formation.N_sats); u_ECI];
        [t_out, rv_ECI_out] = ode45(@(t, rv) rhs_Formation_inertial_LQR_input(t, rv, consts, spacecraft, u_control), [T + t_span_control], rv_ECI(:,end), options_precision);
        rv_ECI_out = rv_ECI_out'; 
        rv_ECI = [rv_ECI, rv_ECI_out(:,2:end)];
        t_vec = [t_vec; t_out(2:end)];    
        T = t_vec(end);
        
        maneuvers_new = vecnorm(u_ECI(:,2:end));
        maneuvers = maneuvers + maneuvers_new';
        
        rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,end), formation.geometry, consts);
       
        for i = 1:formation.N_sats
            rv_orb(:,i) = ECI2orb(rv_ECI(1:6,end), rv_ECI(6*i-5:6*i,end), consts);
        end
        e = rv_orb - rv_orb_required;
        e_r_mag = vecnorm(e(1:3,:));
        e_v_mag = vecnorm(e(4:6,:));
        u_switch = e_r_mag(:) > e_r_max | e_v_mag(:) > e_v_max;
        if size(u_switch,1) ~= formation.N_sats
            u_switch = u_switch';
        end
        convergence1 = u_switch.*maneuvering_satellites;
        convergence2 = vecnorm(rv_orb(1:3,:) - rv_orb_required(1:3,:)) > formation.tracking_error;

        if size(convergence2,1) ~= formation.N_sats
            convergence2 = convergence2';
        end
        if T < T_rude_control
            for i = 1:formation.N_sats
                if u_switch(i) == 0
                    maneuvering_satellites(i) = 0;
                end
            end
        end
        
        if T < T_rude_control
        if (sum(convergence1) == 0 & sum(convergence2) == 0)
            convergence = 1;
        else
            if sum(convergence2 > maneuvering_satellites) > 0
                logical = convergence2 > maneuvering_satellites;
                maneuvering_satellites = maneuvering_satellites + logical;
            for i = 1:length(maneuvering_satellites)
                if maneuvering_satellites(i) ~= 0
                    maneuvering_satellites(i) = 1;
                end
            end  
            end
        end
        else
            if (sum(convergence1) == 0 & sum(convergence2) == 0)
                convergence = 1;
            end
        end
    end  
    
    maneuvering_sats = 1:formation.N_active_sats;
    formation_fuel_level = write_off_fuel(formation.fuel_level, maneuvering_sats, maneuvers, spacecraft, consts);
    
    if mode == 2
        
        omega = sqrt(consts.muEarth/vecnorm(rv_ECI(1:3,1))^3);

        for j = 1:size(rv_ECI,2)

            for i = 1:formation.N_sats
                rv_orb_final(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
            end

            rv_orb_required_final(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), formation.geometry, consts);

            for i = 1:formation.N_sats

                rv_ECI_required_final(:,i,j) = orb2ECI(rv_ECI(1:6,j), rv_orb_required_final(:,i,j), consts);

                deltaR(i,j) = vecnorm(rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j));
                deltaV(i,j) = vecnorm(rv_orb_required_final(4:6,i,j) - rv_orb_final(4:6,i,j));
                c1(i,j) = rv_orb_final(4,i,j)/omega + 2*rv_orb_final(3,i,j);

            end
        end
                
        figure('Name', 'Errors during correction', 'NumberTitle','off');
        for i = 1:formation.N_active_sats
%         for i = 1
            subplot(1,3,1);
            plot(t_vec/60, deltaR(i+1,:));
            grid on;
            hold on;
            xlabel('t, min');
            ylabel('\deltar, meters');
            legend('\delta r');
            subplot(1,3,2);
            plot(t_vec/60, deltaV(i+1,:));
            grid on;
            hold on;
            xlabel('t, min');
            ylabel('\deltav, m/s');
            legend('\delta v');
            subplot(1,3,3);
            plot(t_vec/60, c1(i+1,:));
            grid on;
            hold on;
            xlabel('t, min');
            ylabel('drift constant c1');
            legend('c1');
        end
        
    end    
end