function [t, rv_ECI, maneuvers_out, formation_fuel] = multisatellite_orbit_correction_4_impulse(rv_ECI, consts, spacecraft, formation, mode)

% 1. Correcting \delta sma and \delta ecc
% 2. Correction difference in equinoctial oe expcept for sma

    T_local = 0;
    t = 0;
    maneuvers = [];

    % mode: 1 - operational, 2 - debug, 3 - reconfiguration matrix computation

    if mode == 2
        sat_to_demonstrate = 2;
    end
    options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);

    rv_orb_required_initial = get_rv_from_analytic_HCW_solution(rv_ECI(1:6), formation.geometry, consts);

    for i = 1:formation.N_sats
        rv_orb(:,i) = ECI2orb(rv_ECI(1:6,1), rv_ECI(i*6-5:i*6,1), consts);
        rv_ECI_required(:,i) = orb2ECI(rv_ECI(1:6,1), rv_orb_required_initial(:,i), consts);
    end

    if mode == 2
        %  Displaying Required and current orbital configuration
        if formation.N_sats == 2 
            figure('Name','Pre-correction orbital configuration','NumberTitle','off');
            plot3(rv_orb_required_initial(1,2), rv_orb_required_initial(2,2), rv_orb_required_initial(3,2), 'or');
            hold on
            plot3(rv_orb(1,2), rv_orb(2,2), rv_orb(3,2), 'ok');
            xlabel('tangential,m');
            ylabel('normal,m');
            zlabel('radial,m');
            axis equal;
            legend('Required configuration', 'Current relative position');
        else
            figure('Name','Pre-correction orbital configuration','NumberTitle','off');
            plot3(rv_orb_required_initial(1,:), rv_orb_required_initial(2,:), rv_orb_required_initial(3,:), 'or');
            hold on;
            plot3(rv_orb(1,:), rv_orb(2,:), rv_orb(3,:), 'ok');
            xlabel('tangential,m');
            ylabel('normal,m');
            zlabel('radial,m');
            axis equal;
            legend('Required configuration', 'Current relative position');
        end
    end

corrections = vecnorm(rv_orb(1:3,:) - rv_orb_required_initial(1:3,:))./vecnorm(rv_orb_required_initial(1:3,:)) > formation.tracking_error;
corrections = [1:formation.N_sats; corrections];
disp('Satellites performing correction');
disp(num2str(corrections));

    %% Control
    for i = 1:formation.N_sats
        if corrections(2,i) == 1

             [t_span, dV_perigee, dV_apogee] = calculate_sma_ecc_correction(rv_ECI(i*6-5:i*6,1), rv_ECI_required(:,i), consts);

             maneuvers_new = [corrections(1,i); 0; t_span; dV_perigee ;dV_apogee; zeros(6,1)];
             maneuvers = [maneuvers, maneuvers_new];
        end
    end
 
    maneuvers = maneuvers';
    maneuvers = sortrows(maneuvers,3,'ascend');
    maneuvers = maneuvers';
    t_maneuvering = maneuvers(3,1);

    %% Correction cycle
    while sum(maneuvers(2,:) == 4) < sum(corrections(2,:))

        if t_maneuvering ~= T_local
            
            rv_init = rv_ECI(:,end);
            [t_vec_local, rv_ECI_local] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[T_local t_maneuvering], rv_init, options_precision);

            t = [t; t_vec_local];
            rv_ECI_local = rv_ECI_local';
            rv_ECI = [rv_ECI, rv_ECI_local];
            T_local = t(end);

        end

        % identifying what satellite is maneuvering
        active_sat = maneuvers(1,1);

        % applying impulse
        switch maneuvers(2,1)
            case 0 

             active_sat_oe_before_maneuvering = rv2oe(rv_ECI(6*active_sat-5:6*active_sat,end), consts);
             rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(4:6,1)], consts);
    %          rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; 0;0;0], consts);
             active_sat_oe = rv2oe(rv_ECI(6*active_sat-5:6*active_sat,end), consts);

             if active_sat_oe(6) < 2*pi && active_sat_oe(6) > pi
                 t_span = (3*pi - active_sat_oe(6)) * sqrt(active_sat_oe(1)^3/consts.muEarth);
             else 
                 t_span = (pi - active_sat_oe(6)) * sqrt(active_sat_oe(1)^3/consts.muEarth);
             end

             maneuvers(3,1) = T_local + t_span;
             maneuvers(2,1) = 1;

            case 1

            % Impulse located at apogee or perigee and aimed to adjust sma
            rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(7:9,1)], consts);
            %rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; 0;0;0], consts);

            rv_ECI_current = rv_ECI(6*active_sat-5:6*active_sat,end);
            rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,end), formation.geometry(:,active_sat), consts);
            rv_ECI_required = orb2ECI(rv_ECI(1:6,end), rv_orb_required, consts);

            [t_span, dV1_Vaddi, dV2_Vaddi] = calculate_maneuvers_Vaddi(rv_ECI_current, rv_ECI_required, consts);

            maneuvers(3,1) = T_local + t_span;
            maneuvers(10:12,1) = dV1_Vaddi;
            maneuvers(13:15,1) = dV2_Vaddi;

            maneuvers(2,1) = 2;

            case 2

            rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(10:12,1)], consts);

            active_sat_oe = rv2oe(rv_ECI(6*active_sat-5:6*active_sat,end), consts);
            t_span = pi * sqrt(active_sat_oe(1)^3/consts.muEarth);
            maneuvers(3,1) = t_span + T_local;

            maneuvers(2,1) = 3;

            case 3

            rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(13:15,1)], consts);
            maneuvers(2,1) = 4;
            maneuvers(3,1) = NaN;

            disp(['satellite ', num2str(active_sat), ' has finished orbit correction!']);

        end

        maneuvers = maneuvers';
        maneuvers = sortrows(maneuvers,3,'ascend');
        maneuvers = maneuvers';

        t_maneuvering = maneuvers(3,1);
        
    end
%% Calculating spent fuel    
    maneuvers_out(1,:) = maneuvers(1,:);
    maneuvers_out(2,:) = vecnorm(maneuvers(4:6,:));
    maneuvers_out(3,:) = vecnorm(maneuvers(7:9,:));
    maneuvers_out(4,:) = vecnorm(maneuvers(10:12,:));
    maneuvers_out(5,:) = vecnorm(maneuvers(13:15,:));
    maneuvers_out = maneuvers_out';
    maneuvers_out = sortrows(maneuvers_out,1,'ascend');
    maneuvers_out = maneuvers_out';
    max_dV(1) = max(maneuvers_out(2,:));
    max_dV(2) = max(maneuvers_out(3,:));
    max_dV(3) = max(maneuvers_out(4,:));
    max_dV(4) = max(maneuvers_out(5,:));

    if mode ~= 3
        formation.fuel_level = write_off_fuel(formation.fuel_level, maneuvers_out, spacecraft, consts);
    end

    formation_fuel = formation.fuel_level;    
%% Visualization    

    if mode == 2

    t_finish = T_local;
    t_finish_step = length(t);

    figure('Name', 'Fuel level','NumberTitle','off');            
    for i = 1:size(formation.fuel_level,2)                
        plot(formation.fuel_level(:,i)*1000);
        xlabel('satellite');
        ylabel('Fuel level, gramms');
        hold on;
    end

    [t_vec_local, rv_ECI_local] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), T_local:T_local + 6.0293e+03, rv_ECI(:,end), options_precision);
    rv_ECI_local = rv_ECI_local';
    rv_ECI = [rv_ECI, rv_ECI_local];
    t = [t; t_vec_local];
    T_local = t(end);

    for j = 1:size(rv_ECI,2)

        for i = 1:formation.N_sats
            rv_orb_final(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
        end

        rv_orb_required_final(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), formation.geometry, consts);

        for i = 1:formation.N_sats

            rv_ECI_required_final(:,i,j) = orb2ECI(rv_ECI(1:6,j), rv_orb_required_final(:,i,j), consts);

            deltaR_vec(:,i,j) = rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j);
            deltaR_magn(i,j) = vecnorm(rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j));
            epsilon(i,j) = deltaR_magn(i,j) / vecnorm(rv_orb_required_final(1:3,i,j));

        end
    end

    epsilon_mean = mean(epsilon(2:end,end));        
    delta_r_mean = mean(deltaR_magn(:,end));
    delta_r_sigma = std(deltaR_magn(:,end));

    figure('Name','Orbital configuration after correction','NumberTitle','off');
    for i = 1:formation.N_active_sats
        r_relative = squeeze(rv_orb_final(1:3,i+1,(t_finish_step:end)));
        if i+1 == sat_to_demonstrate
            plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
            hold on;
        else                      
            plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
            hold on;
        end
        plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.1); 

        hold on;
    end

    xlabel('x, meters');
    ylabel('y, meters');
    zlabel('z, meters');
    legend('first orbital configuration', 'relative trajectories', 'q');

    deltaR_x = squeeze(deltaR_vec(1,sat_to_demonstrate,:));
    deltaR_y = squeeze(deltaR_vec(2,sat_to_demonstrate,:));
    deltaR_z = squeeze(deltaR_vec(3,sat_to_demonstrate,:));
    deltaR = vecnorm(squeeze(deltaR_vec(:,sat_to_demonstrate,:)));
    
    figure('Name', 'Relative position error during correction', 'NumberTitle','off');
    plot(t/60,  deltaR_x); 
    hold on;
    plot(t/60,  deltaR_y); 
    hold on;
    plot(t/60,  deltaR_z); 
    grid on;
    plot(t/60, deltaR);
    xlabel('t, min');
    ylabel('\deltar, meters');
    legend('\deltar_x', '\deltar_y', '\deltar_z', '|\delta r|');

%     axes('position',[.15 .50 .25 .25])
%     box on % put box around new pair of axes
% 
%     plot(t(3500:6000)/60,  deltaR_x(3500:6000)); 
%     hold on;
%     plot(t(3500:6000)/60,  deltaR_y(3500:6000)); 
%     hold on;
%     plot(t(3500:6000)/60,  deltaR_z(3500:6000)); 
%     grid on;
%     axis tight

    end        

end