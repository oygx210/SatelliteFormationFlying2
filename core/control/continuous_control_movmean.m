function [t_vec, rv_ECI, formation_fuel_level] = continuous_control_movmean(rv_ECI, consts, spacecraft, formation, mode, MaxControlTime)

    options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
    global Rdiag
    control_cycle = seconds(minutes(10));
    global T
    t_vec = [T];
    convergence = 0;
    while convergence == 0 && t_vec(end) - T < MaxControlTime
        [t_out, rv_ECI_out] = ode45(@(t, rv) rhs_Formation_LQR(t, rv, consts, spacecraft, formation), [T:T+control_cycle], rv_ECI(:,end), options_precision);
        rv_ECI_out = rv_ECI_out';
        rv_ECI = [rv_ECI, rv_ECI_out];
        t_vec = [t_vec; t_out];
        T = t_vec(end);
        
        for j = 1:size(rv_ECI_out,2)
            rv_orb_required_final(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI_out(1:6,j), formation.geometry, consts);
            for i = 1:formation.N_sats
                rv_orb_final(:,i,j) = ECI2orb(rv_ECI_out(1:6,j), rv_ECI_out(6*i-5:6*i,j), consts);
                deltaR_vec(i,j) = vecnorm(rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j));
                deltaV_vec(i,j) = vecnorm(rv_orb_required_final(4:6,i,j) - rv_orb_final(4:6,i,j));
                
            end        
        end
        for i = 1:formation.N_sats
            C(i,:) = movmean(deltaR_vec(i,:), size(rv_ECI_out,2));
        end
        if sum(C(:,end) < 10) == formation.N_sats
            convergence = 1;
        end
        figure('Name','Continuous control', 'NumberTitle', 'off');
        for i = 1:formation.N_sats
            subplot(1,2,1);
            plot(t_out, deltaR_vec(i,:));
            hold on;
            xlabel('t');
            ylabel('\delta r');
            subplot(1,2,2);
            plot(t_out, C(i,:));
        end
            xlabel('t');
            ylabel('\delta r smoothed');
            legend(['Thrust = ' num2str(spacecraft.thrust) ' [N], Rdiag = ' num2str(Rdiag')]);
       
    end

        
%     end
%     
%     disp('LQR-based control converged in ');
%     
%     t_vec = t_out;    
%     rv_ECI = rv_ECI_out;
%     
%     
%     
%     
%     
%     if mode == 2
%         for j = 1:size(rv_ECI,2)
% 
%             for i = 1:formation.N_sats
%                 rv_orb_final(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
%             end
% 
%             rv_orb_required_final(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), formation.geometry, consts);
% 
%             for i = 1:formation.N_sats
% 
%                 rv_ECI_required_final(:,i,j) = orb2ECI(rv_ECI(1:6,j), rv_orb_required_final(:,i,j), consts);
% 
%                 deltaR_vec(:,i,j) = rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j);
%                 deltaR_magn(i,j) = vecnorm(rv_orb_required_final(1:3,i,j) - rv_orb_final(1:3,i,j));
%                 epsilon(i,j) = deltaR_magn(i,j) / vecnorm(rv_orb_required_final(1:3,i,j));
% 
%             end
%         end
% 
%         epsilon_mean = mean(epsilon(2:end,end));        
%         delta_r_mean = mean(deltaR_magn(:,end));
%         delta_r_sigma = std(deltaR_magn(:,end));
%         sat_to_demonstrate = 25;
% 
% 
%         figure('Name', 'Relative position error during correction', 'NumberTitle','off');
%         for i = 1:formation.N_active_sats
%             deltaR(i,:) = vecnorm(squeeze(deltaR_vec(:,i+1,:)));
%             plot(t_vec/60, deltaR(i,:));
%             grid on;
%             hold on;
%             xlabel('t, min');
%             ylabel('\deltar, meters');
%             legend('\delta r');
%         end
%     end    
formation_fuel_level = 0;
end