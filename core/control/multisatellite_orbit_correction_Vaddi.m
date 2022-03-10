function [t, rv_ECI, maneuvers_out] = multisatellite_orbit_correction_Vaddi(rv_ECI, consts, spacecraft, formation)

% Orbit correction with equinoctial orbital elements difference compensation
% The impulsive maneuvers are adopted from Vaddi Formation Establishment and Recon?guration Using Impulsive Control

T_local = 0;
t = 0;
maneuvers = [];

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);

rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6), formation.geometry, consts);

for i = 1:formation.N_sats
    rv_orb(:,i) = ECI2orb(rv_ECI(1:6,1), rv_ECI(i*6-5:i*6,1), consts);
    rv_ECI_required(:,i) = orb2ECI(rv_ECI(1:6,1), rv_orb_required(:,i), consts);
end

% Display Required and current orbital configuration
figure(1);
plot3(rv_orb_required(1,:), rv_orb_required(2,:), rv_orb_required(3,:), 'or');
hold on
plot3(rv_orb(1,:), rv_orb(2,:), rv_orb(3,:), 'ok');
xlabel('tangential,m');
ylabel('normal,m');
zlabel('radial,m');
axis equal;
legend('Required configuration', 'Current relative position');

corrections = vecnorm(rv_orb(1:3,:) - rv_orb_required(1:3,:))./vecnorm(rv_orb_required(1:3,:)) > formation.tracking_error;
corrections = [1:formation.N_sats; corrections];
disp('Satellites performing correction');
disp(num2str(corrections));

for i = 1:formation.N_sats
    if corrections(2,i) == 1
         maneuvers_new = calculate_maneuvers_Vaddi(rv_ECI(i*6-5:i*6,1),rv_ECI_required(:,i), consts, formation.geometry(:,i));
         maneuvers_new = [corrections(1,i); maneuvers_new];
         maneuvers = [maneuvers, maneuvers_new];
    end
end

maneuvers_out(1,:) = maneuvers(1,:);
maneuvers_out(2,:) = vecnorm(maneuvers(3:5,:));
maneuvers_out(3,:) = vecnorm(maneuvers(6:8,:));

maneuvers = maneuvers';
maneuvers = sortrows(maneuvers,2,'ascend');
maneuvers = maneuvers';
t_maneuvering = maneuvers(2,1);

%% Correction cycle
while size(maneuvers,2) ~= 0
    
    if t_maneuvering ~= T_local
        rv_init = rv_ECI(:,end);

        [t_vec_local, rv_ECI_local] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft),[T_local t_maneuvering], rv_init, options_precision);

        t = [t; t_vec_local];
        rv_ECI_local = rv_ECI_local';
        rv_ECI = [rv_ECI, rv_ECI_local];
        T_local = t(end);
    
    end
    
    % identifying what satellite is maneuvering
    active_sat = maneuvers(1,1);

    % applying impulse
    if ~isnan(maneuvers(3,1))
        
        % this means that the impulse is to be applied at theta critical
        rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(3:5,1)], consts);
        maneuvers(3:5,1) = NaN;

        active_sat_oe = rv2oe(rv_ECI(6*active_sat-5:6*active_sat,1), consts);
        
        maneuvers(2,1) = pi * sqrt(active_sat_oe(1)^3/consts.muEarth);
        maneuvers(2,1) = maneuvers(2,1) + T_local;
        
    else 
        
        rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(6:8,1)], consts);
        maneuvers(6:8,1) = NaN;
        disp(['satellite ', num2str(active_sat), ' has finished orbit correction!']);
        maneuvers(:,1) = [];
        
    end
    
    if size(maneuvers,2) ~= 0    
        maneuvers = maneuvers';
        maneuvers = sortrows(maneuvers,2,'ascend');
        maneuvers = maneuvers';

        t_maneuvering = maneuvers(2,1);
    end   
        
end

        % formation flying
%         [t_vec_ff, rv_ECI_ff] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft),[T_local T_local+consts.day2sec/12], rv_ECI(:,end), options_precision);
% 
%         t = [t; t_vec_ff];
%         rv_ECI_ff = rv_ECI_ff';
%         rv_ECI = [rv_ECI, rv_ECI_ff];
%         T_local = t(end);

%         for i = 1:formation.N_sats
%             rv_ECI_required_column(6*i-5:6*i,1) = rv_ECI_required(1:6,i);
%         end
        
%         [t_vec_ff_required, rv_ECI_ff_required] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft),[0 T_local], rv_ECI_required_column, options_precision);
%         rv_ECI_ff_required = rv_ECI_ff_required';
        
        
    for j = 1:size(rv_ECI,2)
        for i = 1:formation.N_sats
            rv_orb_final(:,i,j) = ECI2orb(rv_ECI(1:6,j), rv_ECI(6*i-5:6*i,j), consts);
        end
    end
%     
%     coe = rv2oe(rv_ECI(1:6,1), consts);
%     mean_motion = sqrt(consts.muEarth/ coe(1)^3);
%     for i = 1:length(t)
%         rv_orb_required_HCW(:,:,i) = rv_orb_HCW(formation.geometry, mean_motion, t(i));
%     end
%     
% %     for j = 1:size(rv_ECI_ff_required,2)
% %         for i = 1:formation.N_sats
% %             rv_orb_required_ode(:,i,j) = ECI2orb(rv_ECI_ff_required(1:6,j), rv_ECI_ff_required(6*i-5:6*i,j), consts);
% %         end
% %     end
% %     
% %     rv_orb_required_final(:,:,j) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,j), formation.geometry, consts);
%            
%     fig = figure;
%     xlabel('tangential,m');
%     ylabel('normal,m');
%     zlabel('radial,m');
%     xlim([-3000 3000]);
%     ylim([-3000 3000]);
%     zlim([-3000 3000]);
%     pbaspect([1 1 1]);
%     title('Orbital reference frame');
%     view(-45,-75);
%     hold on;
%     
%     for i = 1:size(rv_orb_final,3)/20
%         
%         a = plot3(rv_orb_required_HCW(1,:,i*20), rv_orb_required_HCW(2,:,i*20), rv_orb_required_HCW(3,:,i*20), 'or');
%         hold on;
%         b = plot3(rv_orb_final(1,:,i*20), rv_orb_final(2,:,i*20), rv_orb_final(3,:,i*20), 'ok');
%         legend('Required configuration', 'Current relative position');
%         
%         drawnow;
%         
% %         Take a Snapshot
%         movieVector(i) = getframe(fig);   %manually specify getframe region    
% 
%         delete(b);
%         delete(a)
%     end
%     plot3(rv_orb_final(1,:,end), rv_orb_final(2,:,end), rv_orb_final(3,:,end), 'ok');
%     plot3(rv_orb_required_HCW(1,:,end), rv_orb_required_HCW(2,:,end), rv_orb_required_HCW(3,:,end), 'or');
%     
%     myWriter = VideoWriter('Test','MPEG-4');   %create an .mp4 file
%     myWriter.FrameRate = 40;
% 
% %     Open the VideoWriter object, write the movie, and close the file
%     open(myWriter);
%     writeVideo(myWriter, movieVector);
%     close(myWriter);
% % 
% %     figure;
% %     subplot(1,2,1);
% %     plot3(rv_orb_required(1,:), rv_orb_required(2,:), rv_orb_required(3,:), 'or');
% %     hold on
% %     plot3(rv_orb(1,:), rv_orb(2,:), rv_orb(3,:), 'ok');
% %     xlabel('tangential,m');
% %     ylabel('normal,m');
% %     zlabel('radial,m');
% %     axis equal;
% %     legend('Required configuration', 'Current relative position');
% % 
% %     subplot(1,2,2);
% %     plot3(rv_orb_required_final(1,:), rv_orb_required_final(2,:), rv_orb_required_final(3,:), 'or');
% %     hold on;
% %     plot3(rv_orb_final(1,:), rv_orb_final(2,:), rv_orb_final(3,:), 'ok');
% %     xlabel('tangential,m');
% %     ylabel('normal,m');
% %     zlabel('radial,m');
% %     axis equal;
% %     legend('Required configuration', 'Current relative position');
% 

    for i = 1:size(rv_ECI,2)
        rv_orb_required(:,:,i) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,i), formation.geometry, consts);
    end
    
    figure;
    plot3(rv_orb_required(1,:,end), rv_orb_required(2,:,end), rv_orb_required(3,:,end), 'or');
    hold on;
    plot3(rv_orb_final(1,:,end), rv_orb_final(2,:,end), rv_orb_final(3,:,end), 'ok');
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    axis equal;
    legend('Required configuration', 'Current relative position');

end