function [t, rv_ECI, maneuvers_out] = multisatellite_orbit_correction(rv_ECI, consts, spacecraft, formation)

% check what rv_ECI is used to 

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
         maneuvers_new = calculate_maneuvers(rv_ECI(i*6-5:i*6,1),rv_ECI_required(:,i), consts);
         maneuvers_new = [corrections(1,i); maneuvers_new];
         maneuvers = [maneuvers, maneuvers_new];
    end
end

maneuvers_out(1,:) = maneuvers(1,:);
maneuvers_out(2,:) = vecnorm(maneuvers(3:5,:));
maneuvers_out(3,:) = vecnorm(maneuvers(6:8,:));
maneuvers_out(4,:) = vecnorm(maneuvers(9:11,:));

maneuvers = maneuvers';
maneuvers = sortrows(maneuvers,2,'ascend');
maneuvers = maneuvers';
t_maneuvering = maneuvers(2,1);

%% Correction cycle
for i = 1:size(maneuvers,2)*3
    
    rv_init = rv_ECI(:,end);

    [t_vec_local, rv_ECI_local] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft),[T_local t_maneuvering], rv_init, options_precision);

    t = [t; t_vec_local];
    rv_ECI_local = rv_ECI_local';
    rv_ECI = [rv_ECI, rv_ECI_local];
    T_local = t(end);
    
    % identifying what satellite is maneuvering
    active_sat = maneuvers(1,1);

    % applying impulse
    if ~isnan(maneuvers(3,1))
        
        % this means that the impulse is to be applied at theta critical
        rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(3:5,1)], consts);
        maneuvers(3:5,1) = NaN;
        active_sat_oe = vecRV2OE(rv_ECI(6*active_sat-5:6*active_sat,1), consts);
        active_sat_moe = rv2moe(rv_ECI(6*active_sat-5:6*active_sat,1), consts);
        
        maneuvers(2,1) = (2*pi - active_sat_oe(6)) * sqrt(active_sat_moe(1)^3/consts.muEarth);
        maneuvers(2,1) = maneuvers(2,1) + T_local;
        
    elseif ~isnan(maneuvers(6,1))

        % this means that the impulse is to be pericenter
        rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(6:8,1)], consts);
        maneuvers(6:8,1) = NaN;
        active_sat_oe = vecRV2OE(rv_ECI(6*active_sat-5:6*active_sat,1), consts);
        active_sat_moe = rv2moe(rv_ECI(6*active_sat-5:6*active_sat,1), consts);
 
        if active_sat_oe(6) > pi
            maneuvers(2,1) = (3*pi - active_sat_oe(6)) * sqrt(active_sat_moe(1)^3/consts.muEarth);
        else
            maneuvers(2,1) = (pi - active_sat_oe(6)) * sqrt(active_sat_moe(1)^3/consts.muEarth);
        end
        
        maneuvers(2,1) = maneuvers(2,1) + T_local;

    else 
        
        rv_ECI(6*active_sat-5:6*active_sat,end) = orb2ECI(rv_ECI(6*active_sat-5:6*active_sat,end), [0; 0; 0; maneuvers(9:11,1)], consts);
        maneuvers(9:11,1) = NaN;
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
    
    sat2.rv_orb_required = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,:), formation.geometry(:,2), consts);
    for i = 1:size(rv_ECI,2)
        sat2.rv_ECI_required(:,i) = orb2ECI(rv_ECI(1:6,i), sat2.rv_orb_required, consts);
    end
    
%     sat2.moe_required = rv2moe(sat2.rv_ECI_required, consts);
%     sat2.moe = rv2moe(rv_ECI(7:12,:), consts);
%     sat2.delta_moe = sat2.moe_required - sat2.moe;
    
    plot_mean_oe_diff(t, sat2.rv_ECI_required, rv_ECI(7:12,:), consts, 'm');
        
  
    for i = 1:formation.N_sats
        rv_orb_final(:,i) = ECI2orb(rv_ECI(1:6,end), rv_ECI(6*i-5:6*i,end), consts);
    end
    
    % we should see what happened to deltas in oe
    % the problem is what values should we compare? we can recover target
    % trajectory based on leader satellite and formation geometry constants
    
    rv_orb_required_final = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,end), formation.geometry, consts);
  
    figure;
    subplot(1,2,1);
    plot3(rv_orb_required(1,:), rv_orb_required(2,:), rv_orb_required(3,:), 'or');
    hold on
    plot3(rv_orb(1,:), rv_orb(2,:), rv_orb(3,:), 'ok');
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    axis equal;
    legend('Required configuration', 'Current relative position');

    subplot(1,2,2);
    plot3(rv_orb_required_final(1,:), rv_orb_required_final(2,:), rv_orb_required_final(3,:), 'or');
    hold on;
    plot3(rv_orb_final(1,:), rv_orb_final(2,:), rv_orb_final(3,:), 'ok');
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    axis equal;
    legend('Required configuration', 'Current relative position');

end