% Scipt is used to process the formation dynamics and control simulation results

%% Fuel consumption
sats(:,1) = 1:formation.N_active_sats;
figure('Name', 'Fuel consumption', 'NumberTitle', 'off');
b = barh(sats,fuel_consumption*1e3, 'stacked'); 
xlabel('Cumulative fuel consumption, g');
ylabel('Satellite');
yticks([1:formation.N_active_sats]); 
grid on; 
legend('Reconfiguration', 'Maintenance'); 
for i = 1:size(t_events,1) 
    if formation_state(i) == 1
        b(i).FaceColor = [0.8500 0.3250 0.0980]; 
    elseif formation_state(i) == 2 
        b(i).FaceColor = [0.9290 0.6940 0.1250]; 
    elseif formation_state(i) == 3 
        b(i).FaceColor = [0 0 0]; 
    end
end
for i = 1:4
    b(i).LineWidth = 0.8;
    b(i).LineStyle = '-';
end
for i = 5:8
    b(i).LineWidth = 0.8;
    b(i).LineStyle = '--';
end

[fartherst_relative_trajectory_1,index_traj_1] = max(HCW_constants_assigned(1,2:end,1));
[fartherst_relative_trajectory_2,index_traj_2] = max(HCW_constants_assigned(1,2:end,2));

% [expensivest_traj_1,index_maint_1] = max(maintenance_cost_1.fuel_consumption(:,4));
% [expensivest_traj_2,index_maint_2] = max(maintenance_cost_2.fuel_consumption(:,8));
% mean_F_m1 = mean(maintenance_cost_1.fuel_consumption(:,4));
% mean_F_m2 = mean(maintenance_cost_2.fuel_consumption(:,8));
% sigma_F_m1 = std(maintenance_cost_1.fuel_consumption(:,4));
% sigma_F_m2 = std(maintenance_cost_2.fuel_consumption(:,8));
 
% figure;
% plot(HCW_constants_assigned(1,2:end,1), 'ok');
% hold on;
% plot(maintenance_cost_1.fuel_consumption(:,4)*9000, '*k');
% 
% figure;
% plot(HCW_constants_assigned(1,2:end,2), 'ok');
% hold on;
% plot(maintenance_cost_2.fuel_consumption(:,8)*9000, '*k');
% 
%% Estimated dV consumption vs. simulated
% Simulated dV consumption
figure('Name', 'dV comparison', 'NumberTitle', 'off');
subplot(1,2,1);
b = barh(sats,[maneuvers(:,5),sum(maneuvers(:,6:8),2)], 'stacked'); 
xlabel('Cumulative dV, m/s');
ylabel('Satellite');
yticks([1:formation.N_active_sats]); 
grid on; 
legend('Reconfiguration', 'Maintenance'); 
b(1).FaceColor = [0.8500 0.3250 0.0980]; 
b(2).FaceColor = [0.9290 0.6940 0.1250]; 
title('Simulation results');
for i = 1:formation.N_active_sats
    predicted_maneuvers_reconf(i,1) = Cost_matrix_reconfiguration(match_matrix(i,1,2), match_matrix(i,2,2));
    predicted_maneuvers_maint(i,1) = cost_matrix_dV.maintenance_cost(match_matrix(i,2,2));
end

subplot(1,2,2);
c = barh(sats,[predicted_maneuvers_reconf,predicted_maneuvers_maint], 'stacked'); 
xlabel('Cumulative dV, m/s');
ylabel('Satellite');
yticks([1:26]); 
grid on; 
legend('Reconfiguration', 'Maintenance'); 
c(1).FaceColor = [0.8500 0.3250 0.0980];  
c(2).FaceColor = [0.9290 0.6940 0.1250]; 
title('Prediction according to cost matrix');

Fuel_level_before_reconf = (formation.fuel_level - sum(fuel_consumption(:,1:4),2))*1000;

% for i = 1:formation.N_active_sats
%     spacecraft_wet_mass_updated = (spacecraft.dry_mass + Fuel_level_before_reconf(i))*exp(-Cost_matrix_dV(i,:)/spacecraft.thruster_Isp/consts.g);
%     Cost_matrix_fuel(i,:) = (spacecraft.dry_mass + Fuel_level_before_reconf(i))*ones(1, formation.N_active_sats) - spacecraft_wet_mass_updated;
% end

mean_error_reconf = mean(maneuvers(:,5) - predicted_maneuvers_reconf);
mean_dV_reconf = mean(maneuvers(:,5));
error_reconf_per = mean_error_reconf/mean_dV_reconf*100;
mean_error_maintenance = mean(sum(maneuvers(:,6:8),2) - predicted_maneuvers_maint);
mean_dV_maint = mean(maneuvers(:,5));
error_maint_per = mean_error_maintenance/mean_dV_maint*100;

%% Maximin optimization
load('maximin_optimization.mat');
maximin_table = [maximin.mean_fuel_consumption, maximin.Fmin];
figure('Name', 'Pareto front');
x = maximin_table(:, 1);
y = maximin_table(:, 2);
plot(x*1e3, y*1e3, 'o-','MarkerSize', 12);
xlabel('mean fuel consumption, g');
ylabel('F_{min}, g');
grid on;

ax = [0.2 0.2];
ay = [0.5 0.5];
annotation('textarrow',ax, ay, 'String', 'Last iteration of maximin optimization; Satellite X is assigned to trajectory Y');

%% Position errors
% demo 1
figure('Name', 'Position error during demo 1', 'NumberTitle', 'off');
for i = 1:formation.N_active_sats
    plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaR(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
    hold on;
    grid on;
    xlabel('t,  UTC+2');
    ylabel('\delta\rho, m');
end
xlim([demonstration{1,1}.deployment_time demonstration{1,1}.reconfiguration_time]);
ylim([0 max(deltaR(:,1:index_events(4,2)),[],'all')]);
xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
xline(formation.orbit_epoch + seconds(t_events(1,2)));
xline(formation.orbit_epoch + seconds(150.7*60));

ax = [0.2 0.3];
ay = [0.4 0.6];

annotation('textarrow',ax, ay, 'String','end of deployment');
annotation('textarrow',ax, ay, 'String','end of impulsive control');
annotation('textarrow',ax, ay, 'String','Image 1 demonstration');

% axes('position',[.15 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaR(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(1,2)));
% 
% 
% axes('position',[.5 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaR(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(1,2)));

% demo 2
figure('Name', 'Position error during demo 2', 'NumberTitle', 'off');
for i = 1:formation.N_active_sats
    plot(formation.orbit_epoch +  seconds(t_vec(index_events(5,1):index_events(8,2))), deltaR(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
    hold on;
    grid on;
    xlabel('t,  UTC+2');
    ylabel('\delta\rho, m');
end
xlim([demonstration{2,1}.deployment_time demonstration{2,1}.reconfiguration_time]);
ylim([0 max(deltaR(:,index_events(5,1):index_events(8,2)),[],'all')]);
xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
xline(formation.orbit_epoch + seconds(t_events(5,2)));
xline(formation.orbit_epoch + seconds(t_events(5,2)) - seconds(23.4*60));

annotation('textarrow',ax, ay, 'String','end of reconfiguration');
annotation('textarrow',ax, ay, 'String','end impulsive control');
annotation('textarrow',ax, ay, 'String','Image 2 demonstration');

% axes('position',[.15 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(index_events(5,1):index_events(8,2))), deltaR(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(5,2)));
% 
% axes('position',[.5 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(index_events(5,1):index_events(8,2))), deltaR(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(5,2)));
% 
%% Velocity errors
% demo 1
figure('Name', 'Velocity error during demo 1', 'NumberTitle', 'off');
for i = 1:formation.N_active_sats
    plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaV(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
    hold on;
    grid on;
    xlabel('t,  UTC+2');
    ylabel('\delta\rho, m');
end
xlim([demonstration{1,1}.deployment_time demonstration{1,1}.reconfiguration_time]);
ylim([0 max(deltaV(:,1:index_events(4,2)),[],'all')]);
xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
xline(formation.orbit_epoch + seconds(t_events(1,2)));
xline(formation.orbit_epoch + seconds(150.7*60));

annotation('textarrow',ax, ay, 'String','end of deployment');
annotation('textarrow',ax, ay, 'String','Image 1 demonstration');
annotation('textarrow',ax, ay, 'String','end of impulsive control');

% axes('position',[.15 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaV(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(1,2)));
% 
% 
% axes('position',[.5 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(1:index_events(4,2))), deltaV(i+1,1:index_events(4,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{1,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(1,2)));
% 
% demo 2
figure('Name', 'Velocity error during demo 2', 'NumberTitle', 'off');
for i = 1:formation.N_active_sats
    plot(formation.orbit_epoch +  seconds(t_vec(index_events(5,1):index_events(8,2))), deltaV(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
    hold on;
    grid on;
    xlabel('t,  UTC+2');
    ylabel('\delta\rho, m');
end
xlim([demonstration{2,1}.deployment_time demonstration{2,1}.reconfiguration_time]);
ylim([0 max(deltaV(:,index_events(5,1):index_events(8,2)),[],'all')]);
xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
xline(formation.orbit_epoch + seconds(t_events(5,2)));
xline(formation.orbit_epoch + seconds(t_events(5,2)) - seconds(23.4*60));

annotation('textarrow',ax, ay, 'String','end of reconfiguration');
annotation('textarrow',ax, ay, 'String','Image 2 demonstration');
annotation('textarrow',ax, ay, 'String','end impulsive control');

% axes('position',[.15 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(index_events(5,1):index_events(8,2))), deltaV(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(5,2)));
% 
% axes('position',[.5 .50 .25 .25])
% box on % put box around new pair of axes
% for i = 1:formation.N_active_sats
%     plot(formation.orbit_epoch + seconds(t_vec(index_events(5,1):index_events(8,2))), deltaV(i+1,index_events(5,1):index_events(8,2)), 'k', 'LineWidth', 1);
%     hold on;
%     grid on;
%     axis tight;
% end
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(1)));
% xline(formation.orbit_epoch + seconds(demonstration{2,1}.demo_time(2)));
% xline(formation.orbit_epoch + seconds(t_events(5,2)));


%% Deployed orbital configurations and relative trajectories

sat_to_demonstrate = [55,55];
[~, image_index1] = min(abs(t_vec - (t_events(3,1)+2*pi*sqrt(formation.coe(1)^3/consts.muEarth))));

figure('Name','Orbital configuration after deployment','NumberTitle','off');
for i = 1:formation.N_active_sats
    r_relative = squeeze(rv_orb(1:3,i+1,index_events(3,1):image_index1));
    if i+1 == sat_to_demonstrate(1)
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 14);
        hold on;
    elseif i+1 == sat_to_demonstrate(2)
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'oy', 'MarkerFaceColor', 'y', 'MarkerSize', 14);
        hold on;        
    else    
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
        hold on;
    end
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 

    hold on;
end
view(180,-90);
axis equal;
grid on;
axis square;

xlabel('x, m');
ylabel('y, m');
zlabel('z, m');
legend('positions of formation satellites', 'relative trajectories');

sat_to_demonstrate = [55,55];
[~, image_index2] = min(abs(t_vec - (t_events(7,1)+2*pi*sqrt(formation.coe(1)^3/consts.muEarth))));

figure('Name','Orbital configuration after reconfiguration','NumberTitle','off');
for i = 1:formation.N_active_sats
    r_relative = squeeze(rv_orb(1:3,i+1,index_events(7,1):image_index2));
    if i+1 == sat_to_demonstrate(1)
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 14);
        hold on;
    elseif i+1 == sat_to_demonstrate(2)
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'oy', 'MarkerFaceColor', 'y', 'MarkerSize', 14);
        hold on;
    else        
        plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
        hold on;
    end
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 

    hold on;
end
view(180,-90);
axis equal;
grid on;
axis square;
xlabel('x, m');
ylabel('y, m');
zlabel('z, m');
legend('positions of formation satellites', 'relative trajectories');

%% Formation view during demonstrations

rv_orb_ima1 = squeeze(rv_orb(1:6,2:end,demo1_index));
rv_orb_ima2 = squeeze(rv_orb(1:6,2:end,demo2_index));

figure('Name', 'Formation view during demonstrations');
subplot(1,2,1);
plot3(rv_orb_ima1(1,:), rv_orb_ima1(2,:), rv_orb_ima1(3,:), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 14);
view(180,-90);
axis equal;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

subplot(1,2,2);
plot3(rv_orb_ima2(1,:), rv_orb_ima2(2,:), rv_orb_ima2(3,:), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
view(180,-90);
axis equal;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

%% Animation
% preparation 
% 1. interpolating data according to mission scenario

rv_orb_interpolated = [];
rv_orb_required_interpolated = [];
t_vec_interpolated = [];

for i = 1:length(formation_state)
    [~,event_index(i,1)] = min(abs(t_vec - t_events(i,1)));
    [~,event_index(i,2)] = min(abs(t_vec - t_events(i,2)));

    if formation_state(i) == 1
        interp_dt(i) = 40;
    elseif formation_state(i) == 2
        interp_dt(i) = 80;
    elseif formation_state(i) == 3
        interp_dt(i) = 4;
    end
    
    t_interp = [t_events(i,1):interp_dt(i):t_events(i,2)];
    rv_orb_interpolated_local = spline(t_vec(event_index(i,1):event_index(i,2)), rv_orb(:,:,event_index(i,1):event_index(i,2)), t_interp);
    rv_orb_required_interpolated_local = spline(t_vec(event_index(i,1):event_index(i,2)), rv_orb_required(:,:,event_index(i,1):event_index(i,2)), t_interp);
    if ~isempty(rv_orb_interpolated)
        l = size(rv_orb_interpolated,3);
        m = length(t_interp);
        rv_orb_interpolated(:,:,l:(l+m-1)) = rv_orb_interpolated_local;
        rv_orb_required_interpolated(:,:,l:(l+m-1)) = rv_orb_required_interpolated_local;
        t_vec_interpolated = [t_vec_interpolated, t_interp(2:end)];
    else
        rv_orb_interpolated = rv_orb_interpolated_local;
        rv_orb_required_interpolated = rv_orb_required_interpolated_local;
        t_vec_interpolated = t_interp;
    end
end

Formation_state_matrix = [t_events, formation_state, interp_dt'];

% Making animation
formation_flying_animation(t_vec_interpolated/60, rv_orb_interpolated, rv_orb_required_interpolated, Formation_state_matrix, formation, 1, 'Paris_mission_new');
