function plot_ISD(t_vec, sat1_rv_ECI, sat2_rv_ECI, consts, formation, time_scale)

ISD = sqrt(dot(sat1_rv_ECI(1:3,:)-sat2_rv_ECI(1:3,:),sat1_rv_ECI(1:3,:)-sat2_rv_ECI(1:3,:)));

if time_scale == 's'
    time_res = 1;
elseif time_scale == 'm'
    time_res = 60;
elseif time_scale == 'h'
    time_res = 60*60;    
elseif time_scale == 'd'
    time_res = 60*60*24;
end

figure;
title('ISD(t)');
plot(t_vec/time_res,ISD/1000);
title('Intersatellite distance(t)');
xlabel(['time, ', time_scale]);
ylabel('ISD, km ');
grid on;
hold on;
plot(t_vec/time_res,formation.ISD/1000*ones(length(t_vec),1), '-.g');
hold on;
plot(t_vec/time_res,(formation.ISD + formation.ISD_acceptable_error)/1000*ones(length(t_vec),1), '-.r');
hold on;
plot(t_vec/time_res,(formation.ISD - formation.ISD_acceptable_error)/1000*ones(length(t_vec),1), '-.r');
legend('ISD(t)','Target ISD','Max acceptable ISD', 'Min acceptable ISD');

end