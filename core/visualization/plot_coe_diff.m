function plot_coe_diff(t_vec, sat1_rv_ECI, sat2_rv_ECI, consts, time_scale)

for i = 1:size(sat1_rv_ECI,2)
    sat1.coe = rv2oe(sat1_rv_ECI(:,i), consts);
    sat2.coe = rv2oe(sat2_rv_ECI(:,i), consts);
    delta_coe(:,i) = sat1.coe - sat2.coe;
end
if time_scale == 's'
    time_res = 1;
elseif time_scale == 'm'
    time_res = 60;
elseif time_scale == 'h'
    time_res = 60*60;    
elseif time_scale == 'd'
    time_res = 60*60*24;
end

figure(1);
subplot(2,3,1);
plot(t_vec/time_res, delta_coe(1,:));
ylabel('\deltasma');
xlabel(['time, ', time_scale]);
subplot(2,3,2);
plot(t_vec/time_res, delta_coe(2,:));
ylabel('\deltaecc');
xlabel(['time, ', time_scale]);
title('Orbital elements differents');
subplot(2,3,3);
plot(t_vec/time_res, delta_coe(3,:));
ylabel('\deltainc');
xlabel(['time, ', time_scale]);
subplot(2,3,4);
plot(t_vec/time_res, delta_coe(4,:));
ylabel('\deltaRAAN');
xlabel(['time, ', time_scale]);
subplot(2,3,5);
plot(t_vec/time_res, delta_coe(5,:));
ylabel('\deltaAOP');
xlabel(['time, ', time_scale]);
subplot(2,3,6);
plot(t_vec/time_res, sin(delta_coe(7,:)));
ylabel('\delta Mean anomaly');
xlabel(['time, ', time_scale]);

end