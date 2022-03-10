function plot_mean_oe_diff(t_vec, sat1_rv_ECI, sat2_rv_ECI, consts, time_scale)

sat1.coe = vecRV2OE(sat1_rv_ECI, consts);
sat2.coe = vecRV2OE(sat2_rv_ECI, consts);
for i = 1:length(t_vec)
    sat1.mean_oe(:,i) = osc2mean(sat1.coe(:,i), consts);
    sat2.mean_oe(:,i) = osc2mean(sat2.coe(:,i), consts);
end   
delta_mena_oe = sat1.mean_oe - sat2.mean_oe;

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
subplot(2,3,1);
plot(t_vec/time_res, delta_mena_oe(1,:));
ylabel('\deltasma');
xlabel(['time, ', time_scale]);
subplot(2,3,2);
plot(t_vec/time_res, delta_mena_oe(2,:));
ylabel('\deltaecc');
xlabel(['time, ', time_scale]);
title('Orbital elements differents');
subplot(2,3,3);
plot(t_vec/time_res, delta_mena_oe(3,:));
ylabel('\deltainc');
xlabel(['time, ', time_scale]);
subplot(2,3,4);
plot(t_vec/time_res, delta_mena_oe(4,:));
ylabel('\deltaRAAN');
xlabel(['time, ', time_scale]);
subplot(2,3,5);
plot(t_vec/time_res, delta_mena_oe(5,:));
ylabel('\deltaAOP');
xlabel(['time, ', time_scale]);
subplot(2,3,6);
plot(t_vec/time_res, sin(delta_mena_oe(7,:)));
ylabel('\delta Mean anomaly');
xlabel(['time, ', time_scale]);

end