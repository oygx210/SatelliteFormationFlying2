function plot_eqoe_diff(t_vec, sat1_rv_ECI, sat2_rv_ECI, consts, time_scale)

for i = 1:size(sat1_rv_ECI,2)
    sat1.coe = rv2oe(sat1_rv_ECI(:,i), consts);
    sat1.eqoe = coe2oe_equinoctial(sat1.coe);
    
    sat2.coe = rv2oe(sat2_rv_ECI(:,i), consts);
    sat2.eqoe = coe2oe_equinoctial(sat2.coe);
    
    delta_eqoe(:,i) = sat1.eqoe - sat2.eqoe;
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

figure;
subplot(2,3,1);
plot(t_vec/time_res, delta_eqoe(1,:));
ylabel('\deltasma');
xlabel(['time, ', time_scale]);
subplot(2,3,2);
plot(t_vec/time_res, delta_eqoe(2,:));
ylabel('\deltaq1 = e*cos(\omega)');
xlabel(['time, ', time_scale]);
title('Orbital elements differents');
subplot(2,3,3);
plot(t_vec/time_res, delta_eqoe(3,:));
ylabel('\deltaq1 = e*sin(\omega)');
xlabel(['time, ', time_scale]);
subplot(2,3,4);
plot(t_vec/time_res, delta_eqoe(4,:));
ylabel('\deltai');
xlabel(['time, ', time_scale]);
subplot(2,3,5);
plot(t_vec/time_res, delta_eqoe(5,:));
ylabel('\delta\Omega');
xlabel(['time, ', time_scale]);
subplot(2,3,6);
plot(t_vec/time_res, sin(delta_eqoe(6,:)));
ylabel('\delta\lambda');
xlabel(['time, ', time_scale]);

end