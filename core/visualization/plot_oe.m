function plot_oe(t_vec, oe, time_scale)

if time_scale == 's'
    time_res = 1;
elseif time_scale == 'm'
    time_res = 60;
elseif time_scale == 'h'
    time_res = 60*60;    
elseif time_scale == 'd'
    time_res = 60*60*24;
end

subplot(2,4,1);
plot(t_vec/time_res, oe(1,:));
ylabel('sma, m');
xlabel(['time, ', time_scale]);

subplot(2,4,2);
plot(t_vec/time_res, oe(2,:));
ylabel('ecc, -');
xlabel(['time, ', time_scale]);

title('Orbital elements');

subplot(2,4,3);
plot(t_vec/time_res, oe(3,:));
ylabel('inc, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,4);
plot(t_vec/time_res, oe(4,:));
ylabel('RAAN, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,5);
plot(t_vec/time_res, oe(5,:));
ylabel('AOP, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,6);
plot(t_vec/time_res, (oe(6,:)));
ylabel('\nu, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,7);
plot(t_vec/time_res, (oe(7,:)));
ylabel('M, rad');
xlabel(['time, ', time_scale]);

subplot(2,4,8);
plot(t_vec/time_res, (oe(8,:)));
ylabel('u, rad');
xlabel(['time, ', time_scale]);

end