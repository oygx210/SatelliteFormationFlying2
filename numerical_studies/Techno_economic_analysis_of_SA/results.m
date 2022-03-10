clear all;
consts = startup_formation_control();

km22 = load('22sq_km_1_day_China.mat');
km100 = load('100sq_km_1_day_China.mat');
km200 = load('200sq_km_1_day_China.mat');
km500 = load('500sq_km_1_day_China.mat');

figure;
schedule22 = km22.schedule;  
counter = 0;
for i = 1:size(schedule22,1)
    counter = counter + schedule22(i,5);
    cumulative_cost_22(i) = counter;
end

schedule100 = km100.schedule;
counter = 0;
for i = 1:size(schedule100,1)
    counter = counter + schedule100(i,5);
    cumulative_cost_100(i) = counter;
end

schedule200 = km200.schedule;
counter = 0;
for i = 1:size(schedule200,1)
    counter = counter + schedule200(i,5);
    cumulative_cost_200(i) = counter;
end

schedule500 = km500.schedule;
counter = 0;
for i = 1:size(schedule500,1)
    counter = counter + schedule500(i,5);
    cumulative_cost_500(i) = counter;
end

system_costs = (600e3 + 1.2e6)*50;
mission_profit = (cumulative_cost_200(end)*60 - system_costs)/1e6;

figure;
plot(schedule22(:,3)/60,cumulative_cost_22/1e6);
hold on;
plot(schedule100(:,3)/60,cumulative_cost_100/1e6);
hold on;
plot(schedule200(:,3)/60,cumulative_cost_200/1e6);
hold on;
plot(schedule500(:,3)/60,cumulative_cost_500/1e6);
xlabel('time,min');
ylabel('Cumulative cost, kUSD');
grid on;
legend('A_{fp} = 22 km^2', 'A_{fp} = 100 km^2','A_{fp} = 200 km^2', 'A_{fp} = 500 km^2');