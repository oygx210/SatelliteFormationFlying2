% In Mission_design_step1.m we search for the exact time for both demonstrations at
% different Sun elevation angles

clear all;
consts = startup_formation_control();

% Data about POI and Sun cartesian position generated using STK
POI_STK = readtable('Paris_Place_Cartesian_Position.txt');
POI.t = datetime(join(string(table2cell(POI_STK(:,1:4)))), 'Format', 'dd MMM yyyy HH:mm:ss') + hours*2;
POI.r = table2array(POI_STK(:,5:7))';

Sun_STK = readtable('Sun_J2000eq_position_Paris_demo.txt');
Sun.t = datetime(join(string(table2cell(Sun_STK(:,1:4)))), 'Format', 'dd MMM yyyy HH:mm:ss') + hours*2;
Sun.r = table2array(Sun_STK(:,5:7))';

POI2Sun = Sun.r - POI.r;
Sun_elevation = rad2deg(pi/2 - acos(dot(POI2Sun,POI.r)./vecnorm(POI2Sun)./vecnorm(POI.r)));
sun_elevation_at_demo = [-2:-1:-10];

for i = 1:length(sun_elevation_at_demo)
    for j = 1:(size(Sun.t,1)-1)
        if Sun_elevation(j+1) > sun_elevation_at_demo(i)  && Sun_elevation(j) < sun_elevation_at_demo(i)
            morning_demo(i) = POI.t(j);
            morning_demo_step(i) = j;
        elseif Sun_elevation(j+1) < sun_elevation_at_demo(i)  && Sun_elevation(j) > sun_elevation_at_demo(i) 
            evening_demo(i) = POI.t(j);
            evening_demo_step(i) = j;
        end
    end
    
    time_between_demonstrations(i) = seconds(evening_demo(i) - morning_demo(i));

end

targ_orb = 5;
figure(1);
plot(POI.t, Sun_elevation);
hold on;
plot(POI.t(morning_demo_step(targ_orb)), Sun_elevation(morning_demo_step(targ_orb)), 'sk', POI.t(evening_demo_step(targ_orb)), Sun_elevation(evening_demo_step(targ_orb)), 'ok');
xlabel('Local time at POI (UTC+2)');
ylabel('Sun elevation, degrees');
legend('Sun elevation \theta', 'Sun elevation at morning demos', 'Sun elevation at evening demos');
grid on;

% axes('position',[.15 .50 .25 .25])
% box on % put box around new pair of axes
% plot(POI.t(min(morning_demo_step) - 150: max(morning_demo_step) + 150),Sun_elevation(min(morning_demo_step) - 150: max(morning_demo_step) + 150)) % plot on new axes
% hold on;
% plot(POI.t(morning_demo_step), Sun_elevation(morning_demo_step), 'sk');
% axis tight
% 
% axes('position',[.65 .50 .25 .25])
% box on % put box around new pair of axes
% plot(POI.t(min(evening_demo_step) - 150: max(evening_demo_step) + 150),Sun_elevation(min(evening_demo_step) - 150: max(evening_demo_step) + 150)) % plot on new axes
% hold on;
% plot(POI.t(evening_demo_step), Sun_elevation(evening_demo_step), 'ok');
% axis tight
