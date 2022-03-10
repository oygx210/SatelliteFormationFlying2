figure('Name','Train formation 1km','NumberTitle','off');

for i = 1:3
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satelites');
xlim([-2000 2000]);
ylim([-100 100]);
zlim ([-100 100]);

axis square;
grid on;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

figure('Name','Train formation 5 km','NumberTitle','off');

for i = 4:6
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satelites');
xlim([-6000 6000]);
ylim([-500 500]);
zlim ([-500 500]);

axis square;
grid on;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

figure('Name','Train formation 10 km','NumberTitle','off');

for i = 7:9
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satelites');
xlim([-11000 11000]);
ylim([-1000 1000]);
zlim ([-1000 1000]);

axis square;
grid on;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

orbit_period = round(2*pi*sqrt(oe(1)^3/consts.muEarth));

figure('Name','GCO formation','NumberTitle','off');
subplot(3,1,1);
for i = 1:3
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
end
axis square;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');


subplot(3,1,2);
for i = 4:6
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
end
axis square;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');


subplot(3,1,3);
for i = 7:9
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
    xlim([-12000; 12000]);
    ylim([-12000; 12000]);
    zlim([-12000; 12000]);

end
axis square;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

figure('Name','Tetrahedron formation 1 km','NumberTitle','off');
r_relative = squeeze(rv_orb.tetrahedron(1:3,1:4,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

figure('Name','Tetrahedron formation 5 km','NumberTitle','off');
r_relative = squeeze(rv_orb.tetrahedron(1:3,5:8,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');

figure('Name','Tetrahedron formation 10 km','NumberTitle','off');
r_relative = squeeze(rv_orb.tetrahedron(1:3,9:12,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel('x, meters');
ylabel('y, meters');
zlabel('z, meters');
