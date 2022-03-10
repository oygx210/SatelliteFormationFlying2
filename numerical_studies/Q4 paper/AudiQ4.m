clear all;

consts = startup_formation_control();

% Scipt with simulations for AudiQ4 paper

%% Defining target orbit
oe = zeros(6,1);
oe(1) = 700e3 + consts.rEarth;
oe(2) = 0;
oe(3) = get_SSO_inclination(oe(1), oe(2), consts);

Rsun = [26127801; -132825709.3; -57579560.5]; % on Jan 1, 2022
angular_momentum_x = Rsun(1)/vecnorm(Rsun);
angular_momentum_y = Rsun(2)/vecnorm(Rsun);
node_vector = [-angular_momentum_y; angular_momentum_x; 0];
node_vector = node_vector./vecnorm(node_vector);

if node_vector(2) >= 0
    oe(4) = acos(node_vector(1));
else
    oe(4) = 2*pi - acos(node_vector(1));
end

oe(5) = 0;
oe(6) = 0;

orbit_epoch = datetime(2022, 1, 1, 0, 0, 0);

oe_deg = [oe(1); oe(2); 180/pi*oe(3); 180/pi*oe(4); 180/pi*oe(5); 180/pi*oe(6)];

target_orbit.oe = oe;
target_orbit.rv = oe2rv(oe, consts);
target_orbit.epoch = datetime(2022, 1, 1, 0, 0, 0);

mean_motion = sqrt(consts.muEarth/target_orbit.oe(1)^3);

%% Defining orbital configuraitons
r1 = 1e3; % ISD, meters
r2 = 5e3; % ISD, meters
r3 = 10e3; % ISD, meters

orbital_configuration.train = [[0;0;0;0;0],...
                                 [0;0;r1;0;0],...
                                 [0;0;-r1;0;0],...
                                 [0;0;0;0;0],...
                                 [0;0;r2;0;0],...
                                 [0;0;-r2;0;0],...                                
                                 [0;0;0;0;0],...
                                 [0;0;r3;0;0],...
                                 [0;0;-r3;0;0]];
                                                              
orbital_configuration.GCO = [[0;0;0;0;0],...
                             [r1; sqrt(3)/2*r1; 0;0;0],...
                             [r1; sqrt(3)/2*r1; 0;pi;pi],...
                             [0;0;0;0;0],...
                             [r2; sqrt(3)/2*r2; 0;0;0],...
                             [r2; sqrt(3)/2*r2; 0;pi;pi],...                             
                             [0;0;0;0;0],...
                             [r3; sqrt(3)/2*r3; 0;0;0],...
                             [r3; sqrt(3)/2*r3; 0;pi;pi]];

orbital_configuration.tetrahedron = [[0;0;0;0;0],...
                                      [2*r1/5; 0; 2*r1*sqrt(5/3); 0; 0],...
                                      [2*r1; r1*sqrt(5); r1*sqrt(5/3); -atan(1/sqrt(2)); atan(sqrt(2)) - pi],...
                                      [2*r1; r1*sqrt(5); r1*sqrt(5/3); atan(1/sqrt(2)); -atan(sqrt(2))],...
                                      [0;0;0;0;0],...
                                      [2*r2/5; 0; 2*r2*sqrt(5/3); 0; 0],...
                                      [2*r2; r2*sqrt(5); r2*sqrt(5/3); -atan(1/sqrt(2)); atan(sqrt(2)) - pi],...
                                      [2*r2; r2*sqrt(5); r2*sqrt(5/3); atan(1/sqrt(2)); -atan(sqrt(2))],...
                                      [0;0;0;0;0],...
                                      [2*r3/5; 0; 2*r3*sqrt(5/3); 0; 0],...
                                      [2*r3; r3*sqrt(5); r3*sqrt(5/3); -atan(1/sqrt(2)); atan(sqrt(2)) - pi],...
                                      [2*r3; r3*sqrt(5); r3*sqrt(5/3); atan(1/sqrt(2)); -atan(sqrt(2))]];                         
                         
rv_orb_initial.train = initital_conditions(orbital_configuration.train, mean_motion, consts);
rv_orb_initial.GCO = initital_conditions(orbital_configuration.GCO, mean_motion, consts);
rv_orb_initial.tetrahedron = initital_conditions(orbital_configuration.tetrahedron, mean_motion, consts);

for i = 1:size(rv_orb_initial.train,2)
    rv_ECI_initial.train(i*6-5:i*6,1) = orb2ECI(target_orbit.rv,rv_orb_initial.train(:,i), consts);
    rv_ECI_initial.GCO(i*6-5:i*6,1) = orb2ECI(target_orbit.rv,rv_orb_initial.GCO(:,i), consts);
end

for i = 1:size(rv_orb_initial.tetrahedron,2)
    rv_ECI_initial.tetrahedron(i*6-5:i*6,1) = orb2ECI(target_orbit.rv,rv_orb_initial.tetrahedron(:,i), consts);
end

%% dynamics simulation
global environment
environment = 'point mass';
spacecraft = [];

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec.train, rv_ECI.train] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[1:consts.day2sec], rv_ECI_initial.train, options_precision);        
rv_ECI.train = rv_ECI.train';

[t_vec.GCO, rv_ECI.GCO] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[1:consts.day2sec], rv_ECI_initial.GCO, options_precision);        
rv_ECI.GCO = rv_ECI.GCO';

[t_vec.tetrahedron, rv_ECI.tetrahedron] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),[1:consts.day2sec], rv_ECI_initial.tetrahedron, options_precision);        
rv_ECI.tetrahedron = rv_ECI.tetrahedron';

%% results

for i = 1:size(orbital_configuration.train,2)
    for j = 1:length(t_vec.train)
        rv_orb.train(:,i,j) = ECI2orb(rv_ECI.train(1:6,j), rv_ECI.train(6*i-5:6*i,j), consts);
    end
end

for i = 1:size(orbital_configuration.GCO,2)
    for j = 1:length(t_vec.train)
        rv_orb.GCO(:,i,j) = ECI2orb(rv_ECI.GCO(1:6,j), rv_ECI.GCO(6*i-5:6*i,j), consts);
    end
end

for i = 1:size(orbital_configuration.tetrahedron,2)
    for j = 1:length(t_vec.tetrahedron)
        rv_orb.tetrahedron(:,i,j) = ECI2orb(rv_ECI.tetrahedron(1:6,j), rv_ECI.tetrahedron(6*i-5:6*i,j), consts);
    end
end

ox = [zeros(3,1),[500;0;0]];
oy = [zeros(3,1),[0;500;0]];
oz = [zeros(3,1),[0;0;500]];

xaxis = ['\it{x}', '\rm{, m}'];
yaxis = ['\it{y}', '\rm{, m}'];
zaxis = ['\it{z}', '\rm{, m}'];

figure('Name','Train formation','NumberTitle','off');
subplot(1,3,1);

for i = 1:3
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satellites');
xlim([-2000 2000]);
ylim([-100 100]);
zlim ([-100 100]);

axis square;
grid on;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

subplot(1,3,2);
for i = 4:6
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satellites');
xlim([-6000 6000]);
ylim([-500 500]);
zlim ([-500 500]);

axis square;
grid on;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

subplot(1,3,3);
for i = 7:9
    r_relative = squeeze(rv_orb.train(1:3,i,:));
    plot3(r_relative(1,1), r_relative(2,1), r_relative(3,1),  'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
%     hold on;
%     plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
end
legend('satellites');
xlim([-11000 11000]);
ylim([-1000 1000]);
zlim ([-1000 1000]);

axis square;
grid on;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

orbit_period = round(2*pi*sqrt(oe(1)^3/consts.muEarth));

figure('Name','GCO formation','NumberTitle','off');

subplot(3,1,1);
for i = 1:3
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
end
axis square;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);


subplot(3,1,2);
for i = 4:6
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
end
axis square;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);


subplot(3,1,3);
for i = 7:9
    r_relative = squeeze(rv_orb.GCO(1:3,i,1:orbit_period));
    plot3(r_relative(1,end), r_relative(2,end), r_relative(3,end),  'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 14);
    hold on;
    plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'k', 'LineWidth', 0.3); 
    hold on;
    grid on;
    xlim([-12000; 12000]);
    ylim([-12000; 12000]);
    zlim([-12000; 12000]);

end
axis square;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

figure('Name','Tetrahedron formation','NumberTitle','off');
subplot(1,3,1);
r_relative = squeeze(rv_orb.tetrahedron(1:3,1:4,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

subplot(1,3,2);
r_relative = squeeze(rv_orb.tetrahedron(1:3,5:8,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

subplot(1,3,3);
r_relative = squeeze(rv_orb.tetrahedron(1:3,9:12,1));
data = [r_relative, r_relative(:,1)];
plot3(r_relative(1,:), r_relative(2,:), r_relative(3,:), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 0.3); 
hold on;
grid on;
axis square;
axis equal;
xlabel(xaxis);
ylabel(yaxis);
zlabel(zaxis);

% animation(t_vec.GCO, rv_orb.tetrahedron(:,1:4,1:orbit_period), 'GCO');

%% functions

function animation(t_vec, rv_orb, VideoHeader)

    fig = figure;
    plot3(rv_orb(1,1,1), rv_orb(2,1,1), rv_orb(3,1,1), 'ok', 'MarkerSize', 1);
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    xlim([-5000 5000]);
    ylim([-5000 5000]);
    zlim([-5000 5000]);
    pbaspect([1 1 1]);
    title('Orbital reference frame');
    view(135,45);
    hold on;

    step = 10;
    for i = 1:size(rv_orb,3)/step
        a = plot3(rv_orb(1,:,i*step), rv_orb(2,:,i*step), rv_orb(3,:,i*step), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
        hold on;
        drawnow;
        grid on;
%         T = t_vec(i*step);
%         info = {['Simulation time: ', T]};
%         an = annotation('textbox',[0.2 0.025 0.6 0.15], 'String', info, 'FitBoxToText', 'off', 'HorizontalAlignment', 'center'); 

    %   Take a Snapshot
        movieVector(i) = getframe(fig);   %manually specify getframe region    

        delete(a);
    end
      
    myWriter = VideoWriter(VideoHeader,'MPEG-4');   %create an .mp4 file
    myWriter.FrameRate = 24; %

    %   Open the VideoWriter object, write the movie, and close the file
    open(myWriter);
    writeVideo(myWriter, movieVector);
    close(myWriter); 
end   

function rv_orb = initital_conditions(HCW_constants, mean_motion, consts)

% HCW_constants = [C1, C2, C3, alpha, beta];

% Notation for orbital reference frame: x - along track; y - out-of-plane; z - radial

% Let's assume that arg_of_lat is equal to zero at the initial epoch

% x = c1 * cos(alpha);
% y = c2 * sin(beta);
% z = c1/2 * sin(alpha);
% dx/dt = - n * c1 * sin(alpha);
% dy/dt =  n * c2 * cos(beta);
% dz/dt =  n * c1/2 * cos(alpha);

C1 = HCW_constants(1,:);
C2 = HCW_constants(2,:);
C3 = HCW_constants(3,:);
alpha = HCW_constants(4,:);
beta = HCW_constants(5,:);

r_x = C1.*cos(alpha) + C3;
r_y = C2.*sin(beta);
r_z = (C1/2).*sin(alpha);
v_x = -mean_motion*C1.*sin(alpha);
v_y = mean_motion*C2.*cos(beta);
v_z = mean_motion*C1/2.*cos(alpha);

rv_orb = [r_x ; r_y ; r_z ; v_x ; v_y ; v_z];

end
