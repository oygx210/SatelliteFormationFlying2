 % Run startup_formation_control
% 1. Formation Establishment and Reconfiguration Using Impulsive Control
% Repeating the experiment from the [1]

clear all;
consts = startup_formation_control_func();

oe_in(1) = 700e3 + consts.rEarth; % sma, km
oe_in(2) = 0; % ecc
oe_in(3) = get_SSO_inclination(oe_in(1), oe_in(2), consts); % inc, rad
oe_in(4) = 0; % RAAN, rad
oe_in(5) = 0; % AOP, rad
oe_in(6) = 0; % Mean anomaly, rad

mean_motion = sqrt(consts.muEarth/oe_in(1)^3);

orbit_period = 2 * pi / mean_motion;

target.rv_in = oe2rv(oe_in, consts);
chaser.rv_in = oe2rv(oe_in, consts);

%% Chaser satellite desired orbit

dist = 230;%[m]
[~,~,c1,c2,alpha] = positions(dist,"Sk",mean_motion);
figure;
plot(c1.*cos(alpha), c1.*sin(alpha), 'ok', 'MarkerSize', 12);
axis off
numsats = length(c1);

[dq1,dq2,dinc,dOmega,dlambda] = deltaOE(c1,c2,alpha,oe_in(1),oe_in(3));

%% Impulses
p=1; %basically, ratio between first and second out-of-plane burns. See Vaddi et al., page 264, section called "Aside"
gamma = sqrt(oe_in(1)/consts.muEarth);  %[s/m]

dV1 = [];
dV2 = [];

for sat=1:numsats
    %deltaVs are in m/s 
    dV_h1 = p/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(oe_in(3)).^2);
    dV_h2 = (1-p)/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(oe_in(3)).^2);
    dV_r1 = -sqrt(dq1.^2 + dq2.^2)/2/gamma;
    dV_r2 = sqrt(dq1.^2 + dq2.^2)/2/gamma;

    dV1(:, sat) = [0 ; dV_h1(sat) ; dV_r1(sat)];
    dV2(:, sat) = [0 ; dV_h2(sat) ; dV_r2(sat)];
end
%% Creating the initial condition for chaser satellite after the first impulse
chaser = [];
for sat=1:numsats
    chaser(sat).rv_orb_in = [0; 0; 0; dV1(:,sat)]; 
                     
    chaser(sat).rv_ECI_in = orb2ECI(target.rv_in, chaser(sat).rv_orb_in, consts);
end
                
%% Mission control

model.dt = 20; % s
model.simulation_time_phase1 = orbit_period/2/pi*pi; % s
model.simulation_time_phase2 = orbit_period/2/pi*pi*2; % s
model.simulation_total_time = model.simulation_time_phase1 + model.simulation_time_phase2;

tspan_phase1 = [0 : model.dt : model.simulation_time_phase1];
tspan_phase2 = [model.simulation_time_phase1 + 1 : model.dt : model.simulation_total_time];
model.dt_vec = zeros(length(tspan_phase1) + length(tspan_phase2), 1);

options = odeset('RelTol',1e-12,'AbsTol',1e-12); % random values

target.rv_ECI = zeros(length(model.dt_vec), 6);
for sat=1:numsats
    chaser(sat).rv_ECI = zeros(length(model.dt_vec), 6);
end


%% Propagation after the first impulse

target_init = target.rv_in;
target_tspan = [-200 : model.dt : 3*orbit_period];
[model.dt_vec, target.rv_ECI] = ode45(@(t, rv)central_gravity(t, rv, consts), target_tspan, target_init, options);
target.rv_ECI = target.rv_ECI';

for sat=1:numsats
    chaser(sat).rv_orb = singlesat_integrator(alpha(sat),orbit_period,dV1(:,sat),dV2(:,sat),target.rv_ECI,consts);
end

figh = figure(1);
subplot(1,2,1);
plot_Earth_meters();
hold on;
axis tight manual;
 
% ECI
plot3(target.rv_ECI(1,:), target.rv_ECI(2,:), target.rv_ECI(3,:), 'LineWidth', 1, 'Color', 'k');
xlabel('x, m', 'fontsize',12);
ylabel('y, m', 'fontsize',12);
zlabel('z, m', 'fontsize',12);
title('Earth-centered inetrial');
view(-225,45);
xlim([-7500e3 7500e3]);
ylim([-7500e3 7500e3]);
zlim([-7500e3 7500e3]);
pbaspect([1 1 1]);

numsats = length(c1);
subplot(1,2,2);
axis tight manual;
xlim([-2000 2000]);
ylim([-2000 2000]);
zlim([-2000 2000]);
view([1 -1 -1]);
pbaspect([1 1 1]);
title('Orbital reference frame');
xlabel('x [m](along-track)', 'fontsize',12);
ylabel('y [m](normal)', 'fontsize',12);
zlabel('z [m](radial)', 'fontsize',12);

for i = (1+10):size(chaser(1).rv_orb,2)-2
    
    subplot(1,2,1);
    ECI = plot3(target.rv_ECI(1,i), target.rv_ECI(2,i), target.rv_ECI(3,i),'or', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    xlabel('x, m', 'fontsize',12);
    ylabel('y, m', 'fontsize',12);
    zlabel('z, m', 'fontsize',12);
    title('Earth-centered inetrial');
    view(-225,45);
    xlim([-7500e3 7500e3]);
    ylim([-7500e3 7500e3]);
    zlim([-7500e3 7500e3]);
    pbaspect([1 1 1]);

    subplot(1,2,2);
    for j=1:numsats
        if mod(i,5)==1
            trail = plot3(chaser(j).rv_orb(1,i),chaser(j).rv_orb(2,i),chaser(j).rv_orb(3,i),'og-','MarkerSize', 1,'MarkerFaceColor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.4,0.4,0.4]); 
        end
            hold on;
            if j == 12
                trail2(j) = plot3(chaser(j).rv_orb(1,i),chaser(j).rv_orb(2,i),chaser(j).rv_orb(3,i), 'sr', 'MarkerFaceColor', 'r', 'MarkerSize',8);

            else
                trail2(j) = plot3(chaser(j).rv_orb(1,i),chaser(j).rv_orb(2,i),chaser(j).rv_orb(3,i), 'sk', 'MarkerFaceColor', 'k', 'MarkerSize',7);

            end
    end
    axis tight manual;
    xlim([-2000 2000]);
    ylim([-2000 2000]);
    zlim([-2000 2000]);
    view([1 -1 -1]);
    pbaspect([1 1 1]);
    title('Orbital reference frame');
    xlabel('x [m](along-track)', 'fontsize',12);
    ylabel('y [m](normal)', 'fontsize',12);
    zlabel('z [m](radial)', 'fontsize',12);

%     text(6000,-6000,"Time,seconds")
%     t2=text(5500,-5500,num2str(i*model.dt));
%     text(5000,-5000,"Time, orbits")
%     t4=text(4500,-4500,num2str(i*model.dt/orbit_period,3));

    %% Step 3: Take a Snapshot
    movieVector(i) = getframe(figh);   %manually specify getframe region    

    delete(ECI);
%     delete(t2);
%     delete(t4);
    
    for j = 1:numsats
        delete(trail2(j));
    end
end

%% Step 5: Save Movie

myWriter = VideoWriter('Letters_in_sky','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector(31:end));
close(myWriter);

disp('DONE!')

%%%%%%%%%
% 1. Prepare letter Sk & look at Danil's function
% 2. Plan a movie
    % Flying letters Sk in orbital frame;
    % Leader satellite flying around the Earth;
    % Short description of the process
% 3. Construct a movie
% 4. Prepare visualization function for formation control

% Work on equinoctial element - check how orbits changes along orbit
% Read the Vaddi's paper and think about control derivation
% Create a simple script for maneuvers check for a pair of satellites

% Start thinking on minimax optimization - watch videos on youtube about
% brutforce method and think about its implementation

% Control scheme 
% Define test letters AB and BA and think about reconfiguration method
% Maximin can be used at every reconfigureation actually

% Think how to maneuver only those satellites that are out of threshold
%%%%%%%%

% 
%% Dynamic visualization
% subplot(1,2,1);
% earth_sphere('m');
% hold on;
% 
% % ECI
% plot3(target.rv_ECI(:,1), target.rv_ECI(:,2), target.rv_ECI(:,3), 'LineWidth', 1, 'Color', 'k');
% xlabel('x', 'fontsize',12);
% ylabel('y', 'fontsize',12);
% zlabel('z', 'fontsize',12);
% title('Earth-centered inetrial');
% view(-225,45);

% colourmap = jet(numsats);
% set(gca,'nextplot','replacechildren','visible','off')
% for i = (1+10):(length(chaser(1).rv_orb(1,:))-2)
%     
% %     subplot(1,2,1);
% %     plot3(target.rv_ECI(i,1), target.rv_ECI(i,2), target.rv_ECI(i,3),'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
% %     
% %     subplot(1,2,2);
%     view([1 -1 -1])
%     pbaspect([1 1 1]);
%     xlabel('x [m](along-track)', 'fontsize',12);
%     ylabel('y [m](normal)', 'fontsize',12);
%     zlabel('z [m](radial)', 'fontsize',12);
%     xlim([-4000 4000]);
%     ylim([-4000 4000]);
%     zlim([-4000 4000]);
%     
%     for sat=1:numsats
%         if mod(i,5)==1
%             trail = scatter3(chaser(sat).rv_orb(1,i),chaser(sat).rv_orb(2,i),chaser(sat).rv_orb(3,i),3,'MarkerFaceColor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.4,0.4,0.4]); 
%             trail.MarkerFaceAlpha = .3;
%             trail.MarkerEdgeAlpha = .3;
%         end
%         scatter2(sat) = scatter3(chaser(sat).rv_orb(1,i),chaser(sat).rv_orb(2,i),chaser(sat).rv_orb(3,i),22, 'MarkerFaceColor', 'r','MarkerEdgeColor',colourmap(sat,:));
%         hold on;
%     end
%     pause(0.00001);
%     if i==length(chaser(1).rv_orb(1,:))-3
%         break
%     end
%   
%     text(6000,-6000,"Time,seconds")
%     t2=text(5500,-5500,num2str(i*model.dt))
%     text(5000,-5000,"Time, orbits")
%     t4=text(4500,-4500,num2str(i*model.dt/orbit_period,3))
% 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1 
%         imwrite(imind,cm,filename,'gif89a', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
%     
%     for sat=1:numsats
%         delete(scatter2(sat));
%     end
%     delete(t2)
%     delete(t4)
% end

% subplot(1,2,1);
% plot3(target.rv_ECI(1,1), target.rv_ECI(1,2), target.rv_ECI(1,3));



% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% for n = 1:0.5:5
%     % Draw plot for y = x.^n
%     x = 0:0.01:1;
%     y = x.^n;
%     plot(x,y) 
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
%   end
