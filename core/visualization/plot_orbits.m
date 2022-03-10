function plot_orbits(rv, consts)

T_orb_max = 0;
mu = consts.muEarth;
    for i = 1:(size(rv,1)/6)
        r = vecnorm(rv(6*i-5:6*i-3));
        v = vecnorm(rv(6*i-2:6*i));        
        ksi = v^2/2 - mu/r;
        a = - mu / 2 / ksi;
        T_orb = 2*pi*sqrt(a^3/mu);
        if T_orb > T_orb_max
            T_orb_max = T_orb;
        end
    end
    
    options_precision = odeset('RelTol',1e-6,'AbsTol',1e-6);
    spacecraft = [];
    [~, rv_ECI] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [0 T_orb_max], rv, options_precision);
    rv_ECI = rv_ECI'; 
    
    % Orbit visualization
%     sun_direction = [zeros(3,1), R_sun_unit_initial];
    ox = [[0;0;0],...
          [1;0;0]];
    oy = [[0;0;0],...
          [0;1;0]];
    oz = [[0;0;0],...
          [0;0;1]];

    figure('Name','Orbits visualization', 'NumberTitle', 'Off');
%     earth_sphere('m');
%     hold on;
    for i = 1: (size(rv,1)/6)
        a = plot3(rv_ECI(6*i-5,:), rv_ECI(6*i-4,:), rv_ECI(6*i-3,:), 'k-');
        hold on
        b = plot3(rv_ECI(6*i-5,1), rv_ECI(6*i-4,1), rv_ECI(6*i-3,1), 'o');
        hold on;
    end
        hold on;
        c1 = plot3(ox(1,:), ox(2,:), ox(3,:),'r');
        c2 = plot3(oy(1,:), oy(2,:), oy(3,:),'g');
        c3 = plot3(oz(1,:), oz(2,:), oz(3,:),'b');
        axis equal;
%         legend([a,b,c1, c2, c3], 'Orbits', 'Initial conditions', 'x-axis', 'y-axis', 'z-axis');
        grid on;
