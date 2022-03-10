function coverage_animation(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, cost_matrix, schedule, consts)

    dt = round((t_vec_JD(2) - t_vec_JD(1))*consts.day2sec);
    t_simulation = round((t_vec_JD(end) - t_vec_JD(1))*consts.day2sec);
    % sphere generation
    [M, N] = size(cost_matrix);
    phi_range = linspace(-pi, pi, N+1);
    phi_range = phi_range(1:end-1);
    th_range = linspace(-pi/2, pi/2, M+1);
    th_range = th_range(1:end-1);
    [PHI, TH] = meshgrid(phi_range, th_range);
    
    % trigonometric matrices
    cTH = cos(TH);
    sTH = sin(TH);
    cPHI = cos(PHI);
    sPHI = sin(PHI);
    r_orbit = vecnorm(r_ECEF(:,1));
    b = consts.rEarth/r_orbit;
    r_POI_ij = consts.rEarth*cat(3, cTH.*cPHI, cTH.*sPHI, sTH); % POIs coordinates
    
    % calculating critical angles
    gamma_max = asin(b*cos(theta_sat_min));
    beta_max = pi/2 - gamma_max - theta_sat_min;

    legend_array = [];
    plot_cost_map = 'yes';
    switch plot_cost_map
        case 'yes'
            cost_matrix_mirrored = flip(cost_matrix,1);

            logical = cost_matrix_mirrored < 100;
            cost_matrix_100 = cost_matrix_mirrored;
            cost_matrix_100(~logical) = 0;

            logical = cost_matrix_mirrored > 100 & cost_matrix_mirrored < 500;
            cost_matrix_100t500 = cost_matrix_mirrored;
            cost_matrix_100t500(~logical) = 0;

            logical = cost_matrix_mirrored > 500;
            cost_matrix_b500 = cost_matrix_mirrored;
            cost_matrix_b500(~logical) = 0;

            [cost_matrix_100_index2, cost_matrix_100_index1] = find(cost_matrix_100);
            [cost_matrix_100t500_index2, cost_matrix_100t500_index1] = find(cost_matrix_100t500);
            [cost_matrix_b500_index2, cost_matrix_b500_index1] = find(cost_matrix_b500);
    end
    figure;
    earth = imread('earthmap1k.jpg');
    I = imshow(earth);
    hold on;
    xlabel('longitude \lambda, deg');
    ylabel('latitude \phi, deg');
    axis on;
    xticks([0 125 250 375 500 625 750 875 1000]);
    xticklabels({'-180','-135','-90', '-45', '0', '45','90', '135', '180'});
    yticks([0 125 250 375 500]);
    yticklabels(flip({'-90', '-45', '0', '45', '90'}));
    gif('Coverage_animation.gif','frame',gcf);

    switch plot_cost_map
        case 'yes'
        plot(cost_matrix_100_index1, cost_matrix_100_index2,'o', 'MarkerEdgeColor', '#EDB120', 'MarkerFaceColor', '#EDB120', 'MarkerSize', 1.5);
        hold on;           
        plot(cost_matrix_100t500_index1, cost_matrix_100t500_index2,'o', 'MarkerEdgeColor', '#D95319', 'MarkerFaceColor', '#D95319', 'MarkerSize', 2);
        hold on;
        plot(cost_matrix_b500_index1, cost_matrix_b500_index2,'o', 'MarkerEdgeColor', '#A2142F', 'MarkerFaceColor', '#A2142F', 'MarkerSize', 2.5);
        hold on;
        grid on;
        legend_array = [legend_array,...
                "Demonstration cost C < 100 USD",...
                "Demonstration cost 100 < C < 500 USD",...
                "Demonstration cost C > 500 USD"];
    end
    T_sat = 2*pi*sqrt(vecnorm(r_ECEF(:,1))^3/consts.muEarth);
    i_local = 0;
    T_in_steps = round(T_sat / dt);
    while i_local < t_simulation/dt            
        [azimuth,elevation,~] = cart2sph(r_ECEF(1,i_local+1:i_local+T_in_steps),r_ECEF(2,i_local+1:i_local+T_in_steps),r_ECEF(3,i_local+1:i_local+T_in_steps));
        azimuth = azimuth + pi;
        elevation = elevation + pi/2;
        step_azimuth = 2*pi/size(cost_matrix,2);
        step_elevation = pi/size(cost_matrix,1);

        azimuth_index = round((azimuth - mod(azimuth, step_azimuth))./step_azimuth);
        elevation_index = round((elevation - mod(azimuth, step_azimuth))./step_elevation);
        elevation_index = size(cost_matrix,1) - elevation_index;
    
        GT = plot(azimuth_index, elevation_index, 'ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 1);
        hold on;

        for i = (i_local+1):(i_local+T_in_steps)
    
            r_sat = reshape(r_ECEF(:,i), [1,1,3]);
            r_sat_repmat = repmat(r_sat, size(PHI));
    
            r_sun_ECI = sun(t_vec_JD(i))*consts.AstronomicUnit;
            r_sun_ECEF = ECItoECEF(t_vec_JD(i),r_sun_ECI');
            r_sun_ECEF = reshape(r_sun_ECEF, [1,1,3]);
            r_sun_ECEF_repmat = repmat(r_sun_ECEF, size(PHI));
    
            r_POI2sun = r_sun_ECEF_repmat - r_POI_ij;
            e_POI2sun = r_POI2sun ./ vecnorm(r_POI2sun,2,3);
            e_POI_ij = r_POI_ij ./ vecnorm(r_POI_ij,2,3);
            e_sat_repmat = r_sat_repmat ./ vecnorm(r_sat_repmat, 2, 3);
    
            % calculating access area
            val_view = dot(e_POI_ij, e_sat_repmat, 3);
            access_area_mask = val_view >= cos(beta_max);
            % calculating illuminating conditions
            val_eclipse = dot(e_POI_ij, e_POI2sun, 3);
    
%             eclipse_mask = val_eclipse < 0;
            illumination_constraint_mask = val_eclipse < -sin(abs(theta_sun_max));

            mask = (0.2).^(0.5*access_area_mask + 0.3*illumination_constraint_mask);
            mask = flip(mask); % the mask required to be mirrored in order to fit the image coordinates
            set(I, 'alphadata', mask);
            
            t_local = i*dt;
            logic = t_local > schedule(:,3) & t_local < schedule(:,4);
            index = find(logic);
            if ~isempty(index)
                vec = [[azimuth_index(i-i_local); elevation_index(i-i_local)],[schedule(index,2); 500-schedule(index,1)]];
                los = plot(vec(1,:), vec(2,:), '-k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 1);
                sat_point = plot(azimuth_index(i-i_local), elevation_index(i-i_local), 'ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
                hold on;
                POI_point = plot(schedule(index,2), 500-schedule(index,1) , 'sg', 'MarkerEdgeColor', 'g', 'LineWidth',2, 'MarkerSize', 7);
                hold on;
                drawnow;
                gif;
                delete(sat_point);
                delete(POI_point);  
                delete(los);
            else
                sat_point = plot(azimuth_index(i-i_local), elevation_index(i-i_local), 'ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
                hold on;
                drawnow;
                gif;
                delete(sat_point);                
            end
            index = [];

        end
        delete(GT);
        i_local = i_local + T_in_steps;

    end 
end