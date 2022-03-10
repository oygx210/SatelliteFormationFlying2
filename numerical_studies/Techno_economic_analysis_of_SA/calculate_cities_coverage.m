function [schedule] = calculate_cities_coverage(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, demo_duration_min, population_matrix, consts)
 
    % Purpose: The covered cities are the cities located in the access area of a
    % satellite
    % ref: Biktimirov, et al., Techno-economic analysis of satellite
    % formation-flying missions for space advertising, Journal of Space
    % Systems Modelling (TBD)
    
    % INPUT: 
    % r_ECEF - satellite position in ECEF reference frame given in m, m/s
    % t_vec - array of time corresponding to r_ECEF
    % theta_sat_min - min elevation of satellite required for observation, rad
    % theta_sun_max - max elevation of Sun at POI for demonstration, rad
    % population_matrix - is matrix of size NxM with population assigned to
    % node

    % OUTPUT
    % Schedule - list of possible demonstrations [citi_i, t_start_i, duration_i]

    demo_flag = 'YES'; % YES | NO
    earth = imread('earth_1024_512.jpg');
    dt = round((t_vec_JD(2) - t_vec_JD(1))*consts.day2sec);
    
    % sphere generation
    [M, N] = size(population_matrix);
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
    
    population_mask = population_matrix > 0;
    demonstration_duration_matrix = zeros(size(population_mask));
    demonstration_opportunities = [];
    schedule = [];
    switch demo_flag
        case 'YES'
            figure('Name', 'Earth coverage','NumberTitle', 'off');
            I = imshow(earth);
            hold on;
    end
    
    f = waitbar(0, 'Coverage calculation');
    
    for i = 1:length(t_vec_JD)

        waitbar(i/length(t_vec_JD));
        
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
        illumination_constraint_mask = val_eclipse < -sin(abs(theta_sun_max));

        demonstration_condition = population_mask & access_area_mask & illumination_constraint_mask;
        
        switch demo_flag
            case 'YES'
                eclipse_mask = val_eclipse < 0;
                mask = (0.2).^(0.5*access_area_mask + 0.3*eclipse_mask);
                mask = flip(mask); % the mask required to be mirrored in order to fit the image coordinates
                set(I, 'alphadata', mask);
                drawnow;
        end
        
        condition = demonstration_duration_matrix >= demo_duration_min & ~demonstration_condition;
        condition = demonstration_duration_matrix .* condition;
        
        if sum(condition, 'all') ~= 0
            [index1,index2, values] = find(condition);
            schedule = [schedule;[index1,index2, t_vec_JD(i) - values,  values]];
            disp([index1,index2, t_vec_JD(i) - values,  values]);
        end
        
        demonstration_duration_matrix = demonstration_duration_matrix .* demonstration_condition + demonstration_condition*dt;
        disp(sum(sum(demonstration_condition)));
    end
    close(f);  
end