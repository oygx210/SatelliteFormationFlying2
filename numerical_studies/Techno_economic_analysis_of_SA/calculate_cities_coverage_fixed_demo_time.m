function schedule = calculate_cities_coverage_fixed_demo_time(r_ECEF, t_vec_JD, theta_sat_min, theta_sun_max, demo_duration, cost_matrix, consts)
 
    dt = round((t_vec_JD(2) - t_vec_JD(1))*consts.day2sec);
    
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

    % Generating POIs grid
    r_orbit = vecnorm(r_ECEF(:,1));
    r_POI_ij = consts.rEarth*cat(3, cTH.*cPHI, cTH.*sPHI, sTH); % POIs coordinates
    
    % calculating critical angles
    b = consts.rEarth/r_orbit;
    gamma_max = asin(b*cos(theta_sat_min));
    beta_max = pi/2 - gamma_max - theta_sat_min;
    t_demo_max = round(2 * beta_max * sqrt(vecnorm(r_ECEF(:,1))^3 / consts.muEarth));
    t_demo_max = t_demo_max - mod(t_demo_max,dt) + dt;
    cost_mask = cost_matrix > 0;
    demonstration_duration_matrix = zeros(size(cost_mask));
    schedule = [];
    schedule_local = [];
    
    f = waitbar(0, 'Demonstration schedule calculation');
    t_global = 0;
    t_simulations = (t_vec_JD(end) - t_vec_JD(1))*consts.day2sec;
    i_global = 0;
    i_local = 0;
    t_local = 0; % to see whether max demonstration time has reached

    attitude_maneuver_time = 30; % seconds
    cumulative_cost = 0;

    % big cycle for the whole demonstration mission
    while t_global < t_simulations

        waitbar(t_global / t_simulations);
        t_local = 0; % to avoid repeating demonstrations at the same POI
        penetration_mask = ones(length(th_range),length(phi_range)); % to avoid repeating demonstrations at the same POI

        % middle cycle for the max demonstration time to avoid repeating
        % dem at the same point
        while t_local < 2*t_demo_max

            penetration_mask_local = ones(length(th_range),length(phi_range));
            
            if (i_global + 2*demo_duration/dt) > length(t_vec_JD)
                return
            else
                for i_local = (i_global+1):(i_global + 2*demo_duration/dt)

                    r_sat = reshape(r_ECEF(:,i_local), [1,1,3]);
                    r_sat_repmat = repmat(r_sat, size(PHI));

                    r_sun_ECI = sun(t_vec_JD(i_local))*consts.AstronomicUnit;
                    r_sun_ECEF = ECItoECEF(t_vec_JD(i_local),r_sun_ECI');
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

                    demonstration_condition = cost_mask & access_area_mask & illumination_constraint_mask & penetration_mask_local & penetration_mask;
                    demonstration_duration_matrix = demonstration_duration_matrix .* demonstration_condition + demonstration_condition*dt;
%                     disp(sum(sum(demonstration_condition)));

                    if sum(sum(demonstration_duration_matrix >= demo_duration)) ~= 0
                        condition = demonstration_duration_matrix >= demo_duration;
                        demonstration_duration_matrix = demonstration_duration_matrix .* condition;
                        [index1,index2, values] = find(demonstration_duration_matrix);
                        schedule_local = [schedule_local;[index1,index2, i_global*ones(length(index1),1)*dt,  i_global*ones(length(index1),1)*dt + values]];
                        penetration_mask_local(index1,index2) = 0;
                    end   
                end

                t_local = t_local + 2*demo_duration;

                if isempty(schedule_local)
                    i_global = i_global + 2*demo_duration/dt;
               else
                    local_max_cost = 0;
                    for i = 1:size(schedule_local,1)
                        index1 = schedule_local(i,1);
                        index2 = schedule_local(i,2);
                        if cost_matrix(index1, index2) > local_max_cost
                            local_max_cost = cost_matrix(index1, index2);
                            index_best = i;                
                       end
                    end

                    schedule = [schedule; [schedule_local(index_best,:), local_max_cost]];
                    penetration_mask(schedule(end,1), schedule(end,2)) = 0;
                    i_global = schedule(end, 4)/dt + round(attitude_maneuver_time/dt);
                    schedule_local = [];
                    demonstration_duration_matrix = zeros(size(cost_mask));
                    cumulative_cost = cumulative_cost + local_max_cost;
%                     disp(cumulative_cost / 1e6);
                end
            end
        end
        if (i_global + 2*demo_duration/dt) > length(t_vec_JD)
            t_global = t_simulations;
        else
            t_global = i_global*dt;
        end 
            
    end
    
    close(f);
end
