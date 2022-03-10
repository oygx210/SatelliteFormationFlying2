function [t_vec, rv_ECI, maneuvers_dV, maneuvers_t] = train_formation_flying_ISD_keeping(rv_ECI, simulation_time, consts, spacecraft, formation)

T_local = 0;
t_vec = 0;
maneuvers_dV = zeros(9,1);
maneuvers_t = zeros(3,1);

% Event function
tracking_formation_quality = @(t, rv) train_formation_quality_ISD(t, rv, consts, spacecraft, formation);
options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', tracking_formation_quality);

quality = train_formation_quality_ISD(T_local, rv_ECI, consts, spacecraft, formation);

if quality == 0 
    disp('Orbit correction');
    [t_vec_m, rv_ECI_m, maneuvers_dV, maneuvers_t] = formation_maintenance_3impulses_Mok_moe_based(rv_ECI, consts, spacecraft, formation);
    t_vec = [t_vec; t_vec_m];
    rv_ECI = [rv_ECI, rv_ECI_m];

    T_local = t_vec(end);
end

quality_check = train_formation_quality_ISD(T_local, rv_ECI(:,end), consts, spacecraft, formation);
if quality_check == 0
    disp('Error 1');
end


while T_local < simulation_time
    
    disp('Formation flying');
    rv = rv_ECI(:,end);
    [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), T_local:simulation_time, rv, options);
    rv_ECI_FF = rv_ECI_FF';
    rv_ECI = [rv_ECI, rv_ECI_FF];
    t_vec = [t_vec; t_vec_FF];

    quality = train_formation_quality_ISD(T_local, rv_ECI(:,end), consts, spacecraft, formation);
    T_local = t_vec(end);

    if abs(simulation_time - T_local) < 1
        T_local = simulation_time;
    end

    if T_local < simulation_time && quality == 0

        disp(['Orbit correction, elapsed time = ', num2str(T_local / simulation_time * 100), '%']);
        [t_vec_maneuvering, rv_maneuvering, maneuvers_dV, maneuvers_t] = formation_maintenance_3impulses_Mok_moe_based(rv_ECI(:,end), consts, spacecraft, formation);
        rv_ECI = [rv_ECI, rv_maneuvering];
        t_vec = [t_vec; t_vec_maneuvering + t_vec(end)];        
        rv = rv_ECI(:,end);
        T_local = t_vec(end);
        disp(['Total dV for correction = ', num2str(sum(abs(maneuvers_dV))), ' m/s']);

        maneuvers_dV = [maneuvers_dV, maneuvers_dV];
        maneuvers_t = [maneuvers_t, maneuvers_t];
    end
    
    quality_check = train_formation_quality_ISD(T_local, rv_ECI(:,end), consts, spacecraft, formation);
    
    if quality_check == 0
        disp('Error 2');
    end
    
    if T_local >= simulation_time
        disp('Mission complete!');
    end
        
end