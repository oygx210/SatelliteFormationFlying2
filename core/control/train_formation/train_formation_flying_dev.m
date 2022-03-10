function [t_vec, rv_ECI] = train_formation_flying(rv_ECI, simulation_time, consts, spacecraft, formation)

disp('Formation Flying');
global formation_condition

T_local = 0;
t_vec = 0;

% Event function
tracking_formation_quality = @(t, rv) train_formation_quality(t, rv, consts, spacecraft, formation);
option_tracking = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', tracking_formation_quality);

while T_local < simulation_time
  
    train_formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);

    while formation_condition ~= 0
    
        if formation_condition == 1
            disp(['Flag-', num2str(formation_condition)]);
            [t_vec_m, rv_ECI_m] = formation_configuration_maintenance(rv_ECI(:,end), consts, spacecraft, formation);
            t_vec = [t_vec; t_vec_m + t_vec(end)];
            rv_ECI = [rv_ECI, rv_ECI_m];

            T_local = t_vec(end);
            train_formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);

        elseif formation_condition == 2
            disp(['Flag-', num2str(formation_condition)]);
            [t_vec_m, rv_ECI_m] = formation_orbit_maintenance(rv_ECI(:,end), consts, spacecraft, formation);
            t_vec = [t_vec; t_vec_m + t_vec(end)];
            rv_ECI = [rv_ECI, rv_ECI_m];

            T_local = t_vec(end);
            train_formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
        end  
    end
    
%     train_formation_quality(T_local, rv_ECI(:,end), consts, spacecraft, formation);
    disp(['Flag-', num2str(formation_condition)]);
    rv = rv_ECI(:,end);
    [t_vec_FF, rv_ECI_FF, te, ye, ie] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), T_local:simulation_time, rv, option_tracking);
    rv_ECI_FF = rv_ECI_FF';
    rv_ECI = [rv_ECI, rv_ECI_FF];
    t_vec = [t_vec; t_vec_FF];
    
    T_local = t_vec(end);
    
end
disp('Successful flight');





