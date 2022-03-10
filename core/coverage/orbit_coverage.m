function tau_max = orbit_coverage(target_orbit_oe, studied_orbit_oe, orbit_epoch, n_nodes, spacecraft, consts)

    options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
    T_local = 0;
    u_vec = 0 : 2*pi/n_nodes : (2*pi - 2*pi/n_nodes);
    
    for i = 1:n_nodes
        oe = studied_orbit_oe;
        oe(6) = u_vec(i);
        rv_studied_orbit_grid(6*i-5:6*i,1) = oe2rv(oe, consts);
    end
    dt = 100;

    rv_initial = [oe2rv(target_orbit_oe, consts); rv_studied_orbit_grid];
    
    convergence = 0;
    n_nodes_local = n_nodes;
    while convergence ~= 1
        
        R_sun = sun(juliandate(orbit_epoch + seconds(T_local)))'*consts.AstronomicUnit;
        
        [t_vec, rv_ECI_local] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft), [T_local (T_local+dt)], rv_initial, options_precision);
        rv_ECI_local = rv_ECI_local'; 
        T_local = t_vec(end);
        disp(n_nodes_local);
        covered_nodes = false(n_nodes_local,1);
        drv_matrix = rv_ECI_local(7:end, :) - repmat(rv_ECI_local(1:6,:), n_nodes_local,1);
        
        for i = 1:n_nodes_local
            dr = vecnorm(drv_matrix(6*i-5:6*i-3,:)) <= spacecraft.dr_observation;
            for j = 1:length(t_vec)
                los(j) = sight(rv_ECI_local(6*(i+1)-5:6*(i+1)-3,j), R_sun, 's');
            end
            covered_nodes(i,1) = sum(dr(:) == 1 & los(:) == 1) > 0;
        end
        
        rv_initial = rv_ECI_local(:,end);
        
        covered_nodes = (repelem([false;covered_nodes] ,6));
        rv_initial(covered_nodes) = [];
        n_nodes_local = size(rv_initial,1)/6 - 1;
        los = [];

        if n_nodes_local == 0
            convergence = 1;
        end
    end
    
    tau_max = t_vec(end);
    
end
