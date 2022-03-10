%% Undocking phase

% Experiment_option: 1 - undocking with random deployment errors and uncontrolled drifting
%                    2 - the satellites have the same sma (idealized case)

rv_ECI = [];
rv_orb = [];
t_vec = [];

Undocking_option = 2;

switch Undocking_option
    case 1   
        % Undocking & uncontrolled flying

        disp('Undocking');
        [t_undocking, rv_ECI_undocking] = undocking(formation.rv, deployer, formation, spacecraft, consts);
        rv_ECI = [rv_ECI, rv_ECI_undocking];
        t_vec = [t_vec; t_undocking];
        T = t_vec(end);
        disp([num2str(formation.N_sats),' satellites were successfully released']);
        % look at the ISD
        
        for i = 1:formation.N_sats
            oe_undocking(:,i) = rv2oe(rv_ECI_undocking(i*6-5:i*6,end), consts);
        end
                        
        disp('Uncontrolled flying');
        options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [t_vec_drifting, rv_ECI_drifting] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),0:formation.drifting_time, rv_ECI_undocking(:,end), options_precision);        
        rv_ECI_drifting = rv_ECI_drifting';
        rv_ECI = [rv_ECI, rv_ECI_drifting];
        t_vec = [t_vec; t_vec_drifting];
        T = t_vec(end);

        for i = 1:formation.N_sats
            rv_orb(:,i) = ECI2orb(rv_ECI_drifting(1:6,end), rv_ECI_drifting(i*6-5:i*6,end), consts);
        end
        
        figure;
        plot3(rv_orb(1,2:end), rv_orb(2,2:end),rv_orb(3,2:end), 'ok', 0,0,0,'or');
        title('Orbital configuration after commissioning period');
        xlabel('tangential,m');
        ylabel('normal,m');
        zlabel('radial,m');
        axis equal;
        legend('active satellites', 'virtual satellites');

    case 2        

        for i = 1:formation.N_sats
            rv_ECI_drifting(i*6-5:i*6,1) = formation.rv;
        end
        rv_ECI = [rv_ECI, rv_ECI_drifting];
        disp([num2str(formation.N_sats-1),' active satellites and 1 virtual satellite are initialized']);

end
