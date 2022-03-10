function schedule = schedule_burns(rv_target_ECI, consts, formation)

% schedule_burns returns times when both impulses are delivered
% we spend 2 orbits at maximum for maneuvering
% the manevering starts when formation at the specific argument of latitude

alpha = formation.geometry(4,:);
oe_target = vecRV2OE(rv_target_ECI, consts);
arg_of_lat = oe_target(8);
period = 2*pi*sqrt(oe_target(1)^3/consts.muEarth);

gamma = 2*pi - alpha;
theta = arg_of_lat;

    if gamma > theta
        delta_theta_1 = gamma - theta;
        delta_theta_2 = delta_theta_1 + pi;
        delta_theta_3 = 6*pi*ones(1, formation.N_sats);
        
    else
        delta_theta_1 = 2*pi - gamma + theta;
        delta_theta_2 = delta_theta_1 + pi;
        delta_theta_3 = 6*pi*ones(1, formation.N_sats);
    end
    
    delta_theta = [delta_theta_1 ; delta_theta_2 ; delta_theta_3];
    
    schedule = delta_theta/2/pi .*period;

end