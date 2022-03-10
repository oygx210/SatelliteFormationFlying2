function rv_prime = rhs_Formation_HCW(t, rv_orb, inertial_orbit, consts)

    % rhs for N satellites' relative motion dynamics 
    % 'point mass' - motion in central gravity field
    % 'J2' - taking into account J2 effect
    
    global environment
  
    N_sats = length(rv_orb)/6;
    rv_orb = reshape(rv_orb, [6,N_sats]);

    r_prime = [rv_orb(4:6,:); zeros(3,N_sats)];
    r_norm = vecnorm(rv_orb(1:3,:));

    switch environment 
        case 'point mass'
            
    rv_prime = [rv_orb(4); rv_orb(5); rv_orb(6);
                -2 * mean_motion * rv_orb(6);
                -mean_motion^2 * rv_orb(2);
                2 * mean_motion * rv_orb(4) + 3 * mean_motion^2 * rv_orb(3)];

end