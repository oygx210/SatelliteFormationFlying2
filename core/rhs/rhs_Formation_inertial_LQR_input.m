function rv_prime = rhs_Formation_inertial_LQR_input(T_local, rv, consts, spacecraft, u_matrix)
% rhs for N satellites
% 'point mass' - motion in central gravity field
% 'J2' - taking into account J2 effect
% 'J2 + atmospheric drag' - taking into account J2 effect and atmosperic drag

    global environment

    N_sats = length(rv)/6;
    rv = reshape(rv, [6,N_sats]);
    r_prime = [rv(4:6,:); zeros(3,N_sats)];
    r_norm = vecnorm(rv(1:3,:));

    % Central gravity field
    A_cg = [zeros(3,N_sats); -consts.muEarth./(r_norm.^3).*rv(1:3,:)];
    switch environment
        case 'point mass'
            rv_prime = r_prime + A_cg - u_matrix;
        case 'J2'            
            R = consts.rEarth_equatorial;
            J2 = consts.J2;
            mu = consts.muEarth;
            delta = 3/2*J2*mu*R^2;
    
            % J2 perturbation
            A_j2 = delta./(r_norm.^5).*rv(1:3,:).*(5*rv(3,:).^2./r_norm.^2 - 1) - 2*delta./r_norm.^5.*[zeros(1,N_sats);zeros(1,N_sats);rv(3,:)];
            A_j2 = [zeros(3,N_sats); A_j2];

            rv_prime = r_prime + A_cg + A_j2 - u_matrix;
        case 'J2 + atmospheric drag'
            R = consts.rEarth_equatorial;
            J2 = consts.J2;
            mu = consts.muEarth;
            delta = 3/2*J2*mu*R^2;
            % J2 perturbation
            A_j2 = delta./(r_norm.^5).*rv(1:3,:).*(5*rv(3,:).^2./r_norm.^2 - 1) - 2*delta./r_norm.^5.*[zeros(1,N_sats);zeros(1,N_sats);rv(3,:)];
            A_j2 = [zeros(3,N_sats); A_j2];

            % Atmospheric drag
            vRelativeECI = rv(4:6,:) - cross(ones(1,N_sats).*consts.wEarth, rv(1:3,:));
            rhoAtmo = CIRA72(consts, (vecnorm(rv(1:3,:)) - consts.rEarth));
            A_ad = - 0.5 * spacecraft.Cdrag * spacecraft.DragArea / spacecraft.mass * rhoAtmo .* vRelativeECI .* vecnorm(vRelativeECI); 
            A_ad = [zeros(3,N_sats); A_ad];

            rv_prime = r_prime + A_cg + A_j2 + A_ad;
    end
    
    rv_prime = reshape(rv_prime, [6*N_sats,1]);

end