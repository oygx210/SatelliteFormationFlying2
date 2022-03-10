function rv_chaser_orb = vecECI2orb(rv_target_ECI, rv_chaser_ECI, consts)

% x - along track; y - out-of-plane; z - radial

if size(rv_target_ECI,2) > 1
    a = size(rv_target_ECI,2);
    oe_target = vecRV2OE(rv_target_ECI, consts);
    mean_motion = sqrt(consts.muEarth./oe_target(1,:).^3);

    z_orb = rv_target_ECI(1:3,:) ./ vecnorm(rv_target_ECI(1:3,:));
    y_orb = cross(rv_target_ECI(1:3,:), rv_target_ECI(4:6,:));
    y_orb = y_orb(1:3,:) ./ vecnorm(y_orb(1:3,:));
    x_orb = cross(y_orb, z_orb);
    ECI2orb = zeros(3,3,a);

    ECI2orb(1,:,:) = x_orb;
    ECI2orb(2,:,:) = y_orb;
    ECI2orb(3,:,:) = z_orb;

    r_chaser_ECI = reshape(rv_chaser_ECI(1:3,:), [3,1,a]);
    r_target_ECI = reshape(rv_target_ECI(1:3,:), [3,1,a]);
    r_diff = r_chaser_ECI - r_target_ECI;
    
    v_chaser_ECI = reshape(rv_chaser_ECI(4:6,:), [3,1,a]);
    v_target_ECI = reshape(rv_target_ECI(4:6,:), [3,1,a]);
    v_diff = v_chaser_ECI - v_target_ECI;

    r_chaser_orb = zeros(3,a);
    v_chaser_orb = zeros(3,a);

    for i = 1:a
    r_chaser_orb(:,i) = ECI2orb(1:3,1:3,i) * r_diff(1:3,1,i);
    v_chaser_orb(:,i) = ECI2orb(1:3,1:3,i) * v_diff(1:3,1,i) -  cross([0;mean_motion(1,i);0],r_chaser_orb(:,i));
    end
    
    rv_chaser_orb = [r_chaser_orb ; v_chaser_orb];
 
else
    b = size(rv_chaser_ECI,2);
    
    oe_target = vecRV2OE(rv_target_ECI, consts);
    mean_motion = sqrt(consts.muEarth/oe_target(1,:)^3);

    z_orb = rv_target_ECI(1:3) / norm(rv_target_ECI(1:3));
    y_orb = cross(rv_target_ECI(1:3), rv_target_ECI(4:6));
    y_orb = y_orb(1:3) / norm(y_orb(1:3));
    x_orb = cross(y_orb, z_orb);

    % Conversion matrix (ECI2orb)
    ECI2orb = [x_orb y_orb z_orb]';

    r_chaser_orb = ECI2orb * (rv_chaser_ECI(1:3,:) - rv_target_ECI(1:3));
    v_chaser_orb = ECI2orb * (rv_chaser_ECI(4:6,:) - rv_target_ECI(4:6)) -  cross([0;mean_motion;0].*ones(1,b),r_chaser_orb(1:3,:));

    rv_chaser_orb = [r_chaser_orb ; v_chaser_orb];
end

end