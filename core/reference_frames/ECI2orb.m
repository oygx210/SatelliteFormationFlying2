function rv_chaser_orb = ECI2orb(rv_target_ECI, rv_chaser_ECI, consts)

% x - along track; y - out-of-plane; z - radial
oe_target = rv2oe(rv_target_ECI, consts);
mean_motion = sqrt(consts.muEarth/oe_target(1,:)^3);

% mean_motion = sqrt(consts.muEarth/norm(rv_target_ECI(1:3))^3);

z_orb = rv_target_ECI(1:3) / norm(rv_target_ECI(1:3));
y_orb = cross(rv_target_ECI(1:3), rv_target_ECI(4:6));
y_orb = y_orb(1:3) / norm(y_orb(1:3));
x_orb = cross(y_orb, z_orb);

% Conversion matrix (ECI2orb)
ECI2orb = [x_orb y_orb z_orb]';

r_chaser_orb = ECI2orb * (rv_chaser_ECI(1:3) - rv_target_ECI(1:3));
v_chaser_orb = ECI2orb * (rv_chaser_ECI(4:6) - rv_target_ECI(4:6)) -  cross([0;mean_motion;0],r_chaser_orb);

rv_chaser_orb = [r_chaser_orb ; v_chaser_orb];
end