function dV = calculate_maneuver(rv_target_ECI, rv_ECI_current, rv_ECI_desired, consts)
a = size(rv_ECI_current,2);
% The function return dV1 and dV2 and time for maneuvering which is
% required to reach desired reference trajectory

% sma - m
% all angles are given in radians

current_oe = vecRV2OE(rv_ECI_current, consts);
desired_oe = vecRV2OE(rv_ECI_desired, consts);
target_oe = vecRV2OE(rv_target_ECI, consts);

a_c = current_oe(1,:);
ecc_c = current_oe(2,:);
inc_c = current_oe(3,:);
Omega_c = current_oe(4,:);
omega_c = current_oe(5,:);
M_c = current_oe(7,:);

a_d = desired_oe(1,:);
ecc_d = desired_oe(2,:);
inc_d = desired_oe(3,:);
Omega_d = desired_oe(4,:);
omega_d = desired_oe(5,:);
M_d = desired_oe(7,:);

% oe_equinoctial = [a, q1, q2, inc, Omega, lambda] = [a, ecc*cos(omega), ecc*sin(omega), inc, Omega, omega + M]

current_oe_equinoctial = [a_c; ecc_c.*cos(omega_c); ecc_c.*sin(omega_c); inc_c; Omega_c; omega_c + M_c]; 
desired_oe_equinoctial = [a_d; ecc_d.*cos(omega_d); ecc_d.*sin(omega_d); inc_d; Omega_d; omega_d + M_d]; 

oe_diff = desired_oe_equinoctial - current_oe_equinoctial;

dinc = oe_diff(4,:);
dOmega = oe_diff(5,:);
dq1 = oe_diff(2,:);
dq2 = oe_diff(3,:);

target_SC_inc = target_oe(3);

p =  1; % basically, ratio between first and second out-of-plane burns. 
% See Vaddi et al., page 264, section called "Aside"
gamma = sqrt(target_oe(1)/consts.muEarth);  %[s/m]

% deltaVs are in m/s
dV_h1 = p/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(target_SC_inc).^2);
dV_h2 = (1-p)/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(target_SC_inc).^2);
dV_r1 = -sqrt(dq1.^2 + dq2.^2)/2/gamma;
dV_r2 = sqrt(dq1.^2 + dq2.^2)/2/gamma;

dV1 = [zeros(1,a) ; dV_h1 ; dV_r1];
dV2 = [zeros(1,a) ; dV_h2 ; dV_r2];
dV = [dV1; dV2];
end