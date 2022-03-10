function [dVmag] = simple_control(consts, rv_ECI_desired, rv_ECI_current)

% target is position to be reached
% current is satellite position at the moment

% Inputs are column vectors representing target and current satellite given
% in ECI, m; m/s

desired_oe = rv2oe(rv_ECI_desired(1:3), rv_ECI_desired(4:6));
current_oe = rv2oe(rv_ECI_current(1:3), rv_ECI_current(4:6));

target_oe_eqinoctial = [desired_oe(1); desired_oe(2)*cos(desired_oe(4)); desired_oe(2)*sin(desired_oe(4)); desired_oe(3); desired_oe(5); desired_oe(6)]; 
current_oe_eqinoctial = [current_oe(1); current_oe(2)*cos(current_oe(4)); current_oe(2)*sin(current_oe(4)); current_oe(3); current_oe(5); current_oe(6)]; 

oe_diff = current_oe_eqinoctial - target_oe_eqinoctial;

dVmag = calculate_dv_for_estimates(consts, desired_oe(1), desired_oe(3), oe_diff(2), oe_diff(3), oe_diff(4), oe_diff(5));

end