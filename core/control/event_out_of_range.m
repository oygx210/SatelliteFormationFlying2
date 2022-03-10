function [value_out,isterminal,direction] = event_out_of_range(t,chaser_rv,target_rv_array,t_span,rho_init,alpha_init,n,range)

% event_out_of_range sends event when satellite is out of desired range.

target_rv = zeros(6,1);

for j=1:6
    target_rv(j) = interp1(t_span,target_rv_array(j,:),t,'spline');
end

target_oe = rv2oe(target_rv(1:3), target_rv(4:6));
arg_of_lat = target_oe(5) + target_oe(6);

chaser_mod = ECI2orb(target_rv,chaser_rv,n);

c1 = rho_init;
c2 = sqrt(3)/2* rho_init;
c3 = 0;

chaser_predicted= get_rv_from_analytic_HCW_solution(c1, c2, c3, alpha_init, arg_of_lat, n);

rel_distance = norm(chaser_mod(1:3) - chaser_predicted(1:3))/norm(chaser_predicted(1:3));
value = range - rel_distance;

isterminal = 1;
direction = 1;

if value < 0
    value_out = 0;
else
    value_out = 1;
end

end


