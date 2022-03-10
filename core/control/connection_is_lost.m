function [value_out,isterminal,direction] = connection_is_lost(t, rv, spacecraft)

% The function controls intersatellite distance and stops integrating ODE
% if current inter-satellite distance exceeds the max ISL

isterminal = 1;
direction = -1;

if spacecraft.max_ISL < norm(rv(1:3) - rv(7:9))
    value_out = 0;
else
    value_out = 1;
end

end


