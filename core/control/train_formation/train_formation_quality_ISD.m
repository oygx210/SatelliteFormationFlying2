function [value_out,isterminal,direction] = train_formation_quality_ISD(t, rv, consts, spacecraft, formation)

isterminal = 1;
direction = 0;

if abs(vecnorm(rv(1:3) - rv(7:9)) - formation.ISD) > formation.ISD_acceptable_error 
    value_out = 0;
else
    value_out = 1;
end

end