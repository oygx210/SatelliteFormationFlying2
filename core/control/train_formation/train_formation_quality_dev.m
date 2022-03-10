function [value_out,isterminal,direction] = train_formation_quality(t, rv, consts, spacecraft, formation)

global formation_condition

leader_coe = rv2coe(rv(1:6), consts);
leader_mean_oe = osc2mean(leader_coe, consts);

formation_orbit = abs(leader_coe(1) - formation.target_sma) > formation.max_sma_error;
formation_quality = abs(vecnorm(rv(1:3) - rv(7:9)) - formation.ISD) > formation.ISD_acceptable_error;

if formation_orbit == 1
    formation_condition = 2;
elseif formation_quality == 1
    formation_condition = 1;
elseif formation_quality == 0 && formation_orbit == 0
    formation_condition = 0;
end

isterminal = 1;
direction = 0;

if  formation_condition == 0
    value_out = 1;
else
    value_out = 0;
end

end