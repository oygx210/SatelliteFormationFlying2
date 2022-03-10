function moe = rv2moe(rv_ECI, consts)

% Function converts rv to mean orbital elements

coe = rv2oe(rv_ECI, consts);
for i = 1:size(coe,2)
    moe(:,i) = osc2mean(coe(:,i), consts);
end

end
