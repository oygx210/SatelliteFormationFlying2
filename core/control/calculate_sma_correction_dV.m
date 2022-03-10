function [t_span, dV_sma] = calculate_sma_correction_dV(rv_ECI_current, rv_ECI_required, consts)

    % the function calculated dV to adjust sma
    
    moe = rv2moe(rv_ECI_current, consts);
    coe = rv2oe(rv_ECI_current, consts);

    moe_required = rv2moe(rv_ECI_required, consts);
    % coe_required = rv2oe(rv_ECI_required, consts);
    delta_moe = moe_required - moe;
%     disp(delta_moe);
    n = sqrt(consts.muEarth/moe(1)^3);
    etta = sqrt(1-moe(2)^2);
    da = delta_moe(1);

    if coe(6) < pi 
        t_span = (pi - coe(6))*sqrt(moe(1)^3/consts.muEarth);
    elseif coe(6) > pi
        t_span = (2*pi - coe(6))*sqrt(moe(1)^3/consts.muEarth);
    end

    if abs(da) > 5
        dV_sma = [n*da*etta/2;
                  0;
                  0];
    else
        dV_sma = [0;
                  0;
                  0];
    end
end