function [collision_flag, ISD_min] = collision_check(t, rv, formation)
    mode = 1;
    ISD = [];
    ISD_min = 10000;
    rv = rv(7:end,:);
    
    for i = 1:formation.N_active_sats
        for j = 1:formation.N_active_sats
            if i~=j
               ISD_local = vecnorm(rv(i*6-5:i*6-3,:) - rv(j*6-5:j*6-3,:));
               ISD_min_local = min(ISD_local);
               if ISD_min_local < ISD_min
                   ISD_min = ISD_min_local;
               end
               ISD = [ISD; ISD_local];
            end
        end
    end
    
    if ISD_min > formation.ISD_safe
        collision_flag = 0;
    else
        collision_flag = 1;
        warning('Potential collision is detected');

    end
    if mode == 2
        figure;
        yline(formation.ISD_safe, 'k--');
        xlabel('time, min');
        ylabel('ISD, m');
        grid on;
        hold on;
        for i = 1:size(ISD,1)
            plot(t/60, ISD(i,:));
            hold on;
        end
        legend('ISD safe', 'ISD(t) between satellites i and j');
    end
end