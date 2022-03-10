function plot_dynamics(rv_ECI, formation, consts)

    for i = 1:size(rv_ECI,2)
        rv_orb_required(:,:,i) = get_rv_from_analytic_HCW_solution(rv_ECI(1:6,i), formation.geometry, consts);
        
        for j = formation.N_sats
            rv_orb(:,j,i) = ECI2orb(rv_ECI(1:6,i),rv_ECI(j*6-5:j*6,i), consts);
        end
    end

    figure;
    subplot(1,2,1);
    plot3(rv_orb_required(1,:,end), rv_orb_required(2,:,end), rv_orb_required(3,:,end), 'ok');
    xlabel('x, m', 'fontsize',12);
    ylabel('y, m', 'fontsize',12);
    zlabel('z, m', 'fontsize',12);
    xlim([-1000 1000]);
    ylim([-1000 1000]);
    zlim([-1000 1000]);
    pbaspect([1 1 1]);
    title('Required configuration');
    view(-225,45);

    subplot(1,2,2);
    plot3(rv_orb(1,2,end), rv_orb(2,2,end), rv_orb(3,2,end), 'ok');
    xlabel('x, m', 'fontsize',12);
    ylabel('y, m', 'fontsize',12);
    zlabel('z, m', 'fontsize',12);
    xlim([-1000 1000]);
    ylim([-1000 1000]);
    zlim([-1000 1000]);
    pbaspect([1 1 1]);
    title('Required configuration');
    view(-225,45);

    for i = 1:size(rv_ECI,2)
        subplot(1,2,1);
        plot3(rv_orb_required(1,:,i), rv_orb_required(2,:,i), rv_orb_required(3,:,i), 'ok');
        subplot(1,2,2);
        plot3(rv_orb(1,1,i), rv_orb(2,1,i), rv_orb(3,1,i), 'ok');

        pause;
    end
    
end