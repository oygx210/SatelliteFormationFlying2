function formation_flying_animation_3subplots(t_vec, rv_orb, rv_orb_required, build_reference_trajectories_flag, VideoHeader)

    fig = figure;
    subplot(1,3,1);
    plot3(rv_orb(1,1,1), rv_orb(2,1,1), rv_orb(3,1,1), 'ok', 'MarkerSize', 1);
    xlabel('tangential,m');
    ylabel('normal,m');
    zlabel('radial,m');
    xlim([-5000 5000]);
    ylim([-5000 5000]);
    zlim([-5000 5000]);
    pbaspect([1 1 1]);
    title('Orbital reference frame');
    view(-45,-75);
    hold on;
    subplot(1,3,3);
    plot3(rv_orb(1,1,1), rv_orb(2,1,1), rv_orb(3,1,1), 'ok', 'MarkerSize', 1);
    xlabel('tangential,m');
    ylabel('normal,m');
    xlim([-5000 5000]);
    ylim([-5000 5000]);
    hold on;
    axis square;
    title('Image projection to the horizontal plane');
    view(0,-90);
    step = 40;
    subplot(1,3,2);
    [X,Y,Z]=sphere(50);
    R=6378e3;
    globe= surf(-X*R,Y*R,-Z*R);
    lat = deg2rad(48.856613);
    lon = deg2rad(2.352222);
    POI = zeros(3,1);
    [POI(1),POI(2), POI(3)] = sph2cart(lon,lat,R+1000);
    cdata = imread('earthmap1k.jpg');
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata,  'EdgeColor', 'none');
    set(gcf,'Color','w')
    set(gca, 'visible', 'off')
    axis equal
    ax.Clipping = 'off';
    view (90,0);
    hold on;
    P = plot3(POI(1), POI(2), POI(3), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    T = 86400;
    legend('Earth','POI');
    pbaspect([1 1 1]);
    
    switch build_reference_trajectories_flag
        case 'on'
            for i = 1:size(rv_orb,3)/step
                subplot(1,3,1);
                a = plot3(rv_orb(1,2:end,i*step), rv_orb(2,2:end,i*step), rv_orb(3,2:end,i*step), 'sk', 'MarkerSize', 4);
                b = plot3(rv_orb_required(1,2:end,i*step), rv_orb_required(2,2:end,i*step), rv_orb_required(3,2:end,i*step), '+r', 'MarkerSize', 2);

                subplot(1,3,3);
                c = plot3(-rv_orb(1,2:end,i*step),-rv_orb(2,2:end,i*step), rv_orb(3,2:end,i*step), 'sk', 'MarkerSize', 4);
                d = plot3(-rv_orb_required(1,2:end,i*step), -rv_orb_required(2,2:end,i*step), rv_orb_required(3,2:end,i*step), '+r', 'MarkerSize', 2);

                subplot(1,3,2);
                rotate(globe, [0 0 1], 1);
                rotate(P, [0 0 1], 1);

                drawnow;
                info = string(['T = ', num2str(t_vec(i*step)), ' min']);
                an = annotation('textbox',[0.4 0.0 0.1 0.1], 'String', info); 

            %   Take a Snapshot
                movieVector(i) = getframe(fig);   %manually specify getframe region    

                delete(a);
                delete(b);
                delete(c);
                delete(d);                
                delete(an);
            end
            
        case 'off'
             for i = 1:size(rv_orb,3)/step
                subplot(1,2,1);
                a = plot3(rv_orb(1,2:end,i*step), rv_orb(2,2:end,i*step), rv_orb(3,2:end,i*step), 'sk', 'MarkerSize', 4);
%                 legend('Orbital reference frame origin','Current configuration', 'Required configuration');
                subplot(1,2,2);
                c = plot3(-rv_orb(1,2:end,i*step),-rv_orb(2,2:end,i*step), rv_orb(3,2:end,i*step), 'sk', 'MarkerSize', 4);
%                 legend('Orbital reference frame origin','Current configuration (x-y plane)', 'Required configuration (x-y plane)');

                drawnow;
                info = string(['T = ', num2str(t_vec(i*step)), ' min']);
                an = annotation('textbox',[0.4 0.0 0.1 0.1], 'String', info); 

            %   Take a Snapshot
                movieVector(i) = getframe(fig);   %manually specify getframe region    

                delete(a);
                delete(c);
                delete(an);
            end
    end
    
            myWriter = VideoWriter(VideoHeader,'MPEG-4');   %create an .mp4 file
            myWriter.FrameRate = 24;

            %   Open the VideoWriter object, write the movie, and close the file
            open(myWriter);
            writeVideo(myWriter, movieVector);
            close(myWriter);
    
end
   