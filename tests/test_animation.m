function test_animation(rv_ECI_simulations, rv_ECI_target_orbit, VideoHeader)
    fig = figure;
    step = 1;
    for i = 1:size(rv_ECI_simulations,2)/step
        a = plot3(rv_ECI_target_orbit(1,i*step), rv_ECI_target_orbit(2,i*step), rv_ECI_target_orbit(3,i*step), 'sk', 'MarkerSize', 4);
        b = plot3(rv_ECI_simulations(1,i*step),rv_ECI_simulations(2,i*step), rv_ECI_simulations(3,i*step), 'sr', 'MarkerSize', 4);

        drawnow;

    %   Take a Snapshot
        movieVector(i) = getframe(fig);   %manually specify getframe region    

        delete(a);
        delete(b);
    end
      
    myWriter = VideoWriter(VideoHeader,'MPEG-4');   %create an .mp4 file
    myWriter.FrameRate = 24; %

    %   Open the VideoWriter object, write the movie, and close the file
    open(myWriter);
    writeVideo(myWriter, movieVector);
    close(myWriter); 
end   