function make_movie(chaser)

figh = figure();
% figh.Color = 'k';
% ax = gca;
% c = ax.Color;
% ax.Color = 'black';
% axis off;
% camlight
% set(gca,'color','k')
% set(gcf,'color','k')
axis tight manual;

for k = 1:size(chaser(1).rv_orb,2)-2

    clf
    for i = 1:31
        x_k(i) = chaser(i).rv_orb(1,k);
        y_k(i) = chaser(i).rv_orb(2,k);
        z_k(i) = chaser(i).rv_orb(3,k);
    end    
    
    plot3(x_k, y_k, z_k, 'sk', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
%     ax = gca;
%     c = ax.Color;
% %     ax.Color = 'black';
%     axis off;

    xlim([-4000 4000]);
    ylim([-4000 4000]);
    zlim([-4000 4000]);
    view([1 -1 -2]);

    %Add plotting options
%     grid on;
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title(['t = ',num2str(t_k)])
%     view([30 35])
%     view([30+20*t_k 35])      %show how viewpoint can be manipulated
    
    %% Step 3: Take a Snapshot
    %Save the frame
%   movieVector(k) = getframe;
    movieVector(k) = getframe(figh, [50 0 400 500]);   %manually specify getframe region    
    %% Step 4: Advance Time
    %Happens automatically if using a for loop
end

%% Step 5: Save Movie
%Create a VideoWriter object and set properties
% myWriter = VideoWriter('curve', );            %create an .avi file

myWriter = VideoWriter('form','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 24;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector(156:end));
close(myWriter);


disp('DONE!')