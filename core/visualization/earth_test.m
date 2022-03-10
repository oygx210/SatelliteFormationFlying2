clear all;

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
for i = 1:360

rotate(globe, [0 0 1], 1);
rotate(P, [0 0 1], 1);

drawnow
end