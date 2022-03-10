pic1 = load('olympic_rings.txt');

% pic1(:,2) = -pic1(:,2);
% dl = [0, 6];
% pic1 = dl + pic1;
figure;
plot(pic1(:,1), pic1(:,2), 'ok');
axis equal;
axis off;
% save('C:\GoogleDrive\SatelliteFormationFlying\olympic_rings', 'pic1');

pic2 = load('eiffel_tower.txt');
figure;
plot(pic2(:,1), pic2(:,2), 'ok');
axis equal;
axis off;
