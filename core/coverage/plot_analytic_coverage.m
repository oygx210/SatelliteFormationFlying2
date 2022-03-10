% plotting
cTH = cos(TH);
sTH = sin(TH);
cPHI = cos(PHI);
sPHI = sin(PHI);

% satellite view parameters
Re = consts.rEarth; % km
b = oe(1)/Re;
th_max = asin(1/b);
th_eff = th_max - th_beam;
gam_max = pi/2 - th_max; % angle along earth surface
gam_eff = -th_eff + asin(b*sin(th_eff));
gam_beam = -th_beam + asin(b*sin(th_beam));
cos_gam = cos(gam_max);
cos_eff = cos(gam_eff);
cos_beam = cos(gam_beam);

view_tot = zeros(size(PHI));
beam_tot = zeros(size(PHI));

R_ij = Re*cat(3, sTH.*cPHI, sTH.*sPHI, cTH);
r_sat = reshape(pos', [1,1,3]);
R_sat = repmat(r_sat, size(PHI));
R_beam = R_ij - R_sat;
R_beam = R_beam./repmat(sqrt(sum(R_beam.^2, 3)), [1,1,3]);
        
        % create satellite view
val_view = dot(R_ij, R_sat, 3)/(Re*oe(1));
view = (val_view >= cos_gam);
view_tot = view_tot + view;

r_cen = R_beam(coords(1),coords(2),:);
beam = R_beam(:,:,1)*r_cen(1) + R_beam(:,:,2)*r_cen(2) + R_beam(:,:,3)*r_cen(3) > cos(th_beam);
beam = beam & view;
beam_tot = beam_tot + beam;

% earth_image = imread('earth.jpg');
figure
I = imshow(earth);
hold on;
coord = squeeze(coords);
plot(coords(2), coords(1), 'r.', 'MarkerSize', 5);
mask = (0.2).^(0.1*view_tot + beam_tot);
set(I, 'alphadata', mask);
