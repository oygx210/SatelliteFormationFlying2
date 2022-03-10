r_obs = 20e3;
d_debris = 10e-2;
alpha = d_debris / (2*r_obs);
lambda = 500e-9;
D_aperture1 = 1.22 * lambda / alpha;
%
r_obs = 600e3;
res = 2.8;
D_aperture2 = 2.44 * r_obs * lambda / res;

r_obs = 500e3;
res = 9.6;
D_aperture3 = 2.44 * r_obs * lambda / res;


habble = 0.1/360 /180*pi;
