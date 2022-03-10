
consts = startup_formation_control();

[X,Y,Z] = sphere(100);
r1 = consts.rEarth + 600e3;
r2 = consts.rEarth + 800e3;

plot_orbits([rv_cluster1; rv_cluster2], consts);
hold on;
surf(X*r1, Y*r1, Z*r1, 'FaceAlpha',0.1, 'EdgeColor', 'none', 'FaceColor', '#0072BD');
hold on;
surf(X*r2, Y*r2, Z*r2, 'FaceAlpha',0.1, 'EdgeColor', 'none', 'FaceColor', '#0072BD');
hold on;
r_obs = spacecraft.dr_observation;
surf(X*r_obs + rv_cluster1(1), Y*r_obs + rv_cluster1(2), Z*r_obs + rv_cluster1(3), 'FaceAlpha',0.1, 'EdgeColor', 'none', 'FaceColor', '#A2142F');
hold on;
surf(X*r_obs + rv_cluster2(1), Y*r_obs + rv_cluster2(2), Z*r_obs + rv_cluster2(3), 'FaceAlpha',0.1, 'EdgeColor', 'none', 'FaceColor', '#A2142F');
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
axis off;