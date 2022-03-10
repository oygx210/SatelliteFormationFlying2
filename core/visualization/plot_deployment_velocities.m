function plot_deployment_velocities(rv_orb_deployment, deployer)

dV_deployment_mag = sqrt(dot(rv_orb_deployment(4:6,:),rv_orb_deployment(4:6,:)));
dV_deployment_theta = atand(sqrt(rv_orb_deployment(5,:).^2 + rv_orb_deployment(6,:).^2) ./ rv_orb_deployment(4,:));
dV_deployment_phi = atan2(rv_orb_deployment(6,:),rv_orb_deployment(5,:));

figure();
subplot(1,3,1);
histogram(dV_deployment_mag);
xlabel('deployment velocity, [m/s]', 'FontSize', 14);
ylabel('number of samples', 'FontSize', 14);
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
        '3\sigma (dV_{||}) = ', num2str(deployer.dV_sigma*3, 2), ' [m/s]']);
set(gca, 'FontSize', 14);

subplot(1,3,2);
histogram(dV_deployment_theta);
xlabel('deployment direction error, degrees', 'FontSize', 14);
ylabel('number of samples', 'FontSize', 14);
legend(['dV = ', num2str(deployer.dV_dir(1), 2), ' [m/s]; ',...
       '3\sigma (dV_{\perp}) = ', num2str(deployer.dV_sigma_wide*3, 2), ' [m/s]']);
set(gca, 'FontSize', 14);

title('Target satellite deployment conditions', 'FontSize', 14);
subplot(1,3,3);
polarplot(dV_deployment_phi, dV_deployment_theta, '+m');
ax = gca;
rruler = ax.RAxis;
set(gca, 'FontSize', 14);


end