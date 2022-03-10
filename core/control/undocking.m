function [t,rv] = undocking(rv_ECI, deployer, formation, spacecraft, consts)

options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);

t = 0 ;

for i = 1:formation.N_sats 
    rv(i*6-5:i*6,1) = rv_ECI;
end

% Generating initiali conditions for deployment
rv_orb_deployment = [zeros(1,formation.N_sats); zeros(1,formation.N_sats); zeros(1,formation.N_sats); 
                     deployer.dep2orb * [deployer.dV_dir(1)*ones(1,formation.N_sats) + normrnd(deployer.dV_mean*ones(1,formation.N_sats), deployer.dV_sigma*ones(1,formation.N_sats));
                     normrnd(zeros(1,formation.N_sats), deployer.dV_sigma_wide*ones(1,formation.N_sats));
                     normrnd(zeros(1,formation.N_sats), deployer.dV_sigma_wide*ones(1,formation.N_sats))]];                        


for i = 1:formation.N_sats

[t_vec_deployment, rv_ECI_deployment] = ode45(@(t, rv) rhs_Formation_inertial(t, rv, consts, spacecraft),t(end): t(end) + formation.deployment_interval, rv(:,end), options_precision);

rv_ECI_deployment = rv_ECI_deployment';
t = [t; t_vec_deployment];
rv = [rv, rv_ECI_deployment];

rv(i*6-5:i*6,end) = orb2ECI(rv(i*6-5:i*6,end), rv_orb_deployment(:,i), consts);

end


end