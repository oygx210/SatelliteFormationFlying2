function output = singlesat_integrator(alpha,orbit_period,dV1,dV2,target_rv_ECI,consts)
% singlesat_integrator propagates satelite with all required burns

mean_motion = 2*pi/orbit_period;

%% Mission control

dt = 20; % s

time_start = -200;
times = burntimes(alpha,orbit_period); % s
time_phase1 = times(1);
time_phase2 = times(2);
time_stop = times(3);

tspan_phase0 = [time_start : dt : time_phase1];
tspan_phase1 = [time_phase1 : dt : time_phase2];
tspan_phase2 = [time_phase2 : dt : time_stop];

dt_vec = zeros(1,length(tspan_phase0) + length(tspan_phase1) + length(tspan_phase2));

options = odeset('RelTol',1e-12,'AbsTol',1e-12); 

chaser.rv_ECI = zeros(6, length(dt_vec));


%% Phase 0, propagation before the first impulse

chaser_init = target_rv_ECI(:,1);
[dt_vec_phase0, chaser.rv_ECI_phase0] = ode45(@(t, rv)central_gravity(t, rv, consts), tspan_phase0, chaser_init, options);
chaser.rv_ECI_phase0 = chaser.rv_ECI_phase0';

chaser.rv_ECI(:,(1:length(tspan_phase0))) = chaser.rv_ECI_phase0;
dt_vec(1, 1:length(dt_vec_phase0)) = dt_vec_phase0;


%% First impulse

chaser_init_orb = ECI2orb(target_rv_ECI(:,length(tspan_phase0)), chaser.rv_ECI_phase0(:,end), consts) + [0; 0; 0; dV1];
chaser_init = orb2ECI(target_rv_ECI(:,length(tspan_phase0)), chaser_init_orb, consts);


%% Phase 1, propagation after the first impulse


[dt_vec_phase1, chaser.rv_ECI_phase1] = ode45(@(t, rv)central_gravity(t, rv, consts), tspan_phase1, chaser_init, options);
chaser.rv_ECI_phase1 = chaser.rv_ECI_phase1';

chaser.rv_ECI(:,(length(tspan_phase0)) : length(dt_vec_phase0) + length(dt_vec_phase1)-1) = chaser.rv_ECI_phase1;
dt_vec(1,(length(dt_vec_phase0)+1) : length(dt_vec_phase0) + length(dt_vec_phase1)) = dt_vec_phase1;


%% Second impulse

chaser_init_orb = ECI2orb(target_rv_ECI(:,length(tspan_phase0) + length(dt_vec_phase1)), chaser.rv_ECI_phase1(:,end), consts) + [0; 0; 0; dV2];
chaser_init = orb2ECI(target_rv_ECI(:,length(tspan_phase0) + length(dt_vec_phase1)), chaser_init_orb, consts);

%% Propagation after the second impulse

[dt_vec_phase2, chaser.rv_ECI_phase2] = ode45(@(t, rv)central_gravity(t, rv, consts), tspan_phase2, chaser_init, options);
chaser.rv_ECI_phase2 = chaser.rv_ECI_phase2';

chaser.rv_ECI(:,(length(tspan_phase0)+length(tspan_phase1)-1) : end-2) = chaser.rv_ECI_phase2;
dt_vec(1,(length(tspan_phase0)+length(dt_vec_phase1)-1) : end-2) = dt_vec_phase2;


% plot(chaser.rv_ECI(1,1:564)-target_rv_ECI(1,1:564))
% hold on
% plot(chaser.rv_ECI(2,1:564)-target_rv_ECI(2,1:564))
% hold on
% plot(chaser.rv_ECI(3,1:564)-target_rv_ECI(3,1:564))


for i = 1:length(dt_vec)
    
    chaser.rv_orb(:,i) = ECI2orb(target_rv_ECI(:,i), chaser.rv_ECI(:,i), consts);
    
end

output = chaser.rv_orb;

end
