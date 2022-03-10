function [t, rv_chaser] = singlesat_orbit_maintenance(rv_target_ECI, rv_chaser_ECI, schedule, dV1, dV2, consts, spacecraft)

% singlesat_orbit_maintenance adjust the orbit of a satellite to the required one
% It is know that the longest time for the Bi-Impuslive manuever does not
% exceed two orbit periods

dt = 1;
time_phase1 = ceil(schedule(1));
time_phase2 = ceil(schedule(2));
time_stop = ceil(schedule(3));

options = odeset('RelTol',1e-12,'AbsTol',1e-12); 

% Phase 0, propagation before the first impulse

rv_init = [rv_target_ECI ; rv_chaser_ECI];
[dt_vec_phase0, rv_ECI_phase0] = ode45(@(t, rv)J2_atmo_2sats(t, rv, consts, spacecraft), 0:dt:time_phase1, rv_init, options);
rv_ECI_phase0 = rv_ECI_phase0';

% First impulse

chaser_after_burn1 = ECI2orb(rv_ECI_phase0(1:6,end), rv_ECI_phase0(7:12,end), consts) + [0; 0; 0; dV1];
rv_ECI_phase0(7:12,end) = orb2ECI(rv_ECI_phase0(1:6,end), chaser_after_burn1, consts);

% Phase 1, propagation after the first impulse

rv_init = rv_ECI_phase0(:,end);
[dt_vec_phase1, rv_ECI_phase1] = ode45(@(t, rv)J2_atmo_2sats(t, rv, consts, spacecraft), time_phase1:dt:time_phase2, rv_init, options);
rv_ECI_phase1 = rv_ECI_phase1';

% Second impulse

chaser_after_burn2 = ECI2orb(rv_ECI_phase1(1:6,end), rv_ECI_phase1(7:12,end), consts) + [0; 0; 0; dV2];
rv_ECI_phase1(7:12,end) = orb2ECI(rv_ECI_phase1(1:6,end), chaser_after_burn2, consts);

% Propagation after the second impulse
rv_init = rv_ECI_phase1(:,end);
[dt_vec_phase2, rv_ECI_phase2] = ode45(@(t, rv)J2_atmo_2sats(t, rv, consts, spacecraft), time_phase2:dt:time_stop, rv_init, options);
rv_ECI_phase2 = rv_ECI_phase2';

t = [dt_vec_phase0; dt_vec_phase1; dt_vec_phase2];
rv_chaser = [rv_ECI_phase0(7:12,:), rv_ECI_phase1(7:12,:), rv_ECI_phase2(7:12,:)];

end