function [t, rv] = formation_maintenance_after_undocking(rv_ECI, consts, spacecraft, formation)

% here the Bi-Impulsive maneuvers are to be done in order to adjust
% reference orbit for satellites in formation

target_oe = vecRV2OE(rv_ECI(1:6), consts);

C1 = formation.geometry(1,:);
C2 = formation.geometry(2,:);
alpha = formation.geomety(4,:);

[dq1,dq2,dinc,dOmega,dlambda] = deltaOE(C1,C2,alpha,target_oe(1),target_oe(3));

%% Impulses
p=1; %basically, ratio between first and second out-of-plane burns. See Vaddi et al., page 264, section called "Aside"
gamma = sqrt(target_oe(1)/consts.muEarth);  %[s/m]

dV1 = [];
dV2 = [];

for sat=1:formation.N_sats
    %deltaVs are in m/s 
    dV_h1 = p/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(target_oe(3)).^2);
    dV_h2 = (1-p)/gamma * sqrt(dinc.^2 + dOmega.^2 .* sin(target_oe(3)).^2);
    dV_r1 = -sqrt(dq1.^2 + dq2.^2)/2/gamma;
    dV_r2 = sqrt(dq1.^2 + dq2.^2)/2/gamma;

    dV1(:, sat) = [0 ; dV_h1(sat) ; dV_r1(sat)];
    dV2(:, sat) = [0 ; dV_h2(sat) ; dV_r2(sat)];
end

schedule = schedule_burns(rv_ECI(:,1), consts, formation);

for i = 1:formation.N_sats
    [t, rv(i*6-5:i*6,:)] = singlesat_orbit_maintenance(rv_ECI(:,1), rv_ECI(:,i), schedule(:,i), dV(1:3,i), dV(4:6,i), consts, spacecraft);
end

end