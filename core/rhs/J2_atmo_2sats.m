function rv_prime = J2_atmo_2sats(t, Targer_Chaser_rv, consts, spacecraft)

%% The function simulates two satellite orbital dynamics. 
% The satellites are the same meaning its physical properties. This might be improved at the further studies.

% Right-hand side for the equation of motion of a satellite taking into
% account:
% 1. Oblated Earth
% 2. Atmospheric drag (piecewise exponential atmospheric density function)

R = consts.rEarth_equatorial;
J2 = consts.J2;
mu = consts.muEarth;
delta = 3/2*J2*mu*R^2;

Target.rv = Targer_Chaser_rv(1:6);
Chaser.rv = Targer_Chaser_rv(7:12);

target_altitude = vecnorm(Target.rv(1:3));
chaser_altitude = vecnorm(Chaser.rv(1:3));
target_altitude_3 = target_altitude^3;
chaser_altitude_3 = chaser_altitude^3;

% J2 effect
Target.acceleration_J2 = delta*Target.rv(1:3)/target_altitude^5*(5*Target.rv(3)^2/target_altitude^2 - 1)- 2*delta/target_altitude^5*[0; 0; Target.rv(3)];
Chaser.acceleration_J2 = delta*Chaser.rv(1:3)/chaser_altitude^5*(5*Chaser.rv(3)^2/chaser_altitude^2 - 1)- 2*delta/chaser_altitude^5*[0; 0; Chaser.rv(3)];

% Atmospheric drag
Target.vRelativeECI = Target.rv(4:6) - cross(consts.wEarth, Target.rv(1:3));
Target.rhoAtmo = CIRA72(consts, (vecnorm(Target.rv(1:3)) - consts.rEarth));
Target.acceleration_AD = - 0.5 * spacecraft.Cdrag * spacecraft.DragArea / spacecraft.mass * Target.rhoAtmo * Target.vRelativeECI * vecnorm(Target.vRelativeECI); 

Chaser.vRelativeECI = Chaser.rv(4:6) - cross(consts.wEarth, Chaser.rv(1:3));
Chaser.rhoAtmo = CIRA72(consts, (vecnorm(Chaser.rv(1:3)) - consts.rEarth));
Chaser.acceleration_AD = - 0.5 * spacecraft.Cdrag * spacecraft.DragArea / spacecraft.mass * Chaser.rhoAtmo * Chaser.vRelativeECI * vecnorm(Chaser.vRelativeECI); 

rv_prime = [Target.rv(4:6);...
            -consts.muEarth*Target.rv(1:3)/target_altitude_3 + Target.acceleration_J2 + Target.acceleration_AD;
            Chaser.rv(4:6);...
            -consts.muEarth*Chaser.rv(1:3)/chaser_altitude_3 + Chaser.acceleration_J2 + Chaser.acceleration_AD
            ];
end