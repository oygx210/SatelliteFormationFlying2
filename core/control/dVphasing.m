function dVphasing = dVphasing(consts, h, l, N)

% The function returns dV required to phase a daughter satellite in a way
% to ensure required triangulation base l within N revolutions of daughter
% satellite

% The mother satellite is in circular LEO at altitude h

Tmother = 2*pi*sqrt((consts.rEarth + h)^3 / consts.muEarth);

alpha = 2 * asin(l/2/(consts.rEarth + h));
% add LOS check

Tdaughter = Tmother * (1 - alpha/2/pi/N);

a = (Tdaughter/2/pi)^(2/3) * consts.muEarth^(1/3);

dVphasing = sqrt(consts.muEarth*(2/(consts.rEarth + h) - 1/a)) - sqrt(consts.muEarth/(consts.rEarth + h));

end