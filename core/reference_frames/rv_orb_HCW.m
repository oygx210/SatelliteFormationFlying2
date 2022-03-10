function rv_orb = rv_orb_HCW(formation_geometry, n, t)

c1 = formation_geometry(1,:);
c2 = formation_geometry(2,:);
c3 = formation_geometry(3,:);
alpha = formation_geometry(4,:);

x = c1.*cos(mod(n*t+alpha, 2*pi)) + c3;
y = c2.*sin(mod(n*t+alpha, 2*pi));
z = c1/2.*sin(mod(n*t+alpha, 2*pi));

vx = -n*c1.*sin(mod(n*t+alpha, 2*pi));
vy = n*c2.*cos(mod(n*t+alpha, 2*pi));
vz = n*c1/2.*cos(mod(n*t+alpha, 2*pi));

rv_orb = [x; y; z; vx; vy; vz];

end


