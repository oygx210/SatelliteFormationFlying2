% generating control input

clear all;
consts = startup_formation_control();

r_ref = 7000e3;
inc = deg2rad(35);
n = sqrt(consts.muEarth/r_ref^3);
s = 3*consts.J2*consts.rEarth_equatorial^2/2/r_ref^2*(1+3*cos(2*inc));
c = sqrt(1 + s);

C = [0 0 0
     0 -(3*c^2-2)*n^2 0 
     0 0 (5*c^2-2)*n^2];
D = [0 0 -2*n*c
     0 0 0
     2*n*c 0 0];
 
A = [[zeros(3), eye(3)]; [C,D]];
  
B = [zeros(3); eye(3)];

Q = eye(6);
Rdiag = [1e-13; 1e-14; 1e-14];
R = diag(Rdiag);

[K] = lqr(A,B,Q,R);

e = [1 ; 0; 0; 0; 0; 0]; 
Kpaper = [14.1421 -0.0049 0 6.1874 0 0
          -0.0049 14.1421 0 0 6.1874 0 
          0 0 12.2474 0 0 5.8732];

u = -Kpaper*e;

% 
% options_precision = odeset('RelTol',1e-12,'AbsTol',1e-12);
% tic;
% [t_out, rv_orb_out] = ode45(@(t, rv_orb) rhs(t, rv_orb, A, B, K, spacecraft), [0 consts.day2sec], rv_orb, options_precision);
% toc;
% rv_orb_out = rv_orb_out';
% 
% e = vecnorm(rv_orb_out(7:9,:) - rv_orb_out(1:3,:));
% % figure;
% % plot(t_out/60, e);
% % xlabel('time, min');
% % ylabel('position erro, min');
% % 
% % 
% % figure;
% % plot(t_out/60, rv_orb_out(7,:)- rv_orb_out(1,:));
% % hold on;
% % plot(t_out/60, rv_orb_out(8,:)- rv_orb_out(2,:));
% % hold on;
% % plot(t_out/60, rv_orb_out(9,:)- rv_orb_out(3,:));
% % legend('x, m', 'y, m', 'z, m');
% % xlabel('time, min');
% % ylabel('position erro, min');
% 
% figure;
% plot3(rv_orb_out(7,:), rv_orb_out(8,:), rv_orb_out(9,:));
% xlabel('x, m');
% ylabel('y, m');
% zlabel('z, m');
% 
% 
% function rv_orb_prime = rhs(t, rv_orb, spacecraft)
% 
% global A B K
% e = rv_orb(7:12) - rv_orb(1:6);
% u = K*e;
% 
% if norm(u) > spacecraft.unit_thrust
%     u = u./norm(u)*spacecraft.unit_thrust;
% end
% 
% rv_orb_prime(1:6,1) = A*rv_orb(1:6);
% rv_orb_prime(7:12,1) = A*rv_orb(7:12) - B*u;
% 
% end
% 
% 
