function [dq1,dq2,dinc,dOmega,dlambda] = deltaOE(C1,C2,alpha,SMA,inc)
% deltaOE returns deltas of OE between central(target) point and a
% satellite, according to its position in orbital frame.

dq1 = -C1 .* sin(alpha)/2/SMA;
dq2 = -C1 .* cos(alpha)/2/SMA;
dinc = C2 .* cos(alpha)/SMA;
dOmega = -C2 .* sin(alpha)/SMA/sin(inc);
dlambda = -dOmega * cos(inc);

end

