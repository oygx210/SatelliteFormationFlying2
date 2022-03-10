function M = mean_anomaly(e,nu)

sinE = sqrt(1 - e.^2).*sin(nu)./(1 + e.*cos(nu));
cosE = (e + cos(nu))./(1+e.*cos(nu));

E = atan2(sinE,cosE);


M = E - e.*sin(E);
if M < 0 
    M = 2*pi + M;
end
end