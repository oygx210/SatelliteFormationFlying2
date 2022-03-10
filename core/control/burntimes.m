function times = burntimes(alpha,period)
%burntimes returns times when both impulses are delivered
    
    
    %phase 0 before first burn
    
    theta(1) = 2*pi - alpha;% time of first burn
    
    %phase 1 tranfer to formation
    
    theta(2) = theta(1) + pi;% time of second burn
    
    %phase 2, demonstration of formation
    
    theta(3) = 4*pi;% end of first period of demonstration
    
    
    times = theta/2/pi .* period;
end

