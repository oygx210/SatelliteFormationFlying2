function y = floor_rad(x)
%Returns value less than +-pi

x = mod(x,2*pi);

if x > pi
    y = 2*pi - x;
else
    y = x;
end

end

