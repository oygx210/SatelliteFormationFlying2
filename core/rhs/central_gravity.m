function rv_prime = central_gravity(t, rv, consts)

% Right hand side for the equations of motion for a body in central gravity
% field

rv_prime = [rv(4); rv(5); rv(6);...
            -consts.muEarth*rv(1)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2);...
            -consts.muEarth*rv(2)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2);...
            -consts.muEarth*rv(3)/(rv(1)^2 + rv(2)^2 + rv(3)^2)^(3/2);
            ];
end