function rv_prime_hcw = hcw(t, rv_orb, mean_motion)

% Right hand side of the Hill-Clohessy-Wiltshire equations of relative
% motion

rv_prime_hcw = [rv_orb(4); rv_orb(5); rv_orb(6);
                -2 * mean_motion * rv_orb(6);
                -mean_motion^2 * rv_orb(2);
                2 * mean_motion * rv_orb(4) + 3 * mean_motion^2 * rv_orb(3)];

end