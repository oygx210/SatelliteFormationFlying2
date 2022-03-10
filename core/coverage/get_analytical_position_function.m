function r = get_analytical_position_function(oe, consts)

    syms t
    
    C = 3 / 2 * (sqrt(consts.muEarth)*consts.J2*consts.rEarth_equatorial^2);

    RAAN_dot =  - C / ((1 - oe(2)^2)^2*oe(1)^(7/2)) * cos(oe(3));

    AOP_dot = - C / ((1 - oe(2)^2)^2 * oe(1)^(7/2)) * (5/2 * sin(oe(3))^2 - 2);

    n = sqrt(consts.muEarth/oe(1)^3);
    arg1 = oe(5) + AOP_dot*t + oe(6) + n*t;

    M1 = [cos(arg1), sin(arg1), 0
        -sin(arg1), cos(arg1), 0
        0, 0, 1];
    M2 = [1, 0, 0
        0, cos(oe(3)), sin(oe(3))
        0, -sin(oe(3)), cos(oe(3))];

    arg2 = oe(4) + RAAN_dot*t;

    M3 = [cos(arg2), sin(arg2),0
        -sin(arg2), cos(arg2), 0
        0, 0, 1];

    Rotation_matrix = M1 * M2 * M3;
    Rotation_matrix = Rotation_matrix';
    
    r = Rotation_matrix * [oe(1); 0; 0];   
    
