function [r,v,C1,C2,alpha] = positions(dist,text,n)

% Returns coordinates of satellites in formation in orbital frame

% Letters

    switch text
        case "ABC"
        letters = [

    %   "A"
        -13 -4;
        -12 -2;
        -11 0;
        -10 2;
        -10 -2;
        -9 4;
        -8 2;
        -8 -2;
        -7 0;
        -6 -2;
        -5 -4;

    %   "B"
        -1 4;
        -1 2;
        -1 0;
        -1 -2;
        -1 -4;
        1 4;
        1 0;
        1 -4;
        2 3;
        2 1;
        3 -1;
        3 -3;

    %   "C"
        7 2;
        7 0;
        7 -2;
        9 4;
        9 -4;
        11 4;
        11 -4;
        13 2;
        13 -2

        ];

        case "DEF"
        letters = [

    %   "D"
        -12 -4;
        -12 -2;
        -12 0;
        -12 2;
        -12 4
        -10 4;
        -10 -4;
        -8 -4;
        -8 4;
        -6 2;
        -6 -2;
        -6 0;

    %   "E"
        -2 4;
        -2 2;
        -2 0;
        -2 -2;
        -2 -4;
        0 4;
        0 0;
        0 -4;
        2 4;
        2 -4;

    %   "F"
        7 4;
        7 2;
        7 0;
        7 -2;
        7 -4;
        9 4;
        9 0;
        11 4;
        11 0;
        13 4;

        ];
    
        case "IAA"
        letters = [

    %   "I"
        -13+2 5;
        -11+2 5;
        -9+2 5;
        -11+2 2.5;
        -11+2 0;
        -11+1 -2.5;
        -13+1 -5;
        -11+1 -5;
        -9+1 -5;

    %   "A"
        -6 -5;
        -5 -2.5;
        -4 0;
        -3 2.5;
        -2 5;
        -1 2.5;
        0 0;
        -3 -2.5;
        -1 -2.5;
        1 -2.5;
        2 -5;
        
        
    %   "A"
        -6+11 -5;
        -5+11 -2.5;
        -4+11 0;
        -3+11 2.5;
        -2+11 5;
        -1+11 2.5;
        0+11 0;
        -3+11 -2.5;
        -1+11 -2.5;
        1+11 -2.5;
        2+11 -5
        
        ];
    
        case "AA"
        letters = [

    %   "A"
        0 0;
        -6 -5;
        -5 -2.5;
        -4 0;
        -3 2.5;
        -2 5;
        -1 2.5;
        -3 -2.5;
        -1 -2.5;
        1 -2.5;
        2 -5;
                
    %   "A"
        -6+11 -5;
        -5+11 -2.5;
        -4+11 0;
        -3+11 2.5;
        -2+11 5;
        -1+11 2.5;
        0+11 0;
        -3+11 -2.5;
        -1+11 -2.5;
        1+11 -2.5;
        2+11 -5
        
        ];
        
        case "Sk"
        letters = [

    %   "S"
        0 0;
        0 4.5;
        -0.5 5;
        -1.5 5;
        -2.5 5;
        -3.5 4.5;
        -3.5 3.5;
        -3.5 2.5;
        -2.5 1.5;
        -1.5 1.5;
        -0.5 1.5;
        0 1;
        0 -1;
        -0.5 -1.5;
        -1.5 -1.5;
        -2.5 -1.5;
        -2.5 -1.5;
        -3.5 -0.5;
                
    %   "k"
        1+1  3.5;
        1+1 2.5;
        1+1 1.5;
        1+1 0.5
        1+1 -0.5;
        1+1 -1.5;
        1+1 0.5;
        2+1 1;
        3+1 2;
        4+1 3;
        2+1 0;
        3+1 -1;
        4+1 -2;
        
        ];
    
        otherwise
            warning('Unexpected text, use ABC, DEF, IAA, AA or Sk')
    end  
    
%     Tagir's version

    [ABC_phi, ABC_rho] = cart2pol(letters(:,1),letters(:,2));
%     ABC_phi = ABC_phi + pi;
    ABC_rho = dist .* ABC_rho;
    disp([ABC_phi, ABC_rho])
    A = ABC_rho/2;
    c1 = A .* sin(ABC_phi);
    c2 = A .* cos(ABC_phi);
    
    C1 = c1;
    C2 = c2;
    alpha = ABC_phi;
    
    x = [];
    y = [];
    z = [];

    x = 2*c1;
    y = sqrt(3)*c2;
    z = c2;
    r = [x'; y'; z'];
    
    vx = 2*n*c2;
    vy = -sqrt(3).*n.*c1;
    vz = n.*c1;
    v = [vx'; vy'; vz'];
    
    disp(size(x));
    scatter3(x,y,z);
    grid on;
    daspect([1 1 1]);

%     Arcticle version

%     [alpha, rho] = cart2pol(letters(:,1),letters(:,2));
%     ang_offset = 0; %rotate text
%     alpha = alpha + ang_offset;
%     rho = dist * rho;
%     C1 = rho;
%     C2 = sqrt(3)/2 * rho;  %(Use C2 = 1/2 * rho; for the PCO)
%     
%     x = C1 .* cos(alpha);
%     y = C2 .* sin(alpha);
%     z = C1/2 .* sin(alpha);
% 
%     vx = -n * C1 .* sin(alpha);
%     vy = n * C2 .* cos(alpha);
%     vz = n * C1/2 .* cos(alpha);
%     
%     r = [x y z];
%     v = [vx vy vz];
end

