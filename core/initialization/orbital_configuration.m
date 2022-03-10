function formation_configuration = orbital_configuration(letters, relative_orbit_type, min_ISD , initial_phase)

switch letters
    case "Sk"
        formation_geometry = [
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
        -3.5 -0.5;
                
    %   "k"
        1+1  3.5;
        1+1 2.5;
        1+1 1.5;
        1+1 0.5
        1+1 -0.5;
        1+1 -1.5;
        2+1 1;
        3+1 2;
        4+1 3;
        2+1 0;
        3+1 -1;
        4+1 -2;
        
           ];
        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
                
        formation_geometry(:,2) = -1*formation_geometry(:,2);
    case "AA"
        % building the letters
        formation_geometry = [
            2.12 2.06
            1.96 1.67
            1.87 1.19
            1.75 0.77
            1.52 0.35
            1.35 0.77
            1.16 1.19
            0.99 1.67
            1.32 1.67
            1.65 1.67
            0.84 2.07
            0.64 2.55
            2.27 2.55
            
            2.12 2.06
            1.96 1.67
            1.87 1.19
            1.75 0.77
            1.52 0.35
            1.35 0.77
            1.16 1.19
            0.99 1.67
            1.32 1.67
            1.65 1.67
            0.84 2.07
            0.64 2.55
            2.27 2.55
            ];

        formation_geometry(14:end,1) = formation_geometry(14:end,1)  + 2.27;
        
        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
                
        formation_geometry(:,2) = -1*formation_geometry(:,2);
    case "MIPT"
        formation_geometry = [
            189.7 501.1;
            227.6 532.6;
            265.6 260.6;
            273.6 559.6;
            281.0 229.0;
            290.6 283.6;
            315.0 577.0;
            334.1 275.7;
            361.6 586.6;
            378.5 268.6;
            413.7 585.1;
            420.1 369.1;
            421.7 261.7;
            431.7 312.1;
            438.6 402.7;
            453.7 573.7;
            462.1 249.6;
            473.6 395.6;
            494.1 554.7;
            506.6 241.6;
            515.6 386.6;
            527.6 530.7;
            550.7 231.6;
            555.0 504.0;
            565.6 379.6;
            574.6 476.6;
            590.1 443.7;
            602.0 224.0;
            604.6 406.6;
%             617.6 366.6;
%             652.7 211.0;

           ];
        
        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
                
        formation_geometry(:,2) = -1*formation_geometry(:,2);
        
    case "SK"
        % building the letters
        formation_geometry = [
            2.12 2.06
            1.96 1.67
            1.87 1.19
            1.75 0.77
            1.52 0.35
            1.35 0.77
            1.16 1.19
            0.99 1.67
            1.32 1.67
            1.65 1.67
            0.84 2.07
            0.64 2.55
            2.27 2.55
            
            2.12 2.06
            1.96 1.67
            1.87 1.19
            1.75 0.77
            1.52 0.35
            1.35 0.77
            1.16 1.19
            0.99 1.67
            1.32 1.67
            1.65 1.67
            0.84 2.07
            0.64 2.55
            2.27 2.55
            ];

        formation_geometry(14:end,1) = formation_geometry(14:end,1)  + 2.27;
        
        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
                
        formation_geometry(:,2) = formation_geometry(:,2);

    case "IAA_logo"
        formation_geometry = [
            9.78 1.36
            9.7 1.88
            10.34 2.04
            10.49 1.45
            9.46 1.04
            9.12 0.87
            8.91 1.18];
        
        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
        formation_geometry(:,1) = -1*formation_geometry(:,1);
                
        phi = 0:2*pi/19:(2*pi - 2*pi/19);
        [x,y] = pol2cart(phi, 1.5*ones(1,length(phi)));
        formation_geometry(9:27,1) = x;
        formation_geometry(9:27,2) = y;
                   
        formation_geometry(:,2) = -1*formation_geometry(:,2);
        
    case "2_sats_test"
        
        formation_geometry = [
            0 0
            min_ISD 0];
    case "AAS_letters"
                
        formation_geometry = [

    %   "A"
        0 0;
        1 2
        2 4;
        3 6;
        4 8
        5 10;
        6 8;
        7 6;
        8 4;
        9 2;
        10 0;
        4 3;
        6 3;
        
    %   "A"
        0 0;
        1 2
        2 4;
        3 6;
        4 8
        5 10;
        6 8;
        7 6;
        8 4;
        9 2;
        10 0;
        4 3;
        6 3;
        
    %   "S"
        0 2;
        1 0;
        3 0;
        5 0;
        6 1;
        6 3;
        5 5;
        3 5;
        1 5;
        0 6;
        0 8;
        1 10;
        3 10;
        5 10;
        6 8;
        ];
     
        formation_geometry(14:26,1) = formation_geometry(14:26,1)  + 13;
        formation_geometry(27:end,1) = formation_geometry(27:end,1)  + 27;

        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
    
    case "IAA_letters"
        
        formation_geometry = [
    %   "I"
        0 0;
        1 0;
        2 0;
        3 0;
        4 0;
        0 10;
        1 10;
        2 10;
        3 10;
        4 10;
        2 1;
        2 3;
        2 5;
        2 7 
        2 9;
       
    %   "A"
        0 0;
        1 2
        2 4;
        3 6;
        4 8
        5 10;
        6 8;
        7 6;
        8 4;
        9 2;
        10 0;
        4 3;
        6 3;
        
    %   "A"
        0 0;
        1 2
        2 4;
        3 6;
        4 8
        5 10;
        6 8;
        7 6;
        8 4;
        9 2;
        10 0;
        4 3;
        6 3;
        ];
     
        formation_geometry(16:28,1) = formation_geometry(16:28,1)  + 10;
        formation_geometry(29:end,1) = formation_geometry(29:end,1)  + 25;

        geometrical_center = geometric_median(formation_geometry');
        geometrical_center = geometrical_center';
        
        formation_geometry = [geometrical_center; formation_geometry];
        formation_geometry = formation_geometry - geometrical_center;
        
     otherwise 
        disp('Contact our sales manager to order a new logo for demonstration');
end

ISD_raw = 99999;
for i = 1:size(formation_geometry,1)-1
    for j = 1:size(formation_geometry,1)-1
        if i ~= j
            ISD_raw_new = sqrt((formation_geometry(i+1,1) - formation_geometry(j+1,1))^2 + (formation_geometry(i+1,2) - formation_geometry(j+1,2))^2);
            if ISD_raw_new < ISD_raw
                ISD_raw = ISD_raw_new;
            end
        end
    end
end
    
scale = min_ISD/ISD_raw;
formation_geometry = formation_geometry*scale;

[phi, rho] = cart2pol(formation_geometry(:,1), formation_geometry(:,2));

switch relative_orbit_type
    case "GCO"
        c1 = rho;
        c2 = sqrt(3)/2*rho;
    case "PCO"
        c1 = rho;
        c2 = rho;
end

c3 = zeros(length(c1),1);
alpha = phi + (phi<0) * 2*pi + initial_phase;


formation_configuration = [c1, c2, c3, alpha]';
formation_configuration(:,1) = zeros(4,1);

% figure;
% plot(formation_geometry(2:end,1)/1000, formation_geometry(2:end,2)/1000, 'sk', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 12);
% axis square;
% hold on;
% plot(formation_geometry(1,1)/1000, formation_geometry(1,2)/1000, 'or', 'MarkerFaceColor',  'r', 'MarkerSize', 6);
% legend('Satellites', 'Geometrical median');
% xlabel('x, km');
% ylabel('y, km');
% image_size = 2*max(vecnorm(formation_geometry(2:end,:)'));

end