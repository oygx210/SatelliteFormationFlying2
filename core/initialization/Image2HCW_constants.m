function HCW_constants = Image2HCW_constants(graphics, relative_orbit_type, IPD_min, alpha_image)

% Step 1. Choosing graphics to demonstrate
    switch graphics
        case "olympic_rings"
            pixels_coordinates = load('olympic_rings.txt');
        case "eiffel_tower"
            pixels_coordinates = load('eiffel_tower.txt');
        case "XY"
            pixels_coordinates = load('XY_letters.txt');
        case "rho_phi"
            pixels_coordinates = load('rho_phi_letters.txt');
        case "AA"            
            % building the letters
            pixels_coordinates = [
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
            pixels_coordinates(14:end,1) = pixels_coordinates(14:end,1)  + 2.27;
            A = [[-1 0];
                [0 -1]];
            pixels_coordinates = pixels_coordinates';
            
            for i = 1:size(pixels_coordinates,2)
                pixels_coordinates(:,i) = A*pixels_coordinates(:,i);
            end
            pixels_coordinates = pixels_coordinates';
            
        case "IAA_logo" 
            
            pixels_coordinates = [
            9.78 1.36
            9.7 1.88
            10.34 2.04
            10.49 1.45
            9.46 1.04
            9.12 0.87
            8.91 1.18];
            
            geometrical_center = geometric_median(pixels_coordinates');

            geometrical_center = geometrical_center';
 
            pixels_coordinates = pixels_coordinates - geometrical_center;

            phi = 0:2*pi/19:(2*pi - 2*pi/19);
            [x,y] = pol2cart(phi, 1.5*ones(1,length(phi)));
            pixels_coordinates(8:26,1) = x;
            pixels_coordinates(8:26,2) = y;
            A = [[1 0];
                [0 -1]];
            pixels_coordinates = pixels_coordinates';
            
            for i = 1:size(pixels_coordinates,2)
                pixels_coordinates(:,i) = A*pixels_coordinates(:,i);
            end
            pixels_coordinates = pixels_coordinates';

    end
    figure('Name', 'Orbital configuration design', 'NumberTitle', 'off');
    subplot(2,2,1);
    plot(pixels_coordinates(:,1), pixels_coordinates(:,2), 'or', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 10);
    axis equal;
    grid on;
% Step 2. Building the image wrt to the pixels geometrical median
    geometrical_center = geometric_median(pixels_coordinates');
    geometrical_center = geometrical_center';

    pixels_coordinates = [geometrical_center; pixels_coordinates];
    pixels_coordinates = pixels_coordinates - geometrical_center;

    subplot(2,2,2);
    plot(pixels_coordinates(1,1), pixels_coordinates(1,2), 'ok', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 5);
    hold on;
    plot(pixels_coordinates(2:end,1), pixels_coordinates(2:end,2), 'or', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 10);
    axis equal;
    grid on;

% Step 3. Scaling image to satisfy requirements on the minimum IPD
    IPD_raw = 99999;
    for i = 1:size(pixels_coordinates,1)-1
        for j = 1:size(pixels_coordinates,1)-1
            if i ~= j
                ISD_raw_new = sqrt((pixels_coordinates(i+1,1) - pixels_coordinates(j+1,1))^2 + (pixels_coordinates(i+1,2) - pixels_coordinates(j+1,2))^2);
                if ISD_raw_new < IPD_raw
                    IPD_raw = ISD_raw_new;
                end
            end
        end
    end

    scale = IPD_min/IPD_raw;
    pixels_coordinates = pixels_coordinates*scale;

    subplot(2,2,3);
    plot(pixels_coordinates(1,1), pixels_coordinates(1,2), 'ok', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 5);
    hold on;
    plot(pixels_coordinates(2:end,1), pixels_coordinates(2:end,2), 'or', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 10);
    axis equal;
    grid on;
    xlabel("x-axis, meters");
    ylabel("y-axis, meters");
% Step 4. Mirroring image wrt to y axes to get the proper view during demonstration    
    
    A = [[-1 0];
        [0 1]];
    pixels_coordinates = pixels_coordinates';

    for i = 1:size(pixels_coordinates,2)
        pixels_coordinates(:,i) = A*pixels_coordinates(:,i);
    end
    
    pixels_coordinates = pixels_coordinates';

    subplot(2,2,4);
    plot(pixels_coordinates(1,1), pixels_coordinates(1,2), 'ok', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 5);
    hold on;
    plot(pixels_coordinates(2:end,1), pixels_coordinates(2:end,2), 'or', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 10);
    axis equal;
    grid on;

 % Step 5. Defining HCW constants for an orbital configuration
    [alpha0, rho] = cart2pol(pixels_coordinates(:,1), pixels_coordinates(:,2));
    alpha = alpha0 + alpha_image;
    
    switch relative_orbit_type
        case "GCO"
            c1 = rho;
            c2 = sqrt(3)/2*rho;
        case "PCO"
            c1 = rho;
            c2 = rho;
    end
    
    c3 = zeros(length(c1),1);
    
    HCW_constants = [c1, c2, c3, alpha]';
% Uncomment to make table with initial conditions for latex

%     for i = 1:(size(HCW_constants,2)-1)
%         Appendix_table(i,1) = i;
%         Appendix_table(i,2) = HCW_constants(1, i+1);
%         alpha = rad2deg(alpha0(i+1));
%         if alpha < 0 
%             alpha = 360 + alpha;
%         end
%         Appendix_table(i,3) = round(alpha,1);
%     end
%     % 1 & 1000 & 123 & 1 & 1000 & 123 & 1 & 1000 & 123 
%     Appendix_table(51,:) = [0,0,0];
%     for i = 1:17
%         Appendix_table_latex(i,1) = Appendix_table(i,1);
%         Appendix_table_latex(i,7) = Appendix_table(i+17,1);
%         Appendix_table_latex(i,13) = Appendix_table(i+34,1);
% 
%         Appendix_table_latex(i,3) = round(Appendix_table(i,2));
%         Appendix_table_latex(i,9) = round(Appendix_table(i+17,2));
%         Appendix_table_latex(i,15) = round(Appendix_table(i+34,2));
% 
%         Appendix_table_latex(i,5) = Appendix_table(i,3);
%         Appendix_table_latex(i,11) = Appendix_table(i+17,3);
%         Appendix_table_latex(i,17) = Appendix_table(i+34,3);
% 
%         Appendix_table_latex(i,2) = NaN;
%         Appendix_table_latex(i,4) = NaN;
%         Appendix_table_latex(i,6) = NaN;
%         Appendix_table_latex(i,8) = NaN;
%         Appendix_table_latex(i,10) = NaN;
%         Appendix_table_latex(i,12) = NaN;
%         Appendix_table_latex(i,14) = NaN;
%         Appendix_table_latex(i,16) = NaN;
%     end           
end