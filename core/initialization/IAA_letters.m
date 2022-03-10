function formation_configuration = IAA_letters(letters, relative_orbit_type, min_ISD)

switch letters
    case "AA"
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
        plot(formation_geometry(:,1), formation_geometry(:,2), 'ok', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 12);
        axis square;

        formation_geometry(:,2) = -1*formation_geometry(:,2);
        formation_geometry = formation_geometry - formation_geometry(1,:);
        formation_geometry(14:end,1) = formation_geometry(14:end,1)  + 2.27;
        % add geometrical center
    case "IAA_logo"
        formation_geometry = [
            9.78 1.36
            9.7 1.88
            10.34 2.04
            10.49 1.45
            9.46 1.04
            9.12 0.87
            8.91 1.18];
        formation_geometry(:,2) = -1*formation_geometry(:,2);
        formation_geometry = formation_geometry - formation_geometry(1,:);
        
        phi = 0:2*pi/19:(2*pi - 2*pi/19);
        [x,y] = pol2cart(phi, 1.5*ones(1,length(phi)));
        formation_geometry(8:26,1) = x;
        formation_geometry(8:26,2) = y;
    
    case "2_sats_test"
        
        formation_geometry = [
            0 0
            min_ISD 0];
            
     otherwise 
        disp('Contact our sales manager to order a new logo for demonstration');
end

ISD_raw = 99999;
for i = 1:size(formation_geometry,1)
    for j = 1:size(formation_geometry,1)
        if i ~= j
            ISD_raw_new = sqrt((formation_geometry(i,1) - formation_geometry(j,1))^2 + (formation_geometry(i,2) - formation_geometry(j,2))^2);
            if ISD_raw_new < ISD_raw
                ISD_raw = ISD_raw_new;
            end
        end
    end
end
    
scale = min_ISD/ISD_raw;
formation_geometry = formation_geometry*scale;
formation_geometry(:,1) = -1*formation_geometry(:,1);
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
alpha = phi + (phi<0) * 2*pi;

formation_configuration = [c1, c2, c3, alpha]';

% plot(formation_geometry(:,1), formation_geometry(:,2), 'ok', 'MarkerFaceColor',  [1 1 1], 'MarkerSize', 12);
% axis square;

end