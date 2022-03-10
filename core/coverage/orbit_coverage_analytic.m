function [i_opt, j_opt, cov_tot, r_sat] = orbit_coverage_analytic(r_ECEF, t_vec, Population_map, consts)

    % sphere generation
    [M, N] = size(Population_map);
    phi_range = linspace(-pi, pi, N+1);
    phi_range = phi_range(1:end-1);
    th_range = linspace(pi, 0, M+1);
    th_range = th_range(1:end-1);
    [PHI, TH] = meshgrid(phi_range, th_range);

    % trigonometric matrices
    cTH = cos(TH);
    sTH = sin(TH);
    cPHI = cos(PHI);
    sPHI = sin(PHI);
    oe(1) = vecnorm(r_ECEF(:,1));
    b = oe(1)/consts.rEarth;
    R_ij = consts.rEarth*cat(3, sTH.*cPHI, sTH.*sPHI, cTH); % coordinates of grid nodes
    
    % beam and satellite view parameters
    th_beam = deg2rad(bw/2);
    th_max = asin(1/b);
    th_eff = th_max - th_beam;
    gam_max = pi/2 - th_max; % angle along earth surface
    gam_eff = -th_eff + asin(b*sin(th_eff));
    gam_beam = -th_beam + asin(b*sin(th_beam));
    cos_gam = cos(gam_max);
    cos_eff = cos(gam_eff);
    cos_beam = cos(gam_beam);

    rs = r_ECEF;
    % optimization
    
    % worst-case pre-allocation of dynamic arrays for speed
    view_crop = zeros(size(Population_map));
    search_view = zeros(size(Population_map));
    dist2 = zeros(size(Population_map));
    slant = zeros(size(Population_map));
    R_beam = zeros([size(Population_map), 3]);
    C_crop = zeros(size(Population_map));
    beam = zeros(size(Population_map));     
    
    cov_tot = 0;
    for i = 1:
       
        r_sat = reshape(rs(:,i), [1,1,3]);
        R_sat = repmat(r_sat, size(PHI));
        R_beam = R_ij - R_sat;
        
        % create satellite view
        val_view = dot(R_ij, R_sat, 3)/(consts.rEarth*oe(1));
        view = (val_view >= cos_gam);
        
        % array cropping 
        idx1 = find(sum(view, 2));
        idx2 = find(sum(view, 1));
        view_crop = view(idx1, idx2);
        search_view = (val_view(idx1, idx2) >= cos_eff);
        dist2 = sum(R_beam(idx1, idx2, :).^2, 3);
        slant = dist2/(oe(1) - consts.rEarth)^2;
        R_beam = R_beam(idx1, idx2, :)./repmat(sqrt(dist2), [1,1,3]); % normalization
        C_crop = Population_map(idx1, idx2);
        
        % The piece of code is used for checking whether a POI is covered
        % by spacecraft beam
%         cov_opt = 0;
%         i_opt = 0;
%         j_opt = 0;
%         % searching for a beam direction with the biggest number of views
%         % - to be changed
%         for i_map = 1:length(idx1)
%             for j_map = 1:length(idx2)
%                 if search_view(i_map, j_map)               
%                     r_cen = R_beam(i_map, j_map,:);
%                     beam = view_crop & ...
%                         (R_beam(:,:,1)*r_cen(1) + R_beam(:,:,2)*r_cen(2) + R_beam(:,:,3)*r_cen(3) > cos(th_beam));
%                     cov = sum(sum(C_crop(beam)./slant(beam)));
%                     if cov >= cov_opt
%                         cov_opt = cov;
%                         i_opt = idx1(i_map);
%                         j_opt = idx2(j_map);
%                     end
%                 end
%             end
%         end
%         cov_tot = cov_tot + cov_opt;
%     end
    
end