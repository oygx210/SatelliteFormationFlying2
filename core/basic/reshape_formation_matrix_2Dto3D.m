function Formation_rv_3D = reshape_formation_matrix_2Dto3D(Formation_rv_2D)

[a, b] = size(Formation_rv_2D);
N_sats = a/6;
t_vec = b;
Formation_rv_3D = zeros(6, t_vec, N_sats);


for i = 1:N_sats
    Formation_rv_3D(:,:,i) = Formation_rv_2D(i*6-5:i*6,:);
end

end