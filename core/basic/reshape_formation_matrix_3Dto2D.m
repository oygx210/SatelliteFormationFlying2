function rv = reshape_formation_matrix_3Dto2D(rv_in)

[a,b,c] = size(rv_in);

for i = 1:c
    rv(i*6-5:i*6,:) = rv_in(:,:,i);
end
end