clear all;

F=[5; 5; 5; 5; 5]; %each of the 5 sats has this amount of fuel

load('Reconfiguration_matrix.mat')

for k = 1:3
    [A(:,:,k), Fleft] = maneuverAssignment(Reconfiguration_matrix, F);
    disp(A(:,:,k));
    
    [~, ind] = min(Fleft);
    
    Reconfiguration_matrix(A(ind,1,k),A(ind,2,k)) = sum(F);
    disp(Reconfiguration_matrix);
end