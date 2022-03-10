function [m, idx] = medoid(X)
% 
% Calculate medoid
%
% Input:
%   X: d x n data matrix
%
% Output:
%   m: medoid
%   idx: index of medoids
%
% Written by Detang Zhong (detang.zhong@canada.ca).
%
% Demo:
%  X = [12,  9, 61, 76,  2, 17, 12, 11, 26,  0;
%        65, 72,  7, 64, 21, 92, 51, 48,  9, 65;
%        39,  7, 50, 56, 29, 79, 47, 45, 10, 52;
%        70, 12, 23, 97, 86, 14, 42, 90, 15, 16;
%        13,  7,  2, 47, 80, 53, 23, 59,  7, 15;
%        83,  2, 40, 12, 22, 75, 69, 61, 28, 53]
% 
% X =
% 
%     12     9    61    76     2    17    12    11    26     0
%     65    72     7    64    21    92    51    48     9    65
%     39     7    50    56    29    79    47    45    10    52
%     70    12    23    97    86    14    42    90    15    16
%     13     7     2    47    80    53    23    59     7    15
%     83     2    40    12    22    75    69    61    28    53
% 
% m = medoid(X)
% 
% m =
% 
%     12
%     51
%     47
%     42
%     23
%     69
%

%% Determine the size of X
[d,n] = size(X);
v = dot(X,X,1);

%% Calculate Euclidean distance matrix
% D = v+v'-2*(X'*X);                  
D = -2*(X'*X); 
D = bsxfun(@plus, D, v); 
D = bsxfun(@plus, D, v');

%% Reduce chance of numerical problems
D(sub2ind([n,n],1:n,1:n)) = 0; 

%% Get the medoid and its index
md = mean(D,2);
[~,idx] = min(md);
m = X(:,idx);

end


