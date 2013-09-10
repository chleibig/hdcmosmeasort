function [ ED ] = euclideandistance(row_ind, col_ind, d_row, d_col)
% [ ED ] = euclideandistance(row_ind, col_ind, d_row, d_col)
% Computes pairwise euclidean distances of the given points
%   
% Input
% =====  
%
% row_ind, col_ind: vectors of sensor coordinates. Both must be of equal
%                   length which is the number of points
% d_row, d_col: distance between neighbouring sensors in \mum 
%
% Output
% ======
%
% ED: matrix with entry ED(i,j) specifying the distance
%           between point i and j in \mum

% christian.leibig@g-node.org, 06.09.2013

N = length(row_ind);
if N ~= length(col_ind)
    error('Lists of row and column indices must be of equal length.');
end

%Compute distance matrix:
ED = zeros(N,N);
for i=1:N
    for j=i+1:N
        ED(i,j) = sqrt(((row_ind(i) - row_ind(j))*d_row)^2 + ...
                  ((col_ind(i) - col_ind(j))*d_col)^2);
        ED(j,i) = ED(i,j);
    end
end

end