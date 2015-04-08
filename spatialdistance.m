function [ SD ] = spatialdistance(A_tau, d_row, d_col, N_row, N_col)
% [ SD ] = spatialdistance(A_tau, d_row, d_col, N_row, N_col)
% Computes the pairwise euclidean distance between all components K, 
% based on the extrema of the mixing matrices A_tau(:,k,:)
%   
%   Input: 
%
%       A_tau: (N_sensors x K x L+1) one or several (lagged) mixing
%              matrices for each component K
%
%       d_row, d_col: distance between neighbouring sensors in \mum 
%
%       N_row, N_col: N_row * N_col = N_sensors
%
%   Output:
%
%       SD: (K x K) matrix with entry SD(i,j) specifying the distance
%           between component i and j in \mum

K = size(A_tau,2);

%Get extrema positions for all components in data coordinates:
[row_ind, col_ind] = getfilterextrema(A_tau,1:N_row,1:N_col,0);

%Compute distance matrix:
SD = zeros(K,K);
for i=1:K
    for j=i+1:K
        SD(i,j) = sqrt(((row_ind(i) - row_ind(j))*d_row)^2 + ...
                  ((col_ind(i) - col_ind(j))*d_col)^2);
        SD(j,i) = SD(i,j);
    end
end

end