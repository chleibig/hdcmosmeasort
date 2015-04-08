function [row_ind, col_ind] = getfilterextrema(A, rowlist, collist, sensor_coord)
% [row_ind, col_ind] = getfilterextrema(A, rowlist, collist, sensor_coord)
% Get the spatial coordinates of the absolute extremum of the spatiotemporal
% mixing matrix A(:,i,:) for each source i 

% Input:

% A : (N,D,L+1) multidimensional array containing the zero-lag plus
% L lagged mixing matrices for N sensors and D sources respectively

% rowlist, collist: the original sensor coordinates

D = size(A,2);
N_row = length(rowlist);
N_col = length(collist);
row_ind = zeros(D,1);
col_ind = zeros(D,1);
for i = 1:D
    [value,frame_ind] = max(max(abs(A(:,i,:))));
    [row,col] = find(reshape(abs(A(:,i,frame_ind)), [N_row N_col]) == value );
    if sensor_coord
        row_ind(i) = rowlist(row);
        col_ind(i) = collist(col);
    else
        row_ind(i) = row;
        col_ind(i) = col;
    end 
end 















end