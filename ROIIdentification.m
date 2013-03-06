function [ X_ROI, sensor_rows_ROI, sensor_cols_ROI ] = ROIIdentification(data,...
                                       sensor_rows, sensor_cols,thr_factor)
%ROIIDENTIFICATION Extract region of interest from data
%
%Input:
%
%   data: (N_row x N_col x N_samples) array
%
%Output:
%
%   X: (N_sensors x N_samples) array containing a single rectangular ROI 
%   sensor_rows_ROI, sensor_cols_ROI: list of respective sensor coordinates
%                                      in ROI

[N_row, N_col, N_samples] = size(data);
X = reshape(data, [N_row*N_col N_samples]);

% one-dim-thresholding:
stds_hat = median(abs(X)/0.6745,2);
tcs = bsxfun(@(x,y)(abs(x) > thr_factor * y), X, stds_hat);
activity = reshape(sum(tcs,2), [N_row N_col]);
[rows, cols] = find(activity > 1);

%ROI selection in space:
data_ROI = data(min(rows):max(rows),min(cols):max(cols),:);
X_ROI = reshape(data_ROI, [size(data_ROI,1)*size(data_ROI,2) size(data_ROI,3)]);
sensor_rows_ROI = [sensor_rows(min(rows):max(rows))];
sensor_cols_ROI = [sensor_cols(min(cols):max(cols))];

%Plotting
figure;imagesc(activity);axis('image');colorbar;
hold on; plot(cols, rows,'go');
plot([min(cols)-0.5 max(cols)+0.5], [min(rows)-0.5 min(rows)-0.5],'w');
plot([min(cols)-0.5 max(cols)+0.5], [max(rows)+0.5 max(rows)+0.5],'w');
plot([min(cols)-0.5 min(cols)-0.5], [min(rows)-0.5 max(rows)+0.5],'w');
plot([max(cols)+0.5 max(cols)+0.5], [min(rows)-0.5 max(rows)+0.5],'w');


end