function [ X_ROI, sensor_rows_ROI, sensor_cols_ROI, frames_ROI, act_chs ] = ...
                ROIIdentification(data, sensor_rows, sensor_cols,options,plotting)
%ROIIDENTIFICATION Extract region of interest from data
%
%Input:
%
%   data: (N_row x N_col x N_samples) array
%   sensor_rows: vector of sensor row IDs
%   sensor_cols: vector of sensor col IDs
%   options: struct containing parameters for hyperspherical threshold
%            detection:
%            options.thr_factor threshold in multiples of std
%            options.n_rows, .n_cols, .n_frames  environment in data coor-
%                                                dinates
%            options.horizon  frames to the past and to the future of detected 
%                             activity is taken for the temporal ROI
%   plotting - flag to switch graphical output on or of
%
%Output:
%
%   X: (N_sensors x N_samples) array containing a single rectangular ROI 
%   sensor_rows_ROI, sensor_cols_ROI: list of respective sensor coordinates
%                                      in ROI

min_act = 3;

if ~strcmp(class(data),'double');
    V = double(data);
else
    V = data;
end

[N_row, N_col, N_frames] = size(V);

thr_factor = options.thr_factor;
n_rows = options.n_rows;
n_cols = options.n_cols;
n_frames = options.n_frames;
dim_env = n_rows * n_cols * n_frames;
if ~all([mod(n_rows,2) mod(n_cols,2) mod(n_frames,2)])
    error('Environment must be odd in each dimension.');
end
drow = (n_rows - 1)/2;
dcol = (n_cols - 1)/2;
dframe = (n_frames - 1)/2;

t1 = clock;

%mean-correction and normalization:
V = bsxfun(@minus,V,mean(V,3));
V = bsxfun(@rdivide,V,median(abs(V)/0.6745,3));
V_quad = V.^2;
clear V;

if dim_env == 1
    tcs = (V_quad > thr_factor^2);
    clear V_quad
else
    %Expensive computation only necessary for multi-dim environment...
    
    V_env = zeros(N_row, N_col, N_frames);
    
    fprintf('Performing activity detection in %g-dim environment...',dim_env);

    % for loops are comparably performant as vectorization with
    % arrayfun (see below):
    % for row = 1+drow:N_row-drow
    %     for col = 1+dcol:N_col-dcol
    %         for frame = 1+dframe:N_frames-dframe
    %             V_env(row,col,frame) = sum(reshape(V_quad(row-drow:row+drow,col-dcol:col+dcol,...
    %                         frame-dframe:frame+dframe),dim_env,1));
    %
    %         end
    %     end
    % end
    
    [rowg,colg,frameg] = ndgrid(1+drow:N_row-drow,1+dcol:N_col-dcol,...
        1+dframe:N_frames-dframe);
    V_env(1+drow:N_row-drow,1+dcol:N_col-dcol,1+dframe:N_frames-dframe) = ...
        arrayfun(@(row,col,frame)(sum(reshape(V_quad(row-drow:row+drow,...
        col-dcol:col+dcol,frame-dframe:frame+dframe),dim_env,1))),rowg,colg,frameg);

    clear V_quad
    
    tcs = (V_env > thr_factor^2);
    V_env_pks = V_env(tcs);
    if plotting
       figure;hist(sqrt(V_env_pks(:)),floor(sqrt(length(V_env(:)))));
       title('hist(|V_{env}(tcs)|)');
       xlabel('V_{env}');
       ylabel('counts');
    end
end

t2 = clock;
fprintf('...done in %g seconds\n',etime(t2,t1));

activity = sum(tcs,3);
[rows, cols] = find(activity >= min_act);
%tcs_per_frame = squeeze(sum(sum(tcs)));
[ frames_ROI ] = fillvector( squeeze(sum(sum(tcs))) >= 1, options.horizon );

%ROI selection in space:
%first set those channels to zero that do not reach the activity threshold:
%data( repmat(activity < min_act,[1 1 N_frames]) ) = 0;
%mask for active channels: 
act_chs = (activity >= min_act);
act_chs = act_chs(min(rows):max(rows),min(cols):max(cols));
act_chs = reshape(act_chs, [size(act_chs,1)*size(act_chs,2) 1]);
fprintf('Found activity on %g of %g input channels.\n',nnz(act_chs),N_row*N_col);
%select rectangular box around activity channels:
data_ROI = data(min(rows):max(rows),min(cols):max(cols),:);
X_ROI = reshape(data_ROI, [size(data_ROI,1)*size(data_ROI,2) size(data_ROI,3)]);
sensor_rows_ROI = sensor_rows(min(rows):max(rows));
sensor_cols_ROI = sensor_cols(min(cols):max(cols));


if plotting
    figure;imagesc(activity);axis('image');colorbar;
    hold on; plot(cols, rows,'go');
    plot([min(cols)-0.5 max(cols)+0.5], [min(rows)-0.5 min(rows)-0.5],'w');
    plot([min(cols)-0.5 max(cols)+0.5], [max(rows)+0.5 max(rows)+0.5],'w');
    plot([min(cols)-0.5 min(cols)-0.5], [min(rows)-0.5 max(rows)+0.5],'w');
    plot([max(cols)+0.5 max(cols)+0.5], [min(rows)-0.5 max(rows)+0.5],'w');
end


end