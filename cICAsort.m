function [ X_ROI,sensor_rows_ROI, sensor_cols_ROI, units, S_ica, S_cica,...
                S_noise, A_noise ] = cICAsort(filename)
%cICAsort perform spike sorting of high density array data
%based on convolutive ICA

% created by Christian Leibig 12.02.13


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = hdf5load(filename);
data = permute(dataset.Data, [3 2 1]);
if ~strcmp(class(data),'double');data = double(data);end;
sensor_rows = dataset.Metadata.RowList;
sensor_cols = dataset.Metadata.ColumnList;
sr = length(dataset.Metadata.FrameStartTimes)/...
    dataset.Metadata.FrameStartTimes(end);%in kHz
d_sensor_row = double(dataset.Metadata.RowList(2) - dataset.Metadata.RowList(1));
d_row = d_sensor_row * 7.4;
d_sensor_col = double(dataset.Metadata.ColumnList(2) - dataset.Metadata.ColumnList(1));
d_col = d_sensor_col * 7.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ROI identification:
if d_sensor_col == 2
    options.thr_factor = 10.95;
    options.n_rows = 5;
    options.n_cols = 3;
    options.n_frames = 3;
elseif d_sensor_col == 1
    options.thr_factor = 9.5;
    options.n_rows = 3;
    options.n_cols = 3;
    options.n_frames = 3;
end
options.horizon = floor(sr/2);%~0.5 ms to the left and to the right
%of detected activity is taken for the temporal ROI

%fastICA:
nonlinearity = 'pow3';
allframes = 0;

%convolutive ICA:
if (round(sr) <= 12) && (round(sr) >= 11)
    L = 7; M = 0;
elseif (round(sr) <= 24) && (round(sr) >= 23)
    L = 8;M = 12;
end

plotting = 1;
min_skewness = 0.2;
d_max = 1000; %maximal distance in \mum for extrema of component filters
min_corr = 0.2;
approach = 'cluster';
max_cluster_size = 2;
max_iter = 10;
min_no_peaks = 2;
maxlags = 10;

%duplicates:
t_s = 0.5; %ms
t_jitter = 0.5; %ms
coin_thr = 0.5; %fraction

t_total_1 = clock;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing - bandpass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% tic;
% data = double(data);
% Wp = [ 0.3  3] * 2 / sr;
% Ws = [ 0.1 min(5,(sr/2 - 0.1))] * 2 / sr;
% [N,Wn] = buttord( Wp, [.001 .999], 3, 10);
% [B,A] = butter(N,Wn);
% [N_ROW,N_COL,N_SAMPLES] = size(data);
% data = reshape(data, [N_ROW*N_COL N_SAMPLES]);
% data_filt = cell2mat(cellfun(@(x) filtfilt(B,A,x),...
%                             num2cell(data,2),'UniformOutput',0));
% %Due to shape expected by ROIIdentification:
% data_filt = squeeze(reshape(data_filt, [N_ROW N_COL N_SAMPLES]));
% toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ X_ROI, sensor_rows_ROI, sensor_cols_ROI, frames_ROI ] = ROIIdentification(data,...
                                       sensor_rows, sensor_cols,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing with fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract only as many ICs as possible with fastICA, first determine this
%number with a deflation approach and then extract the same number
%of components with the symmetric approach:
if allframes
    t1 = clock;
    fprintf('Estimating number of instantaneous components...\n');
    [S_ica, A, W] = fastica(X_ROI,'g',nonlinearity,...
        'approach','defl','verbose','off');
    numOfIC = round(0.8 * size(S_ica,1));
    clear S_ica A W
    t2 = clock;
    fprintf('numOfIC = %g, (%g sec. elapsed.)\n',numOfIC,etime(t2,t1));

    t1 = clock;
    [S_ica, A, W] = fastica(X_ROI,'g',nonlinearity,...
        'approach','symm','numOfIC',numOfIC);
    t2 = clock;
    fprintf('Symmetric extraction of ICs performed in %g seconds\n',etime(t2,t1));

else
    t1 = clock;
    fprintf('Estimating number of instantaneous components...\n');
    [A, W] = fastica(X_ROI(:,frames_ROI),'g',nonlinearity,...
        'approach','defl','verbose','off');
    numOfIC = round(0.8 * size(W,1));
    clear S_ica A W
    t2 = clock;
    fprintf('numOfIC = %g, (%g sec. elapsed.)\n',numOfIC,etime(t2,t1));

    t1 = clock;
    [A, W] = fastica(X_ROI(:,frames_ROI),'g',nonlinearity,...
        'approach','symm','numOfIC',numOfIC);
    S_ica = W*X_ROI;%convolutive ICA has to be performed on contiguous data
    t2 = clock;
    fprintf('Symmetric extraction of ICs performed in %g seconds\n',etime(t2,t1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutive ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize lagged filters:
A_tau = A;
A_tau(:,:,2:L+1) = 0;
%Perform convolutive ICA:
t1 = clock;
[S_cica, A_tau, S_noise, A_noise] = ConvolutiveICA(S_ica,L,A_tau,sr,...
    d_row,d_col,length(sensor_rows_ROI),length(sensor_cols_ROI),d_max,...
    'M',M,'maxlags',maxlags,...
    'plotting',plotting,'min_skewness',min_skewness,'min_corr',min_corr,...
    'approach',approach,'max_cluster_size',max_cluster_size,...
    'max_iter',max_iter,'min_no_peaks',min_no_peaks,...
    't_s',t_s,'t_jitter',t_jitter, 'coin_thr',coin_thr);
t2 = clock;
fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike time identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t1 = clock;
% [units] = SpikeTimeIdentificationKlustaKwik(S_cica, sr,1);
% t2 = clock;
% fprintf('Spike time identification and clustering with KlustaKwik\n');
% fprintf('performed in %g seconds\n',etime(t2,t1));

t1 = clock;
[units] = SpikeTimeIdentificationHartigan(S_cica, sr,1);
t2 = clock;
fprintf('Spike time identification with Hartigans dip test\n');
fprintf('performed in %g seconds\n',etime(t2,t1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixing filters and spatial positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[boss_rows, boss_cols] = GetFilterExtrema(A_tau,sensor_rows_ROI,...
                                            sensor_cols_ROI,1);

for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    units(k).boss_row = boss_rows(k);
    units(k).boss_col = boss_cols(k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%results_filename = strcat(filename,'.basic_events');
%SaveResults(results_filename, units);   
    
t_total_2 = clock;
fprintf('Total cICAsort performed in %g seconds\n',etime(t_total_2,t_total_1));








