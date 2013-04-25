function [ X_ROI,sensor_rows_ROI, sensor_cols_ROI, units, S_ica, S_cica,...
   S_noise, A_noise, duplicate_pairs, units_dupl, S_dupl, A_dupl ] = ...
                                cICAsort(filename)
%cICAsort perform spike sorting of high density array data
%based on convolutive ICA

% created by Christian Leibig 12.02.13

diary logfile_cICAsort.txt


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
d_sensor_col = double(dataset.Metadata.ColumnList(2) - dataset.Metadata.ColumnList(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pitch = 7.4; %µm
d_row = d_sensor_row * pitch;
d_col = d_sensor_col * pitch;

sensor_rho = 1000000/(d_row * d_col); %per mm²

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tissue specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neuron_rho = input('Please specify the expected neuron density in mm⁻²: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ROI identification:
if d_sensor_col == 2 && d_sensor_row == 1
    options.thr_factor = 10.95;
    options.n_rows = 5;
    options.n_cols = 3;
    options.n_frames = 3;
elseif d_sensor_col == d_sensor_row
    options.thr_factor = 9.5;
    options.n_rows = 3;
    options.n_cols = 3;
    options.n_frames = 3;
else
    error(strcat('Unknown readout configuration.\n',...
    'Please check parameters for ROI identification'));
end
options.horizon = floor(sr/2);%~0.5 ms to the left and to the right
%of detected activity is taken for the temporal ROI

%fastICA:
cpn  = 3; %components per neuron
per_var = 1; %keep that many dimensions such that per_var of the total 
              %variance gets explained
nonlinearity = 'pow3';
approach = 'symm';
allframes = 0;
estimate = false;

%convolutive ICA:

if (round(sr) <= 12) && (round(sr) >= 11)
    L = 7; M = 0;
elseif (round(sr) <= 24) && (round(sr) >= 23)
    L = 8;M = 12;
end

allframes_cica = 1;
plotting = 1;
min_skewness = 0.2;
d_max = 1000; %maximal distance in \mum for extrema of component filters
min_corr = 0.1;
grouping = 'cluster';
max_cluster_size = 4;
max_iter = 10;
min_no_peaks = 3;
%maxlags = 10;
maxlags = L;

%duplicates:
t_s = 0.5; %ms
t_jitter = 1; %ms
coin_thr = 0.5; %fraction of coincident spikes
sim_thr = 0.6; %similarity of average waveforms

sign_lev = 0.05; %for automatic threshold adaptation;


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

[ X_ROI, sensor_rows_ROI, sensor_cols_ROI, frames_ROI, act_chs] = ...
            ROIIdentification(data, sensor_rows, sensor_cols,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing with fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract only as many ICs as possible with fastICA, first determine this
%number with a deflation approach and then extract the same number
%of components with the symmetric approach:
if allframes
    t1 = clock;
    fprintf('Estimating number of instantaneous components...\n');
    [S_ica, A, W] = fastica(X_ROI(act_chs,:),'g',nonlinearity,...
        'approach','defl','verbose','off');
    numOfIC = round(0.8 * size(S_ica,1));
    clear S_ica A W
    t2 = clock;
    fprintf('numOfIC = %g, (%g sec. elapsed.)\n',numOfIC,etime(t2,t1));

    t1 = clock;
    [S_ica, A, W] = fastica(X_ROI(act_chs,:),'g',nonlinearity,...
        'approach','symm','numOfIC',numOfIC);
    t2 = clock;
    fprintf('Symmetric extraction of ICs performed in %g seconds\n',etime(t2,t1));

else
    if estimate
    t1 = clock;
    fprintf('Estimating number of instantaneous components...\n');
    [A, W] = fastica(X_ROI(act_chs,frames_ROI),'g',nonlinearity,...
        'approach','defl','verbose','off','interactivePCA','off');
    numOfIC = round(0.8*size(W,1));
    clear S_ica A W
    t2 = clock;
    fprintf('numOfIC = %g, (%g sec. elapsed.)\n',numOfIC,etime(t2,t1));
    else
        numOfIC = ceil(cpn/(sensor_rho/neuron_rho) * nnz(act_chs));
    end
    %dimensionality reduction based on the percentage of variance
    %explained:
    [pcaE, pcaD] = fastica(X_ROI(act_chs, frames_ROI),'only','pca');
    d = sort(diag(pcaD),1,'descend');
    eigs_to_keep = find(cumsum(d)/sum(d) <= per_var);
    t1 = clock;
    [A, W] = fastica(X_ROI(act_chs,frames_ROI),'g',nonlinearity,...
        'approach',approach,'numOfIC',numOfIC,'lastEig',eigs_to_keep(end));
        %'pcaE', pcaE, 'pcaD', pcaD);
    S_ica = W*X_ROI(act_chs,:);%convolutive ICA has to be performed on contiguous data
    t2 = clock;
    fprintf('Extraction of ICs performed in %g seconds\n',etime(t2,t1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutive ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize lagged filters:
A_tau = zeros(size(X_ROI,1),size(A,2));
A_tau(act_chs,:) = A;
A_tau(:,:,2:L+1) = 0;

if allframes_cica
    frames_ROI_cica = true(size(X_ROI,2),1);
else
    frames_ROI_cica = frames_ROI;
end

%Perform convolutive ICA:
t1 = clock;
[S_cica, A_tau, S_noise, A_noise] = ConvolutiveICA(S_ica,L,A_tau,sr,...
    d_row,d_col,length(sensor_rows_ROI),length(sensor_cols_ROI),d_max,...
    frames_ROI_cica,'M',M,'maxlags',maxlags,...
    'plotting',plotting,'min_skewness',min_skewness,'min_corr',min_corr,...
    'approach',grouping,'max_cluster_size',max_cluster_size,...
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
[units] = SpikeTimeIdentificationHartigan(S_cica, sr,sign_lev,1,1);
t2 = clock;
fprintf('Spike time identification with Hartigans dip test\n');
fprintf('performed in %g seconds\n',etime(t2,t1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting preliminary results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_tmp = reshape(X_ROI,...
    [length(sensor_rows_ROI) length(sensor_cols_ROI) size(X_ROI,2)]);
for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    %Consider to calculate STAs only based on "non-coincident" spikes!
    units(k).STA = GetSTA(data_tmp, units(k).time, sr, 0);
    [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
                        == max(max(max(abs(units(k).STA)))));
    units(k).boss_row = sensor_rows_ROI(row_max);
    units(k).boss_col = sensor_cols_ROI(col_max);                
end
clear data_tmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove duplicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = clock;
fprintf('Checking for duplicates...\n');
[duplicate_pairs] = CheckForDuplicates(units, sr, ...
                                       t_s, t_jitter, coin_thr, sim_thr,1);
N_dupl = size(duplicate_pairs,1);
t2 = clock;
fprintf('found %g duplicates in %g seconds\n',N_dupl,etime(t2,t1));

if ~isempty(duplicate_pairs)
    remove = false(length(units),1);
    remove(duplicate_pairs(:,1)) = true;
    units_dupl = units(remove);
    S_dupl = S_cica(remove,:);
    A_dupl = [];
    A_dupl = cat(2,A_dupl,units(remove).A_tau);
    units = units(~remove);
    S_cica = S_cica(~remove,:);
    A_tau = A_tau(:,~remove,:);
    clear remove
else
    units_dupl = [];
    S_dupl = [];
    A_dupl = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose finally accepted units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Just for visualization purposes:
SpikeTimeIdentificationHartigan(S_cica, sr,sign_lev,1,0);

fprintf('Please examine found units.\n');
keyboard;
remove_IDs = input('Please specify all units to delete by their IDs in a vector:\n');

remove = false(length(units),1);
remove(remove_IDs) = true;

S_noise = [S_noise;S_cica(remove,:)];
A_noise = cat(2,A_noise,units(remove).A_tau);
S_cica = S_cica(~remove,:);
units = units(~remove);
A_tau = A_tau(:,~remove,:);
clear remove

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting final results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(A_tau,2) ~= size(S_cica,1)
    error('Something went wrong!');
end

%This time with interactive thresholding:
sign_lev_man = input('Specify significance level for unimodality:');
[units] = SpikeTimeIdentificationHartigan(S_cica, sr,sign_lev_man,1,1);

data_tmp = reshape(X_ROI,...
    [length(sensor_rows_ROI) length(sensor_cols_ROI) size(X_ROI,2)]);
for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    %Consider to calculate STAs only based on "non-coincident" spikes!
    units(k).STA = GetSTA(data_tmp, units(k).time, sr, 0);
    [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
                        == max(max(max(abs(units(k).STA)))));
    units(k).boss_row = sensor_rows_ROI(row_max);
    units(k).boss_col = sensor_cols_ROI(col_max);                
end
clear data_tmp


fprintf('Another chance to inspect final results before they are saved.\n');
keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResults(filename, units);   
    
t_total_2 = clock;
fprintf('Total cICAsort performed in %g seconds\n',etime(t_total_2,t_total_1));

diary off







