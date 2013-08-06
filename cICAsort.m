function [ ROIs ] = cICAsort(filename, filenameEvents)
%cICAsort perform spike sorting of high density array data
%based on convolutive ICA

% created by Christian Leibig 12.02.13

diary logfile_cICAsort.txt

% dbstop in cICAsort.m at 192

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[dataset] = read_data(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data and array specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sensor_rows = dataset.sensorRows;
sensor_cols = dataset.sensorCols;
sr = dataset.sr;
frameStartTimes = dataset.frameStartTimes;
pitch = dataset.sensorPitch; %in µm

d_sensor_row = double(sensor_rows(2) - sensor_rows(1)); %in coord
d_sensor_col = double(sensor_cols(2) - sensor_cols(1)); %in coord
d_row = d_sensor_row * pitch; %in µm
d_col = d_sensor_col * pitch; %in µm

sensor_rho = 1000000/(d_row * d_col); %per mm²

data = dataset.X;

clear dataset


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%general:
plotting =  0;
interactive = 0;

%Tissue specs:
neuron_rho = 1000; %in mm⁻²
if interactive
    neuron_rho = input('Please specify the expected neuron density in mm⁻²: ');
end
    
%ROI segmentation:
par.roi.method = 'cog';
par.roi.minNoEvents = 3;
if d_sensor_col == 2 && d_sensor_row == 1
    par.roi.thr_factor = 10.95;
    par.roi.n_rows = 5;
    par.roi.n_cols = 3;
    par.roi.n_frames = 3;
elseif d_sensor_col == d_sensor_row
    par.roi.thr_factor = 9.5;
    par.roi.n_rows = 3;
    par.roi.n_cols = 3;
    par.roi.n_frames = 3;
else
    error(strcat('Unknown readout configuration.\n',...
    'Please check parameters for ROI identification'));
end
par.roi.horizon = floor(sr);%~1 ms to the left and to the right
%of detected activity is taken for the temporal ROI

%fastICA parameters:

par.ica.allframes = true; %if true, all frames are used
par.ica.allchannels = false; %if true, all channels are used
%par.ica.frames - boolean array of length T, indicating the frames to
%be used (overwritten, if allframes is true) gets assigned after roi
%identification
%par.ica.channels - boolean array of length M, indicating the channels
%to be used (overwritten, if allchannels is true) gets assigned after
%roi identification
par.ica.nonlinearity = 'pow3';
par.ica.estimate = false;
par.ica.cpn  = 1; %components per neuron for later use to calculate
%par.ica.numOfIC (overwritten, if params.estimate is true)
par.ica.per_var = 1; %keep that many dimensions such that per_var of the
%total variance gets explained
par.ica.approach = 'symm';
par.ica.verbose = 'off';

%convolutive ICA:
if interactive
    do_cICA = input('Do you want to perform convolutive ICA(1)? (0, otherwise)?');
else
   do_cICA = false;
end
 
if (round(sr) <= 12) && (round(sr) >= 11)
    L = 7; M = 0;
elseif (round(sr) <= 24) && (round(sr) >= 23)
    L = 8;M = 12;
else
    L = input('Please specify L: ');
    M = input('Please specify M: ');    
end

allframes_cica = 1;
min_skewness = 0.2;
d_max = 1000; %maximal distance in \mum for extrema of component filters
min_corr = 0.1;
grouping = 'cluster';
max_cluster_size = 4;
max_iter = 10;
min_no_peaks = 5;
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
% ROI segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metaData.sensor_rows = sensor_rows;
metaData.sensor_cols = sensor_cols;
metaData.sr = sr;
metaData.frameStartTimes = frameStartTimes;
metaData.filename_events = filenameEvents;

[ROIs, OL] = roisegmentation(data, metaData, par.roi, plotting);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate treatment of each region of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(ROIs)
    
    fprintf('\nWorking on ROI %g out of %g\n\n',i,length(ROIs));
    
    X = ROIs(i).X;
    sensor_rows_roi = ROIs(i).sensor_rows;
    sensor_cols_roi = ROIs(i).sensor_cols;
    T_mask = ROIs(i).T_mask;
    N_mask = ROIs(i).N_mask;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing with fastICA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    par.ica.frames = T_mask;
    par.ica.channels = N_mask;

    par.ica.numOfIC = ceil(par.ica.cpn/(sensor_rho/neuron_rho) * ...
        nnz(par.ica.channels));

    [S, A, W, par.ica] = fasticanode(X,par.ica);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convolutive ICA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_cICA
        %Initialize lagged filters:
        A_tau = zeros(size(X,1),size(A,2));
        A_tau(N_mask,:) = A;
        A_tau(:,:,2:L+1) = 0;

        if allframes_cica
            frames_ROI_cica = true(size(X,2),1);
        else
            frames_ROI_cica = T_mask;
        end

        %Perform convolutive ICA:
        t1 = clock;
        [S, A_tau, S_noise, A_noise] = ConvolutiveICA(S,L,A_tau,sr,...
            d_row,d_col,length(sensor_rows_roi),length(sensor_cols_roi),d_max,...
            frames_ROI_cica,do_cICA,'M',M,'maxlags',maxlags,...
            'plotting',plotting,'min_skewness',min_skewness,'min_corr',min_corr,...
            'approach',grouping,'max_cluster_size',max_cluster_size,...
            'max_iter',max_iter,'min_no_peaks',min_no_peaks,...
            't_s',t_s,'t_jitter',t_jitter, 'coin_thr',coin_thr);
        t2 = clock;
        fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));
    else
        fprintf('Convolutive ICA is not applied!\n');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for noisy components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if convolutive ICA was not applied, we need to initialize some variables:
    if ~exist('A_noise','var'); A_noise = []; end
    if ~exist('S_noise','var'); S_noise = []; end
    if ~exist('A_tau','var')
        A_tau = zeros(size(X,1),size(A,2));
        A_tau(N_mask,:) = A;
        A_tau(:,:,2:L+1) = 0;
    end

    [keep] = checkfornoisycomponents(S,min_skewness,min_no_peaks,sr,plotting);

    % store noisy stuff away and remove it from components and filters:
    S_noise = [S_noise;S(~keep,:)];
    A_noise = cat(2,A_noise,A_tau(:,~keep,:));
    S = S(keep,:);
    A_tau = A_tau(:,keep,:);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spike time identification
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % t1 = clock;
    % [units] = SpikeTimeIdentificationKlustaKwik(S_cica, sr,1);
    % t2 = clock;
    % fprintf('Spike time identification and clustering with KlustaKwik\n');
    % fprintf('performed in %g seconds\n',etime(t2,t1));

    fprintf('Spike time identification with Hartigans dip test\n');
    t1 = clock;
    [units] = SpikeTimeIdentificationHartigan(S, sr,sign_lev,plotting,interactive);
    t2 = clock;
    fprintf('performed in %g seconds\n',etime(t2,t1));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collecting preliminary results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data_tmp = reshape(X,...
        [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);
    for k = 1:length(units)
        units(k).A_tau = A_tau(:,k,:);
        %Consider to calculate STAs only based on "non-coincident" spikes!
        units(k).STA = GetSTA(data_tmp, units(k).time, sr, 0);
        [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
            == max(max(max(abs(units(k).STA)))));
        units(k).boss_row = sensor_rows_roi(row_max);
        units(k).boss_col = sensor_cols_roi(col_max);
    end
    clear data_tmp


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove duplicates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t1 = clock;
    fprintf('Checking for duplicates...\n');
    [duplicate_pairs] = checkforintraroiduplicates(units, sr, ...
        t_s, t_jitter, coin_thr, sim_thr, plotting, interactive);
    N_dupl = size(duplicate_pairs,1);
    t2 = clock;
    fprintf('found %g intraregional duplicates in %g seconds\n',...
        N_dupl,etime(t2,t1));

    if ~isempty(duplicate_pairs)
        remove = false(length(units),1);
        remove(duplicate_pairs(:,1)) = true;
        units_dupl = units(remove);
        S_dupl = S(remove,:);
        A_dupl = [];
        A_dupl = cat(2,A_dupl,units(remove).A_tau);
        units = units(~remove);
        S = S(~remove,:);
        A_tau = A_tau(:,~remove,:);
        clear remove
    else
        units_dupl = [];
        S_dupl = [];
        A_dupl = [];
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choose finally accepted units
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Just for visualization purposes:
    if plotting
        SpikeTimeIdentificationHartigan(S, sr,sign_lev,1,0);
    end

    if interactive
        fprintf('Please examine found units.\n');
        keyboard;
        remove_IDs = input(...
         'Please specify all units to delete by their IDs in a vector:\n');
        
        remove = false(length(units),1);
        remove(remove_IDs) = true;
        
        S_noise = [S_noise;S(remove,:)];
        A_noise = cat(2,A_noise,units(remove).A_tau);
        S = S(~remove,:);
        units = units(~remove);
        A_tau = A_tau(:,~remove,:);
        clear remove
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collecting final results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if size(A_tau,2) ~= size(S,1)
        error('Something went wrong!');
    end

    %This time with interactive thresholding:
    if interactive
        sign_lev = input('Specify significance level for unimodality:');
    end
    
    [units] = SpikeTimeIdentificationHartigan(S, sr,sign_lev,plotting,interactive);

    data_tmp = reshape(X,...
        [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);
    for k = 1:length(units)
        units(k).A_tau = A_tau(:,k,:);
        %Consider to calculate STAs only based on "non-coincident" spikes!
        units(k).STA = GetSTA(data_tmp, units(k).time, sr, 0);
        [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
            == max(max(max(abs(units(k).STA)))));
        units(k).boss_row = sensor_rows_roi(row_max);
        units(k).boss_col = sensor_cols_roi(col_max);
    end
    clear data_tmp

    if interactive
        fprintf('Another chance to inspect final results before they are saved.\n');
        keyboard;
    end
    
    ROIs(i).A = A;
    ROIs(i).A_dupl = A_dupl;
    ROIs(i).A_noise = A_noise;
    ROIs(i).A_tau = A_tau;
    ROIs(i).S = S;
    ROIs(i).S_dupl = S_dupl;
    ROIs(i).S_noise = S_noise;
    ROIs(i).W = W;
    ROIs(i).duplicate_pairs = duplicate_pairs;
    ROIs(i).units = units;
    ROIs(i).units_dupl = units_dupl;
    
    clear A A_dupl A_noise A_tau S S_dupl S_noise W duplicate_pairs units units_dupl
    
    clear X sensor_rows_roi sensor_cols_roi T_mask N_mask
    

end %end loop over regions of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine results of different ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nCombining results of different regions of interest...\n');
if length(ROIs) > 1
    t1 = clock;
    [ROIs, N_INTER_ROI_DUPL] = combinerois(ROIs, OL, sr, data,...
        sensor_rows, sensor_cols, t_s, t_jitter, coin_thr, sim_thr,...
        plotting, interactive);
    t2 = clock;    
    fprintf('found %g interregional duplicates in %g seconds\n',...
        N_INTER_ROI_DUPL,etime(t2,t1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResults(filename, ROIs);   
    
t_total_2 = clock;
fprintf('Total cICAsort performed in %g seconds\n',etime(t_total_2,t_total_1));

diary off







