function [ ROIs, params, ROIsAsCC ] = HDCMOSMEAsort(filename, filenameEvents)
%[ ROIs, params, ROIsAsCC ] = cICAsort(filename, filenameEvents)
%cICAsort perform spike sorting of high density array data
%based on (convolutive) ICA

% created by Christian Leibig 12.02.13

diary logfile_cICAsort.txt

memory

params = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[dataset] = read_data(filename);
metadata = read_metadata(filename);

dummyData = []; %refactor: remove data from parameter list of roisegmentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data and array specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.filename = filename;
params.sensor_rows = metadata.sensorRows;
params.sensor_cols = metadata.sensorCols;
params.sr = metadata.sr;
params.frameStartTimes = metadata.frameStartTimes;
params.pitch = metadata.sensorPitch; %in µm
params.chipType = metadata.chipType;

%delta sensor coord.
d_sensor_row = double(params.sensor_rows(2) - params.sensor_rows(1));
d_sensor_col = double(params.sensor_cols(2) - params.sensor_cols(1));
%delta \mum.
params.d_row = d_sensor_row * params.pitch;
params.d_col = d_sensor_col * params.pitch;


params.sensor_rho = 1000000/(params.d_row * params.d_col); %per mm²

% data = dataset.X;
% 
% if strcmp(dataset.chipType,'NCA')
%     data = 1000 * data; %conversion to mV scale
% end
% 
% clear dataset



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%general:
params.plotting =  0;
params.interactive = 0;

%Tissue specs:
params.neuron_rho = 2029; %in mm�?�² only used if params.ica.estimate = 'none'
if params.interactive
    params.neuron_rho = input(['Please specify the expected neuron '...
                               ' density in mm�?�²: ']);
end
    
%ROI segmentation:
params.roi.method = 'cog';
params.roi.maxSensorsPerEvent = 100;
params.roi.maxSensorsPerROI = 400;
params.roi.minNoEvents = 3 * ... %multiplying factor in Spikes / second
    (params.frameStartTimes(end) - params.frameStartTimes(1))/1000;
params.roi.mergeThr = 0.8;
if d_sensor_col == 2 && d_sensor_row == 1
    params.roi.thr_factor = 10.95;
    params.roi.n_rows = 5;
    params.roi.n_cols = 3;
    params.roi.n_frames = 3;
elseif d_sensor_col == d_sensor_row
    params.roi.thr_factor = 9.5;
    params.roi.n_rows = 3;
    params.roi.n_cols = 3;
    params.roi.n_frames = 3;
else
    error(strcat('Unknown readout configuration.\n',...
    'Please check parameters for ROI identification'));
end
params.roi.horizon = 2*floor(params.sr);%~1 ms to the left and to the right
%of detected activity is taken for the temporal ROI

%fastICA parameters:
params.ica.allframes = true; %if true, all frames are used
params.ica.allchannels = false; %if true, all channels are used
%params.ica.frames - boolean array of length T, indicating the frames to
%be used (overwritten, if allframes is true) gets assigned after roi
%identification
%params.ica.channels - boolean array of length M, indicating the channels
%to be used (overwritten, if allchannels is true) gets assigned after
%roi identification
params.ica.nonlinearity = 'pow3';
%Methods for estimating the number of ICs:
%{'none','svdSpectrum','eigSpectrum','icaDeflation'}
params.ica.estimate = 'eigSpectrum'; 
params.ica.cpn  = 1; %components per neuron for later use to calculate
%params.ica.numOfIC (overwritten, if params.estimate is true)
params.ica.per_var = 1; %keep that many dimensions such that per_var of the
%total variance gets explained
params.ica.approach = 'symm';
params.ica.verbose = 'off';
params.ica.renorm = false; %if true renormalize W and S such that only noise
                       %instead of all signal is of unit variance

%convolutive ICA:
if params.interactive
    params.do_cICA = input(['Do you want to perform convolutive ICA(1)? '...
                            ' (0, otherwise)?']);
else
   params.do_cICA = false;
end

if params.do_cICA
    
    if (round(params.sr) <= 12) && (round(params.sr) >= 11)
        params.L = 7; params.M = 0;
    elseif (round(params.sr) <= 24) && (round(params.sr) >= 23)
        params.L = 8;params.M = 12;
    else
        params.L = input('Please specify L: ');
        params.M = input('Please specify M: ');
    end
    
else
    params.L = 0;params.M = 0;
end

params.allframes_cica = 1;

params.d_max = 35; %maximal distance in \mum for extrema of filters
params.min_corr = 0.1;
params.grouping = 'cluster';
params.max_cluster_size = 4;
params.max_iter = 1;
params.maxlags = params.L;

%Noise components.
params.min_no_peaks = 3 * ... %multiplying factor in Spikes / second
    (params.frameStartTimes(end) - params.frameStartTimes(1))/1000;
params.min_skewness = 0.05;

%Peak identification.
params.thrFactor = 5;
%Upsampling factor used for spike time identification
params.upsample = floor(100/params.sr);

%mixture units:
params.maxRSTD = 0.5;

%duplicates:
params.t_s = 0.5; %ms
params.t_jitter = 1; %ms
params.coin_thr = 0.3; %fraction of coincident spikes
params.sim_thr = 0.5; %similarity of average waveforms

params.sign_lev = 0.05; %for automatic threshold adaptation;


t_total_1 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing - bandpass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% tic;
% Wp = [ 0.3  3] * 2 / params.sr;
% Ws = [ 0.1 min(5,(params.sr/2 - 0.1))] * 2 / params.sr;
% [N,Wn] = buttord( Wp, [.001 .999], 3, 10);
% [B,A] = butter(N,Wn);
% [N_ROW,N_COL,N_SAMPLES] = size(data);
% data = reshape(data, [N_ROW*N_COL N_SAMPLES]);
% data_filt = cell2mat(cellfun(@(x) filtfilt(B,A,x),...
%                             num2cell(data,2),'UniformOutput',0));
% %Due to shape expected by ROIIdentification:
% data = squeeze(reshape(data_filt, [N_ROW N_COL N_SAMPLES]));
% toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metaData.sensor_rows = params.sensor_rows;
metaData.sensor_cols = params.sensor_cols;
metaData.sr = params.sr;
metaData.frameStartTimes = params.frameStartTimes;
metaData.filename_events = filenameEvents;

fprintf('\nConstructing regions of interest...\n');
t1 = clock;
[ROIs, OL, ROIsAsCC] = roisegmentation(dummyData, metaData, params.roi, params.plotting);
t2 = clock;
fprintf('...prepared %g ROIs in %g seconds\n',length(ROIs),etime(t2,t1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel processing of multiple regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrOfROIs = length(ROIs);

if nrOfROIs >= feature('numCores')
    N_SESSIONS = ceil(0.7*feature('numCores')-1);%-1 due to master process
else
    N_SESSIONS = nrOfROIs - 1;%-1 due to master process
end


if N_SESSIONS > 0
    
    multicoreDir = 'multicorefiles';
    mkdir(multicoreDir);
    startmatlabsessions(N_SESSIONS,multicoreDir);

    settings.multicoreDir = multicoreDir;
    settings.nrOfEvalsAtOnce = 1;%floor(nrOfROIs/N_SESSIONS); % default: 4
    settings.maxEvalTimeSingle = 2500;
    settings.masterIsWorker = true;
    settings.useWaitbar = true;

    % Build cell array containing all nrOfROIs parameter sets.
    parameterCell = cell(1, nrOfROIs);
    for k = 1:nrOfROIs
        ROIs(k).k = k;
        parameterCell{1,k} = {ROIs(k), params};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Actual work is done here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nParallel processing of %g different ROIs...\n',nrOfROIs);
    t1 = clock;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ROIs = cell2mat(startmulticoremaster(@processroi,parameterCell,settings));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    fprintf(['Parallel processing of %g different regions ',...
       'of interest performed in %g seconds\n'],nrOfROIs,etime(t2,t1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stopmatlabsessions(multicoreDir);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Actual work is done here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nSerially processing %g ROIs...\n',nrOfROIs);
    t1 = clock;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpROIs = cell(1, nrOfROIs);
    for i = 1:nrOfROIs
        ROIs(i).k = i;
        [ tmpROIs{i} ] = processroi( ROIs(i), params );
    end
    ROIs = cell2mat(tmpROIs);
    clear tmpROIs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    fprintf('done in %g seconds\n',etime(t2,t1));   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine results of different ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if length(ROIs) > 1
%     fprintf('\nCombining results of different regions of interest...\n');
%     t1 = clock;
%     [ROIs, N_INTER_ROI_DUPL] = combinerois(ROIs, OL, params.sr, data,...
%         params.sensor_rows, params.sensor_cols, params.t_s, params.t_jitter,...
%         params.coin_thr, params.sim_thr, params.plotting, params.interactive);
%     t2 = clock;    
%     fprintf('found %g interregional duplicates in %g seconds\n',...
%         N_INTER_ROI_DUPL,etime(t2,t1));
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SaveResults(ROIs, params);   
    
t_total_2 = clock;
fprintf('Total cICAsort performed in %g seconds\n',etime(t_total_2,t_total_1));

diary off







