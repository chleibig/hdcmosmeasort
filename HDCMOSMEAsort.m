function [ ROIs, params, ROIsAsCC ] = HDCMOSMEAsort(filename, filenameEvents)
%[ ROIs, params, ROIsAsCC ] = HDCMOSMEAsort(filename, filenameEvents)
%HDCMOSMEAsort performs spike sorting of high density array data
%based on (convolutive) ICA

% created by Christian Leibig 12.02.13

diary logfile_HDCMOSMEAsort.txt

%memory

params = struct(); %bundles all parameters
params.filename = filename;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data and array specs from metadata of hdf5 file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata = readmetadata(filename);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.plotting = 0;
params.interactive = 0;

%%%%% Event detection is not planned to be done in Matlab - refactor that %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% ROI construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.roi.method = 'cog';
params.roi.maxSensorsPerEvent = 169;
params.roi.minNoEvents = 3 * ... %multiplying factor in Spikes / second
    (params.frameStartTimes(end) - params.frameStartTimes(1))/1000;

params.roi.mergeThr = 0.1;
params.roi.maxSensorsPerROI = 169;

params.roi.horizon = 2*floor(params.sr);%~1 ms to the left and to the right
%of detected activity is taken for the temporal ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Expected number of neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.neuron_rho = 2029; %in mm�?�² only used if params.ica.estimate = 'none'
if params.interactive
    params.neuron_rho = input(['Please specify the expected neuron '...
                               ' density in mm�?�²: ']);
end
% Upper bound for neuron number gets automatically estimated if 
% params.ica.estimate = 'eigSpectrum'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% instantaneous ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                       
%%%%% convolutive ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.do_cICA = false;
params.L = 8;
params.M = 12;
params.allframes_cica = 1;
params.min_corr = 0.02;
params.max_cluster_size = 4;
params.max_iter = 1;
%params.maxlags = params.L;
params.maxlags = params.L + params.M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Spike time identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.thrFactor = 3;
%Upsampling factor used for spike time identification
params.upsample = floor(100/params.sr);
%deprecated:
params.sign_lev = 0.05; %for automatic threshold adaptation;
% was used previously with Hartigans dip test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Automatic removal of noise sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.min_no_peaks = 5 * ... %multiplying factor in Spikes / second
    (params.frameStartTimes(end) - params.frameStartTimes(1))/1000 ...
    + (1 - normcdf(params.thrFactor))*length(params.frameStartTimes);
%previous line adds expected FPs due to Gaussian noise <-> min. Spiking
%frequency less dependent on thrFactor.
params.min_skewness = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Automatic removal of mixture units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.maxRSTD = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Selective features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.minSeparability = 0;
params.minIsoIBg = 0;
params.minIsoINN = 0;
params.minKurtosis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Fusion of results from different regions of interest %%%%%%%%%%%%%%%%
%spike train alignment.
params.t_s = 0.5; %ms
params.t_jitter = 1; %ms
%redundancy reduction parameters, adjustable in GUI.
params.d_max = 1000; %maximal distance in \mum for extrema of average waveforms
params.coin_thr = 0.5; %fraction of coincident spikes
params.sim_thr = 0.5; %similarity of average waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_total_1 = clock;% ⅰVamos!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metaData.sensor_rows = params.sensor_rows;
metaData.sensor_cols = params.sensor_cols;
metaData.sr = params.sr;
metaData.frameStartTimes = params.frameStartTimes;
metaData.filename_events = filenameEvents;

fprintf('\nConstructing regions of interest...\n');
t1 = clock;
dummyData = []; %refactor: remove data from parameter list of roisegmentation
[ROIs, OL, ROIsAsCC] = roisegmentation(dummyData, metaData, params.roi, params.plotting);
t2 = clock;
fprintf('...prepared %g ROIs in %g seconds\n',length(ROIs),etime(t2,t1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel processing of multiple regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrOfROIs = length(ROIs);

if nrOfROIs >= feature('numCores')
%     N_SESSIONS = ceil(0.8*feature('numCores')-1);%-1 due to master process
      N_SESSIONS = feature('numCores')-5;%at least -1 due to master process
else
    N_SESSIONS = nrOfROIs - 1;%-1 due to master process
end


if N_SESSIONS > 0
    
    multicoreDir = 'multicorefiles';
    mkdir(multicoreDir);
    startmatlabsessions(N_SESSIONS,multicoreDir);

    settings.multicoreDir = multicoreDir;
    settings.nrOfEvalsAtOnce = 1;%floor(nrOfROIs/N_SESSIONS); % default: 4
    settings.maxEvalTimeSingle = Inf;
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
    ROIs = cell2mat(startmulticoremaster(@processroi,...
                                         parameterCell,settings));
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
% Combine results from different ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start splitmerge from the command line to do that in a semiautomized
%fashion after termination of HDCMOSMEAsort.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%saveresults(ROIs, params);   <-> save results from GUI !
    
t_total_2 = clock;
fprintf('Total HDCMOSMEAsort performed in %g seconds\n',...
    etime(t_total_2,t_total_1));

diary off







