function [ ROIs, params ] = cICAsort(filename, filenameEvents)
%cICAsort perform spike sorting of high density array data
%based on convolutive ICA

% created by Christian Leibig 12.02.13

diary logfile_cICAsort.txt

params = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[dataset] = read_data(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data and array specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.sensor_rows = dataset.sensorRows;
params.sensor_cols = dataset.sensorCols;
params.sr = dataset.sr;
params.frameStartTimes = dataset.frameStartTimes;
params.pitch = dataset.sensorPitch; %in µm

d_sensor_row = double(params.sensor_rows(2) - params.sensor_rows(1)); %in coord
d_sensor_col = double(params.sensor_cols(2) - params.sensor_cols(1)); %in coord
params.d_row = d_sensor_row * params.pitch; %in µm
params.d_col = d_sensor_col * params.pitch; %in µm

params.sensor_rho = 1000000/(params.d_row * params.d_col); %per mm²

data = dataset.X;

clear dataset



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%general:
params.plotting =  0;
params.interactive = 0;

%Tissue specs:
params.neuron_rho = 1000; %in mm⁻²
if params.interactive
    params.neuron_rho = input('Please specify the expected neuron density in mm⁻²: ');
end
    
%ROI segmentation:
params.roi.method = 'cog';
params.roi.minNoEvents = 15;
params.roi.mergeThr = 0.5;
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
params.roi.horizon = floor(params.sr);%~1 ms to the left and to the right
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
params.ica.estimate = false;
params.ica.cpn  = 1; %components per neuron for later use to calculate
%params.ica.numOfIC (overwritten, if params.estimate is true)
params.ica.per_var = 1; %keep that many dimensions such that per_var of the
%total variance gets explained
params.ica.approach = 'symm';
params.ica.verbose = 'off';
params.ica.renorm = true; %if true renormalize W and S such that only noise
                       %instead of all signal is of unit variance

%convolutive ICA:
if params.interactive
    params.do_cICA = input('Do you want to perform convolutive ICA(1)? (0, otherwise)?');
else
   params.do_cICA = false;
end

if params.do_cICA
    
    if (round(sr) <= 12) && (round(sr) >= 11)
        params.L = 7; params.M = 0;
    elseif (round(sr) <= 24) && (round(sr) >= 23)
        params.L = 8;params.M = 12;
    else
        params.L = input('Please specify L: ');
        params.M = input('Please specify M: ');
    end
    
else
    params.L = 0;params.M = 0;
end

params.allframes_cica = 1;
params.min_skewness = 0.2;
params.d_max = 1000; %maximal distance in \mum for extrema of component filters
params.min_corr = 0.1;
params.grouping = 'cluster';
params.max_cluster_size = 4;
params.max_iter = 10;
params.min_no_peaks = 3;
params.maxlags = params.L;


%mixture units:
params.maxRSTD = 0.5;

%duplicates:
params.t_s = 0.5; %ms
params.t_jitter = 1; %ms
params.coin_thr = 0.5; %fraction of coincident spikes
params.sim_thr = 0.8; %similarity of average waveforms

params.sign_lev = 0.05; %for automatic threshold adaptation;


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

metaData.sensor_rows = params.sensor_rows;
metaData.sensor_cols = params.sensor_cols;
metaData.sr = params.sr;
metaData.frameStartTimes = params.frameStartTimes;
metaData.filename_events = filenameEvents;

[ROIs, OL] = roisegmentation(data, metaData, params.roi, params.plotting);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel processing of multiple regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrOfROIs = length(ROIs);

if nrOfROIs >= feature('numCores')
    N_SESSIONS = round(0.8*feature('numCores')-1);%-1 due to master process
else
    N_SESSIONS = nrOfROIs - 1;%-1 due to master process
end
%N_SESSIONS = 5;



if N_SESSIONS > 0
    
    multicoreDir = 'multicorefiles';
    mkdir(multicoreDir);
    startmatlabsessions(N_SESSIONS,multicoreDir);

    
    settings.multicoreDir = multicoreDir;
    settings.nrOfEvalsAtOnce = 2;%floor(nrOfROIs/N_SESSIONS); % default: 4
    settings.maxEvalTimeSingle = 500;
    settings.masterIsWorker = true;
    settings.useWaitbar = false;

    % Set handle to postprocessing function
    % settings.postProcessHandle   = @postprocessdemo
    % settings.postProcessUserData = datestr(now);

    % Build cell array containing all nrOfROIs parameter sets.
    parameterCell = cell(1, nrOfROIs);
    for k = 1:nrOfROIs
        ROIs(k).k = k;
        parameterCell{1,k} = {ROIs(k), params};
    end


    % Call function STARTMULTICOREMASTER.
    fprintf('\nParallel processing of %g different regions of interest...\n',nrOfROIs);
    t1 = clock;
    ROIs = cell2mat(startmulticoremaster(@processroi, parameterCell, settings));
    t2 = clock;
    fprintf(strcat('Parallel processing of %g different regions ',...
       'of interest performed in %g seconds\n'),nrOfROIs,etime(t2,t1));

    stopmatlabsessions(multicoreDir);
    %rmdir(multicoreDir);
    
else
    fprintf('\nNo need to parallelize, processing %g region of interest...\n',nrOfROIs);
    t1 = clock;
    ROIs(1).k = 1;
    [ ROIs ] = processroi( ROIs, params );
    t2 = clock;
    fprintf('done in %g seconds\n',etime(t2,t1));   
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate treatment of each region of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(ROIs)
%     
%     fprintf('\nWorking on ROI %g out of %g\n\n',i,length(ROIs));
%     
%     X = ROIs(i).X;
%     sensor_rows_roi = ROIs(i).sensor_rows;
%     sensor_cols_roi = ROIs(i).sensor_cols;
%     T_mask = ROIs(i).T_mask;
%     N_mask = ROIs(i).N_mask;
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Preprocessing with fastICA
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     params.ica.frames = T_mask;
%     params.ica.channels = N_mask;
% 
%     params.ica.numOfIC = ceil(params.ica.cpn/(sensor_rho/neuron_rho) * ...
%         nnz(params.ica.channels));
% 
%     [S, A, W, params.ica] = fasticanode(X,params.ica);
%     
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Convolutive ICA
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if params.do_cICA
%         %Initialize lagged filters:
%         A_tau = zeros(size(X,1),size(A,2));
%         A_tau(N_mask,:) = A;
%         A_tau(:,:,2:params.L+1) = 0;
% 
%         if params.allframes_cica
%             frames_ROI_cica = true(size(X,2),1);
%         else
%             frames_ROI_cica = T_mask;
%         end
% 
%         %Perform convolutive ICA:
%         t1 = clock;
%         [S, A_tau, S_noise, A_noise] = ConvolutiveICA(S,params.L,A_tau,params.sr,...
%             params.d_row,params.d_col,length(sensor_rows_roi),length(sensor_cols_roi),params.d_max,...
%             frames_ROI_cica,params.do_cICA,'M',params.M,'maxlags',params.maxlags,...
%             'plotting',params.plotting,'min_skewness',params.min_skewness,'min_corr',params.min_corr,...
%             'approach',params.grouping,'max_cluster_size',params.max_cluster_size,...
%             'max_iter',params.max_iter,'min_no_peaks',params.min_no_peaks,...
%             't_s',params.t_s,'t_jitter',params.t_jitter, 'coin_thr',params.coin_thr);
%         t2 = clock;
%         fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));
%     else
%         fprintf('Convolutive ICA is not applied!\n');
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Check for noisy components
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %if convolutive ICA was not applied, we need to initialize some variables:
%     if ~exist('A_noise','var'); A_noise = []; end
%     if ~exist('S_noise','var'); S_noise = []; end
%     if ~exist('A_tau','var')
%         A_tau = zeros(size(X,1),size(A,2));
%         A_tau(N_mask,:) = A;
%         A_tau(:,:,2:params.L+1) = 0;
%     end
% 
%     [keep] = checkfornoisycomponents(S,params.min_skewness,params.min_no_peaks,params.sr,params.plotting);
% 
%     % store noisy stuff away and remove it from components and filters:
%     S_noise = [S_noise;S(~keep,:)];
%     A_noise = cat(2,A_noise,A_tau(:,~keep,:));
%     S = S(keep,:);
%     A_tau = A_tau(:,keep,:);
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Spike time identification
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     fprintf('Spike time identification and clustering with KlustaKwik\n');
%     t1 = clock;
%     [units] = SpikeTimeIdentificationKlustaKwik(S,0,10, params.sr, params.plotting);
%     t2 = clock;
%     fprintf('performed in %g seconds\n',etime(t2,t1));
%     
% %     fprintf('Spike time identification with Hartigans dip test\n');
% %     t1 = clock;
% %     [units] = SpikeTimeIdentificationHartigan(S, params.sr,params.sign_lev,params.plotting,1);
% %     t2 = clock;
% %     fprintf('performed in %g seconds\n',etime(t2,t1));
% 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Remove mixed units
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if ~exist('A_mix','var'); A_mix = []; end
%     if ~exist('S_mix','var'); S_mix = []; end
%     if ~exist('units_mix','var'); units_mix = []; end
%     
%     if length(units) > 0
% 
%         clear keep
%         keep = ([units.RSTD] <= params.maxRSTD);
%         %dbstop in cICAsort.m at 276 if (nnz(~keep) > 0)
%         
% 
%         S_mix = [S_mix;S(~keep,:)];
%         A_mix = cat(2,A_mix,A_tau(:,~keep,:));
%         units_mix = units(~keep);
% 
%         S = S(keep,:);
%         A_tau = A_tau(:,keep,:);
%         units = units(keep);
% 
%         fprintf('Removed %g units supposed to contain mixtures.\n',nnz(~keep));
% 
%         clear keep
% 
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Collecting preliminary results
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     data_tmp = reshape(X,...
%         [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);
%     for k = 1:length(units)
%         units(k).A_tau = A_tau(:,k,:);
%         %Consider to calculate STAs only based on "non-coincident" spikes!
%         units(k).STA = GetSTA(data_tmp, units(k).time, params.sr, 0);
%         [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
%             == max(max(max(abs(units(k).STA)))));
%         units(k).boss_row = sensor_rows_roi(row_max);
%         units(k).boss_col = sensor_cols_roi(col_max);
%     end
%     clear data_tmp
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Remove duplicates
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     t1 = clock;
%     fprintf('Checking for duplicates...\n');
%     [duplicate_pairs] = checkforintraroiduplicates(units, params.sr, ...
%         params.t_s, params.t_jitter, params.coin_thr, params.sim_thr, params.plotting, params.interactive);
%     N_dupl = size(duplicate_pairs,1);
%     t2 = clock;
%     fprintf('found %g intraregional duplicates in %g seconds\n',...
%         N_dupl,etime(t2,t1));
% 
%     %dbstop in cICAsort.m at 326 if (N_dupl > 0)
%     %Experiment with additional criteria to decide upon which duplicate
%     %partner to remove:
%     for d = 1:N_dupl
% %         SpikeTimeIdentificationKlustaKwik(S(duplicate_pairs(d,:),:),0,10, params.sr, 1);
% %         if (units(duplicate_pairs(d,1)).RSTD > 1.5*units(duplicate_pairs(d,2)).RSTD)
% %             %duplicate_pairs(d,1) is considered to be a mixture and will be
% %             %removed
% %             break;
% %         end
% %         if (units(duplicate_pairs(d,2)).RSTD > 1.5*units(duplicate_pairs(d,1)).RSTD)
% %             %duplicate_pairs(d,1) is considered to be a mixture and will be
% %             %removed
% %             duplicate_pairs(d,:) = duplicate_pairs(d,end:-1:1);
% %             break;
% %         end
%         
%         %No mixture detected - the unit with higher separability will
%         %be kept
%         if units(duplicate_pairs(d,1)).separability <= units(duplicate_pairs(d,2)).separability
%             %remove the first
%         else
%             %remove the second
%             duplicate_pairs(d,:) = duplicate_pairs(d,end:-1:1);
%         end
%     end
%     
%     
%     
%     if ~isempty(duplicate_pairs)
%         remove = false(length(units),1);
%         remove(duplicate_pairs(:,1)) = true;
%         units_dupl = units(remove);
%         S_dupl = S(remove,:);
%         A_dupl = [];
%         A_dupl = cat(2,A_dupl,units(remove).A_tau);
%         units = units(~remove);
%         S = S(~remove,:);
%         A_tau = A_tau(:,~remove,:);
%         clear remove
%     else
%         units_dupl = [];
%         S_dupl = [];
%         A_dupl = [];
%     end
% 
% 
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % Choose finally accepted units
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     
% %     %Just for visualization purposes:
% %     if params.plotting
% %         SpikeTimeIdentificationHartigan(S, sr,sign_lev,1,0);
% %     end
% % 
% %     if params.interactive
% %         fprintf('Please examine found units.\n');
% %         keyboard;
% %         remove_IDs = input(...
% %          'Please specify all units to delete by their IDs in a vector:\n');
% %         
% %         remove = false(length(units),1);
% %         remove(remove_IDs) = true;
% %         
% %         S_noise = [S_noise;S(remove,:)];
% %         A_noise = cat(2,A_noise,units(remove).A_tau);
% %         S = S(~remove,:);
% %         units = units(~remove);
% %         A_tau = A_tau(:,~remove,:);
% %         clear remove
% %     end
% 
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % Collecting final results
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %     if size(A_tau,2) ~= size(S,1)
% %         error('Something went wrong!');
% %     end
% % 
% %     %This time with params.interactive thresholding:
% %     if params.interactive
% %         sign_lev = input('Specify significance level for unimodality:');
% %     end
% %     
% %     [units] = SpikeTimeIdentificationHartigan(S, sr,sign_lev,params.plotting,params.interactive);
% % 
% %     data_tmp = reshape(X,...
% %         [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);
% %     for k = 1:length(units)
% %         units(k).A_tau = A_tau(:,k,:);
% %         %Consider to calculate STAs only based on "non-coincident" spikes!
% %         units(k).STA = GetSTA(data_tmp, units(k).time, sr, 0);
% %         [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
% %             == max(max(max(abs(units(k).STA)))));
% %         units(k).boss_row = sensor_rows_roi(row_max);
% %         units(k).boss_col = sensor_cols_roi(col_max);
% %     end
% %     clear data_tmp
% 
% %     if params.interactive
% %         fprintf('Another chance to inspect final results before they are saved.\n');
% %         keyboard;
% %     end
%     
%     ROIs(i).A = A;
%     ROIs(i).A_dupl = A_dupl;
%     ROIs(i).A_noise = A_noise;
%     ROIs(i).A_mix = A_mix;
%     ROIs(i).A_tau = A_tau;
%     ROIs(i).S = S;
%     ROIs(i).S_dupl = S_dupl;
%     ROIs(i).S_noise = S_noise;
%     ROIs(i).S_mix = S_mix;
%     ROIs(i).W = W;
%     ROIs(i).duplicate_pairs = duplicate_pairs;
%     ROIs(i).units = units;
%     ROIs(i).units_dupl = units_dupl;
%     ROIs(i).units_mix = units_mix;
%     
%     clear A A_dupl A_noise A_tau S S_dupl S_noise W duplicate_pairs units units_dupl
%     
%     clear A_mix S_mix units_mix
%     
%     clear X sensor_rows_roi sensor_cols_roi T_mask N_mask
%     
% 
% end %end loop over regions of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine results of different ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(ROIs) > 1
    fprintf('\nCombining results of different regions of interest...\n');
    t1 = clock;
    [ROIs, N_INTER_ROI_DUPL] = combinerois(ROIs, OL, params.sr, data,...
        params.sensor_rows, params.sensor_cols, params.t_s, params.t_jitter,...
        params.coin_thr, params.sim_thr, params.plotting, params.interactive);
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







