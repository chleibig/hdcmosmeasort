function [ ROI ] = processroicicacicaar( ROI, params )
%[ ROI ] = processroicicacicaar( ROI, params ) process - e.g. sort - 
%a single region of interest. This function acts as a wrapper to be able to
%work on different ROIs in parallel.
%
% Input
% =====
%
% ROI - struct containing the information of a single region of interest
% params - struct containing all necessary parameters.
%
%
% Output
% ======
%
% ROI - input ROI get's changed in place
%
%
% christian.leibig@g-node.org
%
  
fprintf('\nWorking on ROI %g...\n\n',ROI.k);

d_max = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load box shaped data block that contains region of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sensor_rows_roi = ROI.sensor_rows;
sensor_cols_roi = ROI.sensor_cols;
      
firstRow = find(params.sensor_rows == sensor_rows_roi(1));
lastRow = find(params.sensor_rows == sensor_rows_roi(end));
firstCol = find(params.sensor_cols == sensor_cols_roi(1));
lastCol = find(params.sensor_cols == sensor_cols_roi(end));
firstFrame = 1;
lastFrame = length(params.frameStartTimes);

X = readdatablock(params.filename,...
                    firstRow,lastRow,firstCol,lastCol,firstFrame,lastFrame);
X = reshape(X,[size(X,1)*size(X,2) size(X,3)]);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%masks to index into the data block
T_mask = ROI.T_mask;
N_mask = ROI.N_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing with fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.ica.frames = T_mask;
params.ica.channels = N_mask;

[S, A, W, params.ica] = fasticanode(X,params.ica);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutive ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize lagged filters:
A_tau = zeros(size(X,1),size(A,2));
A_tau(N_mask,:) = A;
A_tau(:,:,2:(params.L+params.M+1)) = 0;

if params.allframes_cica
    frames_ROI_cica = true(size(X,2),1);
else
    frames_ROI_cica = T_mask;
end

%Perform convolutive ICA:
t1 = clock;
[S, A_tau, S_noise, A_noise] = convolutiveica(S,params.L,A_tau,params.sr,...
    params.d_row,params.d_col,length(sensor_rows_roi),length(sensor_cols_roi),d_max,...
    frames_ROI_cica,'M',params.M,'maxlags',params.maxlags,...
    'plotting',params.plotting,'min_skewness',params.min_skewness,'min_corr',params.min_corr,...
    'max_cluster_size',params.max_cluster_size,...
    'max_iter',params.max_iter,'thrFactor',params.thrFactor,...
    'min_no_peaks',params.min_no_peaks,...
    't_s',params.t_s,'t_jitter',params.t_jitter, 'coin_thr',params.coin_thr);
t2 = clock;
fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct skewness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Correct skewness such that all spikes are negative deflections.
skewn = skewness(S');
S(skewn > 0,:) = -1 * S(skewn > 0,:);
clear skewn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for noisy components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if convolutive ICA was not applied, we need to initialize some variables:
if ~exist('A_noise','var'); A_noise = []; end
if ~exist('S_noise','var'); S_noise = []; end
if ~exist('A_tau','var')
    A_tau = zeros(size(X,1),size(A,2));
    A_tau(N_mask,:) = A;
    A_tau(:,:,2:params.L+1) = 0;
end

[keep] = checkfornoisycomponents(S,params.min_skewness,params.thrFactor,...
                            params.min_no_peaks,params.sr,params.plotting);

% store noisy stuff away and remove it from components and filters:
S_noise = [S_noise;S(~keep,:)];
A_noise = cat(2,A_noise,A_tau(:,~keep,:));
S = S(keep,:);
A_tau = A_tau(:,keep,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike time identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Spike time identification and clustering with KlustaKwik\n');
t1 = clock;
[units] = spiketimeidentificationklustakwik(S,0,params.upsample, ...
                                            params.sr,params.thrFactor,...
                                            params.plotting);
t2 = clock;
fprintf('performed in %g seconds\n',etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove mixed units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fuse these with noise and include all other criteria from chapter 3 of
%Diss
if ~exist('A_mix','var'); A_mix = []; end
if ~exist('S_mix','var'); S_mix = []; end
if ~exist('units_mix','var'); units_mix = []; end

if length(units) > 0

    clear keep
    keep = ([units.RSTD] <= params.maxRSTD);
    %dbstop in cICAsort.m at 276 if (nnz(~keep) > 0)


    S_mix = [S_mix;S(~keep,:)];
    A_mix = cat(2,A_mix,A_tau(:,~keep,:));
    units_mix = units(~keep);

    S = S(keep,:);
    A_tau = A_tau(:,keep,:);
    units = units(keep);

    fprintf('Removed %g units supposed to contain mixtures.\n',nnz(~keep));

    clear keep

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting preliminary results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_tmp = reshape(X,...
    [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);

fprintf('Computing skewness...\n');
skewn = skewness(S');
fprintf('Computing kurtosis...\n');
kurt = kurtosis(S');

for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    %Consider to calculate STAs only based on "non-coincident" spikes!
    units(k).STA = computetemplate(data_tmp, units(k).time, params.sr, 0);
    extrSTA = max(max(max(abs(units(k).STA))));
    [row_max,col_max] = find(max(abs(units(k).STA),[],3) == extrSTA);
    units(k).boss_row = sensor_rows_roi(row_max);
    units(k).boss_col = sensor_cols_roi(col_max);
    units(k).snr = extrSTA/median(abs(data_tmp(row_max,col_max,:))/0.6745);
    units(k).skewn = skewn(k);
    units(k).kurt = kurt(k);
end
clear data_tmp skewn kurtosis


units_dupl = [];
S_dupl = [];
A_dupl = [];
duplicate_pairs = [];


ROI.A = A;
ROI.A_dupl = A_dupl;
ROI.A_noise = A_noise;
ROI.A_mix = A_mix;
ROI.A_tau = A_tau;
ROI.S = S;
ROI.S_dupl = S_dupl;
ROI.S_noise = S_noise;
ROI.S_mix = S_mix;
ROI.W = W;
ROI.duplicate_pairs = duplicate_pairs;
ROI.units = units;
ROI.units_dupl = units_dupl;
ROI.units_mix = units_mix;

clear A A_dupl A_noise A_tau S S_dupl S_noise W duplicate_pairs units units_dupl

clear A_mix S_mix units_mix

clear X sensor_rows_roi sensor_cols_roi T_mask N_mask




end
