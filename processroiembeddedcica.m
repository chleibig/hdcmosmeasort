function [ ROI ] = processroiembeddedcica( ROI, params )
%[ ROI ] = processroiembeddedcica( ROI, params ) process - e.g. sort - a single region 
%of interest. This function acts as a wrapper to be able to work on 
%different ROIs in parallel.
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
% christian.leibig@g-node.org, 16.02.15
%
  
fprintf('\nWorking on ROI %g...\n\n',ROI.k);


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
% T_mask - not used, as we do not want to bother with different chunks of
%          spatiotemporally embedded data
N_mask = ROI.N_mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatiotemporal embedding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Xbar] = spatiotemporalembedding(X(N_mask,:), params.L+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate maximal number of neurons from covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = clock;
[pcaE, pcaD] = fastica(Xbar,'only','pca','verbose',params.ica.verbose);
t2 = clock;
fprintf('PCA performed in %g seconds.\n',etime(t2,t1));

[d,ind] = sort(diag(pcaD),1,'descend');

if strcmp(params.ica.estimate,'eigSpectrum')
   nEig = estimatenumberofneurons(d,'linearFit', size(Xbar,2)); 
   params.ica.numOfIC = nEig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce dimensionality to maximal number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcaE = pcaE(:,ind(1:nEig));
pcaD = pcaD(ind(1:nEig),ind(1:nEig));

pcsig = pcaE'*Xbar;
fprintf('Dimensionality reduced to %g out of %g PCs.\n',nEig,...
    size(Xbar,1));
clear Xbar

fprintf('Trying to extract %g independent components.\n',params.ica.numOfIC);
t1 = clock;
[A, W] = fastica(pcsig,'g',params.ica.nonlinearity,...
    'approach',params.ica.approach,'numOfIC',params.ica.numOfIC,...
    'verbose',params.ica.verbose);
t2 = clock;
fprintf('Extraction of ICs performed in %g seconds\n',etime(t2,t1));
a
S = W*pcsig;
clear pcsig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disembed mixing/unmixing matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %calculate effective As and Ws in original sensor space
A = pcaE*A; % A is now nnz(N_mask) x M
W = W*pcaE';

% Get A_tau (N_sensors x M x L + 1) from A
A_tau = zeros(length(sensor_rows_roi)*length(sensor_cols_roi),...
              size(A,2), params.L+1);
for k=1:size(A,2)
      A_tau(N_mask,k,:) = reshape(A(:,k), [params.L+1 nnz(N_mask)])';
end

if params.ica.renorm
    error('Reormalization not implemented, see other processroi examples!');
end

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

if ~exist('A_noise','var'); A_noise = []; end
if ~exist('S_noise','var'); S_noise = []; end

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



ROI.A = A;
ROI.A_noise = A_noise;
ROI.A_mix = A_mix;
ROI.A_tau = A_tau;
ROI.S = S;
ROI.S_noise = S_noise;
ROI.S_mix = S_mix;
ROI.W = W;
ROI.units = units;
ROI.units_mix = units_mix;

clear A A_noise A_tau S S_noise W units

clear A_mix S_mix units_mix

clear X sensor_rows_roi sensor_cols_roi T_mask N_mask




end
