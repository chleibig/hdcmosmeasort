function [ ROI ] = processroi( ROI, params )
%[ ROI ] = processroi( ROI, params ) process - e.g. sort - a single region 
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
% christian.leibig@g-node.org, 10.09.13
%
  
fprintf('\nWorking on ROI %g...\n\n',ROI.k);

X = ROI.X;
sensor_rows_roi = ROI.sensor_rows;
sensor_cols_roi = ROI.sensor_cols;
T_mask = ROI.T_mask;
N_mask = ROI.N_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing with fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.ica.frames = T_mask;
params.ica.channels = N_mask;

params.ica.numOfIC = ceil(params.ica.cpn/(params.sensor_rho/params.neuron_rho) * ...
    nnz(params.ica.channels));

[S, A, W, params.ica] = fasticanode(X,params.ica);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutive ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if params.do_cICA
    %Initialize lagged filters:
    A_tau = zeros(size(X,1),size(A,2));
    A_tau(N_mask,:) = A;
    A_tau(:,:,2:params.L+1) = 0;

    if params.allframes_cica
        frames_ROI_cica = true(size(X,2),1);
    else
        frames_ROI_cica = T_mask;
    end

    %Perform convolutive ICA:
    t1 = clock;
    [S, A_tau, S_noise, A_noise] = ConvolutiveICA(S,params.L,A_tau,params.sr,...
        params.d_row,params.d_col,length(sensor_rows_roi),length(sensor_cols_roi),params.d_max,...
        frames_ROI_cica,params.do_cICA,'M',params.M,'maxlags',params.maxlags,...
        'plotting',params.plotting,'min_skewness',params.min_skewness,'min_corr',params.min_corr,...
        'approach',params.grouping,'max_cluster_size',params.max_cluster_size,...
        'max_iter',params.max_iter,'thrFactor',params.thrFactor,...
        'min_no_peaks',params.min_no_peaks,...
        't_s',params.t_s,'t_jitter',params.t_jitter, 'coin_thr',params.coin_thr);
    t2 = clock;
    fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));
else
    fprintf('Convolutive ICA is not applied!\n');
end

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
[units] = SpikeTimeIdentificationKlustaKwik(S,0,10, params.sr,params.thrFactor,params.plotting);
t2 = clock;
fprintf('performed in %g seconds\n',etime(t2,t1));

%     fprintf('Spike time identification with Hartigans dip test\n');
%     t1 = clock;
%     [units] = SpikeTimeIdentificationHartigan(S, params.sr,params.sign_lev,params.plotting,1);
%     t2 = clock;
%     fprintf('performed in %g seconds\n',etime(t2,t1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove mixed units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    %Consider to calculate STAs only based on "non-coincident" spikes!
    units(k).STA = GetSTA(data_tmp, units(k).time, params.sr, 0);
    [row_max,col_max] = find(max(abs(units(k).STA),[],3)...
        == max(max(max(abs(units(k).STA)))));
    units(k).boss_row = sensor_rows_roi(row_max);
    units(k).boss_col = sensor_cols_roi(col_max);
end
clear data_tmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove duplicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = clock;
fprintf('Checking for duplicates...\n');
[duplicate_pairs] = checkforintraroiduplicates(units, params.sr, ...
    params.t_s, params.t_jitter, params.coin_thr, params.sim_thr, params.plotting, params.interactive);
N_dupl = size(duplicate_pairs,1);
t2 = clock;
fprintf('found %g intraregional duplicates in %g seconds\n',...
    N_dupl,etime(t2,t1));

%dbstop in cICAsort.m at 326 if (N_dupl > 0)
%Experiment with additional criteria to decide upon which duplicate
%partner to remove:
for d = 1:N_dupl
    %         SpikeTimeIdentificationKlustaKwik(S(duplicate_pairs(d,:),:),0,10, sr, 1);
    %         if (units(duplicate_pairs(d,1)).RSTD > 1.5*units(duplicate_pairs(d,2)).RSTD)
    %             %duplicate_pairs(d,1) is considered to be a mixture and will be
    %             %removed
    %             break;
    %         end
    %         if (units(duplicate_pairs(d,2)).RSTD > 1.5*units(duplicate_pairs(d,1)).RSTD)
    %             %duplicate_pairs(d,1) is considered to be a mixture and will be
    %             %removed
    %             duplicate_pairs(d,:) = duplicate_pairs(d,end:-1:1);
    %             break;
    %         end

    %No mixture detected - the unit with higher separability will
    %be kept
    if units(duplicate_pairs(d,1)).separability <= units(duplicate_pairs(d,2)).separability
        %remove the first
    else
        %remove the second
        duplicate_pairs(d,:) = duplicate_pairs(d,end:-1:1);
    end
end



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
