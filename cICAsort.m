function [ X_ROI,sensor_rows_ROI, sensor_cols_ROI, units, S_cica ] = ...
                       cICAsort( data , sensor_rows, sensor_cols, filename)
%cICAsort perform spike sorting of high density array data
%based on convolutive ICA

% created by Christian Leibig 12.02.13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dataset specs:
sr = 11.49; %sampling rate in kHz
d_row = 7.4; %neighbouring sensor distance in \mum
d_col = 7.4;

%ROI identification:
thr_factor = 1;

%fastICA:
nonlinearity = 'pow3';
approach_fastICA = 'defl';

%convolutive ICA:
L = 7;
M = 0;
plotting = 1;
min_skewness = 0.2;
d_max = 1000; %maximal distance in \mum for extrema of component filters
min_corr = 0.2;
approach = 'cluster';
max_cluster_size = 3;
max_iter = 5;
min_no_peaks = 2;
maxlags = 10;

t_total_1 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ X_ROI, sensor_rows_ROI, sensor_cols_ROI ] = ROIIdentification(data,...
                                      sensor_rows, sensor_cols, thr_factor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing with fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = clock;
[S_ica, A, W] = fastica(X_ROI,'g',nonlinearity,...
                        'approach',approach_fastICA);
t2 = clock;
fprintf('fastICA step performed in %g seconds\n',etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutive ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize lagged filters:
A_tau = A;
A_tau(:,:,2:L+1) = 0;
%Perform convolutive ICA:
t1 = clock;
[S_cica, A_tau] = ConvolutiveICA(S_ica,L,A_tau,sr,d_row,d_col,...
     length(sensor_rows_ROI),length(sensor_cols_ROI),d_max,...
    'M',M,'maxlags',maxlags,...
    'plotting',plotting,'min_skewness',min_skewness,'min_corr',min_corr,...
    'approach',approach,'max_cluster_size',max_cluster_size,...
    'max_iter',max_iter,'min_no_peaks',min_no_peaks);
t2 = clock;
fprintf('convolutive ICA step performed in %g seconds\n',etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike time identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = clock;

[units] = SpikeTimeIdentification(S_cica, sr);

t2 = clock;
fprintf('Spike time identification and clustering with KlustaKwik\n');
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
    
SaveResults(filename, units);   
    
t_total_2 = clock;
fprintf('Total cICAsort performed in %g seconds\n',etime(t_total_2,t_total_1));








