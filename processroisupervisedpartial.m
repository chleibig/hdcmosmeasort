function [ ROI ] = processroisupervisedpartial( ROI, params)
%[ ROI ] = processroisupervisedpartial( ROI, params ) process - e.g. sort - a single region 
%of interest. This function acts as a wrapper to be able to work on 
%different ROIs in parallel. It implements a supervised sanity checked
%based on optimal filtering. Only one neuron per ROI is accepted.
%This allows to test whether optimal filtering can be deployed if not all
%templates were found.
%
% Input
% =====
%
% ROI - struct containing the information of a single region of interest
% params - struct containing all necessary parameters, in particular
% params.gtFilename
%
%
% Output
% ======
%
% ROI - input ROI get's changed in place
%
%
% christian.leibig@g-node.org, 23.11.14
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

X = readDataBlock(params.filename,...
                    firstRow,lastRow,firstCol,lastCol,firstFrame,lastFrame);

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%masks to index into the data block
T_mask = ROI.T_mask;
N_mask = ROI.N_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain mixing matrices by using true spike train information for each
% neuron and take the optimal filter response as source activation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Get true spike trains of neurons within ROI %%%%%%%%%%%%%%%%%%%%%%%%%
[neuron, time, sensor_position_col, sensor_position_row] = textread(...
    params.gtFilename,...
    '%*s %s %f %*f %u %u %*f','delimiter','\t','headerlines',1);

%%%%% Figure out which neurons are located within ROI %%%%%%%%%%%%%%%%%%%%%
allNeuronLabels = unique(neuron);
neuronLabels = {};%Neurons to be accepted
for i = 1:length(allNeuronLabels)
    sensor_col_i = sensor_position_col(strcmp(neuron,allNeuronLabels(i)));
    sensor_col_i = sensor_col_i(1);
    sensor_row_i = sensor_position_row(strcmp(neuron,allNeuronLabels(i)));
    sensor_row_i = sensor_row_i(1);
    if (sensor_col_i <= sensor_cols_roi(end) && ...
            sensor_col_i >= sensor_cols_roi(1) && ...
            sensor_row_i <= sensor_rows_roi(end) && ...
            sensor_row_i >= sensor_rows_roi(1))
        %Neuron i is within rectangle around ROI
        [data_row_i, data_col_i] = sensorcoord2datacoord(...
            sensor_row_i,sensor_col_i,...
            sensor_rows_roi, sensor_cols_roi);
        linear_idx = sub2ind(...
            [length(sensor_rows_roi) length(sensor_cols_roi)],...
            data_row_i, data_col_i);
        if N_mask(linear_idx)
            %Neuron is in ROI
            neuronLabels{end+1} = allNeuronLabels(i);
        end
    end
end

clear sensor_position_col sensor_position_row

N_NEURONS = length(neuronLabels);

fprintf('Current ROI contains %g true neurons, but only 1 is accepted.\n',...
    N_NEURONS);
N_NEURONS = 1; %random selection of one neuron only.

%%%%%%%%%%%%%%%% Collect templates (mixing matrices) %%%%%%%%%%%%%%%%%%%%%%

Tf = round(0.5 * params.sr) + round(0.5 * params.sr) + 1;
templates = zeros(nnz(N_mask)*Tf, N_NEURONS);
A_tau = zeros(length(sensor_rows_roi)*length(sensor_cols_roi),N_NEURONS,Tf);

for i = 1:N_NEURONS
    tspk = time(strcmp(neuron,neuronLabels{i}));
    %Calculate STA on rectangular box surrounding the ROI
    [template] = GetSTA(X,tspk,params.sr,0);
    template = reshape(template,...
        [size(template,1)*size(template,2) size(template,3)]);
    %use only values inside ROI for optimal filtering:
    template = template(N_mask,:);
    A_tau(N_mask,i,:) = template;

    [Ntemplate,tmp] = size(template);
    assert(Ntemplate == nnz(N_mask) & Tf == tmp);
    templates(:,i) = reshape(template',[Ntemplate*Tf 1]);
end



%%%%%%%%%%%%%% Template space embedding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = reshape(X,[size(X,1)*size(X,2) size(X,3)]);
T = length(params.frameStartTimes);
N = nnz(N_mask);
Xbar = zeros(N*Tf,T - mod(T,Tf));
for t = 1:T-Tf
    Xbar(:,t) = reshape(X(N_mask,t:t+Tf-1)',[N*Tf 1]);
end

%%%%%%%%%%%%%% Data covariance matrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For notational simplicity evaluate covariance matrix in the
% computationally inefficient way (equation 3.1, F. Franke Diss)
tic;
R = Xbar*Xbar'/(T - Tf + 1);
toc;
invR = inv(R);
%%%%%%%%%%%%%% Optimal filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = invR*templates/(templates'*invR*templates);


S = F'*Xbar;

clear Xbar R templates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct sign such that all spikes are negative deflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = -S;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike time identificatin accounts for crosstalk from other neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Spike time identification and clustering with KlustaKwik\n');
t1 = clock;
[units] = SpikeTimeIdentificationKlustaKwik(S,0,params.upsample, params.sr,params.thrFactor,params.plotting);
t2 = clock;
fprintf('performed in %g seconds\n',etime(t2,t1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting preliminary results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_tmp = reshape(X,...
    [length(sensor_rows_roi) length(sensor_cols_roi) size(X,2)]);
clear X

fprintf('Computing skewness...\n');
skewn = skewness(S');
fprintf('Computing kurtosis...\n');
kurt = kurtosis(S');

for k = 1:length(units)
    units(k).A_tau = A_tau(:,k,:);
    %Consider to calculate STAs only based on "non-coincident" spikes!
    units(k).STA = GetSTA(data_tmp, units(k).time, params.sr, 0);
    extrSTA = max(max(max(abs(units(k).STA))));
    [row_max,col_max] = find(max(abs(units(k).STA),[],3) == extrSTA);
    units(k).boss_row = sensor_rows_roi(row_max);
    units(k).boss_col = sensor_cols_roi(col_max);
    units(k).snr = extrSTA/median(abs(data_tmp(row_max,col_max,:))/0.6745);
    units(k).skewn = skewn(k);
    units(k).kurt = kurt(k);
end
clear data_tmp skewn kurtosis



ROI.A_tau = A_tau;
ROI.S = S;
ROI.units = units;

clear A_tau S units


end