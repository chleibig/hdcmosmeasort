function [units] = SpikeTimeIdentificationKlustaKwik(X,nEig,upsample, sr,...
                                                     thrFactor, show)
% [units] = SpikeTimeIdentificationKlustaKwik(X,nEig,upsample,sr,thrFactor,show)
% 
% Perform spike time identification on each component
% of X (dims x samples) and disambiguate the most prominent neuron
% per component with KlustaKwik
%
% Input
% =====
%
% X - (N_DIMS x N_SAMPLES) array
%
% nEig - if = 0 <-> amplitudes are used as feature, otherwise PCA is
%        performed on aligned waveforms and scores of nEig most prominent
%        eigenvectors are used as features for KlustaKwik
%
% upsample - factor by which signals get upsampled
%            (before threshold detection)
%
% sr - sampling rate of X in kHz
%
% thrFactor - multiple of noise std to be used as peak threshold
%
% show - flag to control graphical output
%
%
% Output
% ======
%
% units - struct array with the fields .time and .amplitude 
%
% SDscore
%
% contact: christian.leibig@g-node.org
%
%

[dims, samples] = size(X);

if show;
    fig1 = figure; fig2 = figure; %fig3 = figure; 
    splt_size = ceil(sqrt(dims));
end
units = struct('time',{},'amplitude',{});

for i=1:dims;
    x = resample(X(i,:),upsample,1);
    N_SAMPLES = length(x);
    if(show);figure(fig1);subplot(splt_size,splt_size,i);plot(x);hold on;end
    noise_std = median(abs(x)/0.6745);
    %Skewness is assumed to be corrected previously such that spikes
    %are negative deflections.
    [indices,pos_amplitudes] = find_peaks(-x,thrFactor*noise_std,ceil(sr*upsample));
    amplitudes = -1*pos_amplitudes;
    if ~isempty(amplitudes)
        %Get waveforms.
        pre = round(0.5 * sr*upsample);
        post = round(0.5 * sr*upsample);
        f_tot = pre + post + 1;
        %take only those peaks for which the desired window is entirely contained
        %in the data:
        amplitudes = amplitudes( ((indices - pre) >= 1) & ((indices + post) <= N_SAMPLES) );
        indices = indices( ((indices - pre) >= 1) & ((indices + post) <= N_SAMPLES) );
        N_pks = length(indices);
        pks = zeros(f_tot,N_pks);
        for j = 1:N_pks
            pks(:,j) = x(indices(j)-pre:indices(j)+post);
        end
        
        switch nEig
            case 0
                %Perform KlustaKwik on amplitudes.
                KluRes = doKlustaKwik(amplitudes');               
            otherwise
                %Perform PCA on waveforms.
                Cpks = pks*pks'/(size(pks,2) - 1);
                [V,D] = eig(Cpks);
                pcasig = V(:,end-nEig+1:end)'*pks;
                %Perform KlustaKwik on nEig most prominent pca scores.
                KluRes = doKlustaKwik(pcasig');
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform potentially a merging
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        K = max(KluRes.dataClass);
        kRange = unique(KluRes.dataClass);
        N_CLU = length(kRange);
        
        clear k_winner k_nearest
        
        if N_CLU == 1
            %there is only one cluster
            k_winner = K;
        end
        if N_CLU > 1
            class_score = zeros(K,1);
            for k = kRange(1):kRange(end)
                class_score(k) = norm(mean(pks(:,KluRes.dataClass == k),2));
            end
            [unused, ind] = sort(class_score, 1, 'descend');
            clear class_score kRange
            %merge potential outliers:
            if  (nnz(KluRes.dataClass == ind(1)) == 1) || ...
                ( ( (std(amplitudes(KluRes.dataClass == ind(1))) <= 0.2) || ...
                    (std(amplitudes(KluRes.dataClass == ind(2))) <= 0.2) ) ...
                && ( abs(mean(abs(amplitudes(KluRes.dataClass == ind(1)))) ...
                    - mean(abs(amplitudes(KluRes.dataClass == ind(2))))) <= 2.5 )) 
                
                %merge second largest unit with largest unit.
                tomerge = ind(2);
                KluRes.dataClass(KluRes.dataClass == tomerge) = ind(1);
                %adopt k's.
                KluRes.dataClass(KluRes.dataClass > tomerge) = ...
                    KluRes.dataClass(KluRes.dataClass > tomerge) - 1;
                ind(2) = ind(1);
                ind = ind - (ind > tomerge);
                
                k_winner = ind(1);
                try
                    k_nearest = ind(3);
                catch
                    N_CLU = 1;
                end
            else
                k_winner = ind(1);
                k_nearest = ind(2);
            end
            clear ind
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fill units
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        
        %SD test:
        residuals_std = std(pks(:,KluRes.dataClass==k_winner),0,2);
        units(i).SDscore = max(residuals_std);
        units(i).amplitudeSD = std(amplitudes(KluRes.dataClass==k_winner));
        meanAmpWinner = mean(abs(amplitudes(KluRes.dataClass == k_winner)));
        units(i).RSTD = units(i).amplitudeSD/meanAmpWinner;
        
        if N_CLU > 1
            units(i).separability = ( meanAmpWinner - ...
                mean(abs(amplitudes(KluRes.dataClass == k_nearest))) )/noise_std;
            units(i).IsoINN = isoi(...
                              amplitudes(KluRes.dataClass == k_winner),...
                              amplitudes(KluRes.dataClass == k_nearest));
            units(i).IsoIBg = isoi(...
                              amplitudes(KluRes.dataClass == k_winner),...
                              amplitudes(KluRes.dataClass ~= k_winner));
        else
            units(i).separability = meanAmpWinner/noise_std; %corresponds to average
            %signal to noise ratio (noise is assumed to have unit variance)
            %there is no nearest neighbour or background
            units(i).IsoINN = NaN;
            units(i).IsoIBg = NaN;
        end
        
        units(i).time = indices(KluRes.dataClass == k_winner)/(sr*upsample);
        units(i).amplitude = x(indices(KluRes.dataClass==k_winner));
        units(i).noise_std = noise_std;
        
        
        
        
        if show
            K = max(KluRes.dataClass);
            kRange = unique(KluRes.dataClass);
            N_CLU = length(kRange);
            colors = hsv(K);
            %colors = hsv(max(KluRes.dataClass));
            figure(fig1);
%             for k=1:max(KluRes.dataClass);
            for k=kRange(1):kRange(end);
                if k == k_winner; marker = 'x';else marker = 'o';end
                plot(indices(KluRes.dataClass==k),...
                    x(indices(KluRes.dataClass==k)),...
                    'Color',colors(k,:),'Marker',marker,'LineStyle','none');
            end
            title(strcat('amplSD = ',num2str(units(i).amplitudeSD,2),'; RSTD = ',...
             num2str(units(i).RSTD,2),'; sep. = ',num2str(units(i).separability,2),...
             '; IsoIBg = ',num2str(units(i).IsoIBg,2),...
             '; IsoINN = ',num2str(units(i).IsoINN,2)));
            
            figure(fig2);subplot(splt_size,splt_size,i);hold on;
            %for k=1:max(KluRes.dataClass);
            for k=kRange(1):kRange(end);
                if k == k_winner;
                    color = 'black';
                else
                    color = [0.75 0.75 0.75];
                end
                plot(pks(:,KluRes.dataClass==k),'Color',color);
            end
%             
%             figure(fig3);subplot(splt_size,splt_size,i);hold on;
%             plot(residuals_std);
%             title(strcat('SDscore = ',num2str(units(i).SDscore),...
%                 '; dip = ',num2str(min(residuals_std)/noise_std)));
%             plot(1:f_tot,noise_std*ones(1,f_tot),'-.');
        end
        
    else
        units(i).time = [];
        units(i).amplitude = [];
        units(i).noise_std = noise_std;
    end
    
end


end