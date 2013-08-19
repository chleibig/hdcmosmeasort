function [units, SDscore] = SpikeTimeIdentificationKlustaKwik(...
                                    X,nEig,upsample, sr, show)
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
%correct skewness:
X(skewness(X') > 0,:) = -1 * X(skewness(X') > 0,:);

if show;
    fig1 = figure; fig2 = figure; fig3 = figure; 
    splt_size = ceil(sqrt(dims));
end
units = struct('time',{},'amplitude',{});
SDscore = zeros(1,dims);
%lenRes = round(0.5 * sr*upsample) + round(0.5 * sr*upsample) + 1;
%fullSDs = zeros(lenRes, dims);
for i=1:dims;
    x = resample(X(i,:),upsample,1);
    N_SAMPLES = length(x);
    if(show);figure(fig1);subplot(splt_size,splt_size,i);plot(x);hold on;end
    noise_std = median(abs(x)/0.6745);
    [indices,pos_amplitudes] = find_peaks(-x,5*noise_std,ceil(sr*upsample));
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
        class_score = zeros(max(KluRes.dataClass),1);
        for k = 1:max(KluRes.dataClass)
            %consider different rating scheme
            %class_score(k) = sum(x(indices(KluRes.dataClass==k)));
            if nnz(KluRes.dataClass == k) <= 2
                class_score(k) = 0;
            else
                class_score(k) = norm(mean(pks(:,KluRes.dataClass == k),2));
            end
        end
        %[value_winner,k_winner] = min(class_score);
        [value_winner,k_winner] = max(class_score);
        %SD test:
        residuals_std = std(pks(:,KluRes.dataClass==k_winner),0,2);
        %fullSDs(:,i) = residuals_std - noise_std;
        SDscore(i) = max(residuals_std)/noise_std;
        if show
            colors = hsv(max(KluRes.dataClass));
            figure(fig1);
            for k=1:max(KluRes.dataClass);
                if k == k_winner; marker = 'x';else marker = 'o';end
                plot(indices(KluRes.dataClass==k),...
                    x(indices(KluRes.dataClass==k)),...
                    'Color',colors(k,:),'Marker',marker,'LineStyle','none');
            end
            figure(fig2);subplot(splt_size,splt_size,i);hold on;
            for k=1:max(KluRes.dataClass);
                if k == k_winner;
                    color = 'black';
                else
                    color = [0.75 0.75 0.75];
                end
                plot(pks(:,KluRes.dataClass==k),'Color',color);
            end
            figure(fig3);subplot(splt_size,splt_size,i);hold on;
            plot(residuals_std);
            title(strcat('SDscore = ',num2str(SDscore(i)),...
                '; dip = ',num2str(min(residuals_std)/noise_std)));
            plot(1:f_tot,noise_std*ones(1,f_tot),'-.');
        end
        units(i).time = indices(KluRes.dataClass == k_winner)/(sr*upsample);
        units(i).amplitude = x(indices(KluRes.dataClass==k_winner));
    else
        units(i).time = [];
        units(i).amplitude = [];
        
    end
end


end