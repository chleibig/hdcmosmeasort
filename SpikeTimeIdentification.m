function [units] = SpikeTimeIdentification(X, sr, show)
% Perform spike time identification on each component
% of X (dims x samples) individually via 
% determining threshold crossing events and clustering
% them based on their amplitude with KlustaKwik

[dims, samples] = size(X);
%correct skewness:
X(skewness(X') > 0,:) = -1 * X(skewness(X') > 0,:);

if(show);figure;splt_size = ceil(sqrt(dims));end
units = struct('time',{},'amplitude',{});
for i=1:dims;
    if(show);subplot(splt_size,splt_size,i);plot(X(i,:));hold on;end
    [indices,pos_amplitudes] = find_peaks(-X(i,:),...
                                    5*median(abs(X(i,:))/0.6745),ceil(sr));
    amplitudes = -1*pos_amplitudes;
    if ~isempty(amplitudes)
        KluRes = doKlustaKwik(amplitudes');
        class_score = zeros(max(KluRes.dataClass),1);
        for k = 1:max(KluRes.dataClass)
            class_score(k) = sum(X(i,indices(KluRes.dataClass==k)));
        end
        [value_winner,k_winner] = min(class_score);
        if show
            colors = hsv(max(KluRes.dataClass));
            for k=1:max(KluRes.dataClass);
                if k == k_winner; marker = 'x';else marker = 'o';end
                plot(indices(KluRes.dataClass==k),...
                    X(i,indices(KluRes.dataClass==k)),...
                    'Color',colors(k,:),'Marker',marker,'LineStyle','none');
            end
        end
        units(i).time = indices(KluRes.dataClass == k_winner)/sr;
        units(i).amplitude = X(i,indices(KluRes.dataClass==k_winner));
    else
        units(i).time = [];
        units(i).amplitude = [];
    end
end



end