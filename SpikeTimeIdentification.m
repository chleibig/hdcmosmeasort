function [units] = SpikeTimeIdentification(X, sr)
% Perform spike time identification on each component
% of X (dims x samples) individually via 
% determining threshold crossing events and clustering
% them based on their amplitude with KlustaKwik

[dims, samples] = size(X);
%correct skewness:
X(skewness(X') > 0,:) = -1 * X(skewness(X') > 0,:);

figure;
splt_size = ceil(sqrt(dims));
units = struct('time',{},'amplitude',{});
for i=1:dims;
    subplot(splt_size,splt_size,i);plot(X(i,:));hold on;...
    [indices,pos_amplitudes] = find_peaks(-X(i,:),...
                                    5*median(abs(X(i,:))/0.6745),ceil(sr));
    amplitudes = -1*pos_amplitudes;
    KluRes = doKlustaKwik(amplitudes');
    colors = hsv(max(KluRes.dataClass));
    class_score = zeros(max(KluRes.dataClass),1);
    for k = 1:max(KluRes.dataClass)
        class_score(k) = sum(X(i,indices(KluRes.dataClass==k)));
    end
    [value_winner,k_winner] = min(class_score);
    for k=1:max(KluRes.dataClass);
        if k == k_winner; marker = 'x';else marker = 'o';end
        plot(indices(KluRes.dataClass==k),...
             X(i,indices(KluRes.dataClass==k)),...
            'Color',colors(k,:),'Marker',marker,'LineStyle','none');
    end
    units(i).time = indices(KluRes.dataClass == k_winner)/sr;
    units(i).amplitude = X(i,indices(KluRes.dataClass==k_winner));
end



end