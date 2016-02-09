function [indices, peaks] = searchpeaks(X,thr,L)
% [indices, peaks] = searchpeaks(X,thr,L)
% Find local maxima in X (vector) that are 
% above threshold thr and within L samples after
% threshold crossing

%Author: Christian Leibig, 29.11.12

onset_ind = find(diff(X > thr) > 0);
if isempty(onset_ind); indices = [];peaks = []; return;
else
    indices = zeros(1,length(onset_ind));
    peaks = zeros(1,length(onset_ind));
    X(end:end+L) = 0;%zero padding to allow peaks to occur at the very end
    for i=1:length(peaks)
        [peaks(i) ind] = max(X(onset_ind(i):onset_ind(i)+L));
        indices(i) = ind + onset_ind(i) - 1;
    end
end

end
    
    
    

