function [keep] = checkfornoisycomponents(X,minSkew,minNoPeaks,sr,show)
%checkfornoisycomponents tests row vectors of X for
%both a minimum skewness and minimum number of peaks
%
% Input
% =====
%
% X - (N_COMPONENTS x N_SAMPLES) array
% minSkew - minimum skewness for components to be kept
% minNoPeaks - minimum number of peaks for components to be kept
% sr - sampling rate in kHz
% show - if true, plots are shown
% 
%
% Output
% ======
%
% keep - boolean vector, positions of components that are at least minSkew
%        skewed and show at least minNoPeaks peaks are assigned true 
%
% christian.leibig@g-node.org, 16.07.13
%

% 1. set mask based on skewness and presence of peaks criteria
skewn = skewness(X');
%assess number of peaks in each component:
n_peaks = zeros(1,size(X,1));

for i = 1:size(X,1)
    [indices, peaks] = find_peaks(abs(X(i,:)),...
        5*median(abs(X(i,:))/0.6745),ceil(sr));
    n_peaks(i) = length(peaks);
end

keep = (abs(skewn) > minSkew) & (n_peaks >= minNoPeaks);

fprintf('%g from the %g input components are considered noise...\n',...
    length(nonzeros(~keep)),size(X,1));

if show
    figure;
    subplot(2,1,1);hist(skewn,floor(sqrt(length(skewn))));title('Skewness');
    xlabel('skewness');ylabel('counts');
    subplot(2,1,2);hist(n_peaks,floor(sqrt(length(n_peaks))));title('Number of peaks per component');
    xlabel('Number of peaks per component');ylabel('counts');
end


end
