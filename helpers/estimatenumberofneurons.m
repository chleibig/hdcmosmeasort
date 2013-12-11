function N = estimatenumberofneurons(spectrum)
% N = estimatenumberofneurons(spectrum) thresholds
% spectrum (either singular values of data matrix
% or eigenvalues of data covariance matrix).
% The idea is to locate the kink in the spectrum.

% christian.leibig@g-node.org, 11.12.2013

N = find(spectrum < 2*median(spectrum) - min(spectrum),1);


end