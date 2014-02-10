function N = estimatenumberofneurons(spectrum,approach)
% N = estimatenumberofneurons(spectrum) thresholds
% spectrum (either singular values of data matrix
% or eigenvalues of data covariance matrix).
% The idea is to locate the kink in the spectrum.

% christian.leibig@g-node.org, 11.12.2013

spectrum = sort(spectrum,'descend');

switch approach
    case 'median'
        N = find(spectrum < 2*median(spectrum) - min(spectrum),1);
    case 'kde'
        [unused,density,xmesh,unused]=kde(spectrum);
        [unused,i] = max(density);
        mode = xmesh(i);
        if mode <= min(spectrum)
            N = length(spectrum);
        else
            N = find(spectrum < 2*mode - min(spectrum),1);
        end
    case 'linearFit'
        spectrum = spectrum(:)';%guarantee spectrum to be a row vector
        N_DIM = length(spectrum);
        [P,S] = polyfit(round(N_DIM/2):N_DIM,spectrum(round(N_DIM/2):N_DIM),1);
        [Y] = polyval(P,1:N_DIM,S);
        N = find(spectrum < (Y+3*S.normr),1);
    otherwise
        error('Please specify a valid estimation mode!')
end

end