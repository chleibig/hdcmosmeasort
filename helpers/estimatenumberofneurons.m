function N = estimatenumberofneurons(spectrum,approach,T)
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
    case 'mean'
        N = find(spectrum <= mean(spectrum),1);
    case 'AIC'
        n = length(spectrum);
        rho = arrayfun(@(m) (prod(spectrum(m+1:n).^(1/(n-m)))/(sum(spectrum(m+1:n))/(n-m))), 1:n);
        AIC = arrayfun(@(m) (-2*T*(n-m)*log(rho(m))+2*m*(2*n-m)),1:n);
        [unused,N] = min(AIC);
%         figure;plotyy(1:n,spectrum,1:n,AIC);
%         title(['minimum at ' num2str(N)]);
        
    case 'MDL'
        n = length(spectrum);
        rho = arrayfun(@(m) (prod(spectrum(m+1:n).^(1/(n-m)))/(sum(spectrum(m+1:n))/(n-m))), 1:n);
        MDL = arrayfun(@(m) (-T*(n-m)*log(rho(m))+0.5*m*(2*n-m)*log(T)),1:n);
        [unused,N] = min(MDL);
%         figure;plotyy(1:n,spectrum,1:n,MDL);
%         title(['minimum at ' num2str(N)]);
    otherwise
        error('Please specify a valid estimation mode!')
end

fprintf(['%g neurons estimated with ' approach '\n'],N);

end