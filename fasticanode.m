function [S_ica, A, W, params] = fasticanode(X,params)
%fasticanode encapsulates all the fastica stuff

% Input
% =====
%
% X - (M channels x T samples) array of mixed signals
% params - struct containing all the parameters in the following fields:
%     params.allframes - if true, all frames are used
%     params.allchannels - if true, all channels are used
%     params.frames - boolean array of length T, indicating the frames to
%     be used (overwritten, if allframes is true)
%     params.channels - boolean array of length M, indicating the channels
%     to be used (overwritten, if allchannels is true)
%     params.nonlinearity - the nonlinearity used by fastica
%     params.estimate - if true, the number of extractable components is
%     estimated
%     params.numOfIC - number of components to be extracted (overwritten,
%     if params.estimate is true)
%     params.per_var - keep that many dimensions such that per_var of the total 
              %variance gets explained
%     params.approach - either "symm" or "defl"
%     params.renorm - if true, renormalize W and S such that noise instead
%     of signal is of variance 1 (ref. Jaeckel, 2012)
%
% Output
% ======
%
% S_ica - source activations
% A - mixing matrix
% W - unmixing matrix
% params - output because some parameters might be overwritten
%
% christian.leibig@g-node.org, 16.07.13

noise_frames = ~params.frames;

if params.allframes; params.frames = true(size(X,2),1); end

if params.allchannels; params.channels = true(size(X,1),1); end

if strcmp(params.estimate,'icaDeflation')
    %estimate params.numOfIC
    t1 = clock;
    fprintf('Estimating number of instantaneous components...\n');
    [A, W] = fastica(X(params.channels,params.frames),...
        'g', params.nonlinearity,'approach','defl','verbose',params.verbose);
    
    %0.8 is empirical, because 1.0 was found to exhibit convergence problems
    params.numOfIC = round(0.8 * size(W,1));
    
    clear A W
    t2 = clock;
    fprintf('...done in %g seconds.\n',etime(t2,t1));
end

if strcmp(params.estimate,'svdSpectrum')
    s = svd(X(params.channels, params.frames));
    params.numOfIC = estimatenumberofneurons(s,'median',nnz(params.frames));
end

%compute spectral decomposition of covariance matrix separately and keep
%only that many dimensions such that params.per_var of the variance
%gets explained:
t1 = clock;
[pcaE, pcaD] = fastica(X(params.channels, params.frames),'only','pca',...
    'verbose',params.verbose);
t2 = clock;
fprintf('PCA performed in %g seconds.\n',etime(t2,t1));

%reduce dimensionality - consider reducing dimensionality directly to 
%numOfIC to save runtime
[d,ind] = sort(diag(pcaD),1,'descend');
eigs_to_keep = find(cumsum(d)/sum(d) <= params.per_var);
nEig = eigs_to_keep(end);
pcaE = pcaE(:,ind(1:nEig));
pcaD = pcaD(ind(1:nEig),ind(1:nEig));
fprintf('%g out of %g eigenvectors are kept.\n',nEig,nnz(params.channels));

if strcmp(params.estimate,'eigSpectrum')
   params.numOfIC = estimatenumberofneurons(d,'median',nnz(params.frames)); 
end

%fastica using pca results
fprintf('Trying to extract %g independent components.\n',params.numOfIC);
t1 = clock;
[A, W] = fastica(X(params.channels, params.frames),...
    'g',params.nonlinearity,'approach',params.approach,...
    'numOfIC',params.numOfIC,'pcaE', pcaE,...
    'pcaD', pcaD,'verbose',params.verbose);

if params.renorm
    %Normalize DCVs such that noise in S = W*X has unit variance
    C_noise = X(params.channels,noise_frames)*...
        X(params.channels,noise_frames)'/nnz(noise_frames);
    W_renorm = W;
    for k = 1:size(W,1)
        W_renorm(k,:) = W(k,:)/sqrt(W(k,:)*C_noise*W(k,:)');
    end
    W = W_renorm;
end
%Compute the source activations over all samples T
S_ica = W*X(params.channels,:);
t2 = clock;
fprintf('Extraction of ICs performed in %g seconds\n',etime(t2,t1));


end