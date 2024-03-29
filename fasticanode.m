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



%compute spectral decomposition of covariance matrix separately and keep
%only that many dimensions such that params.per_var of the variance
%gets explained:
t1 = clock;
[pcaE, pcaD] = fastica(X(params.channels, params.frames),'only','pca',...
    'verbose',params.verbose);
t2 = clock;
fprintf('PCA performed in %g seconds.\n',etime(t2,t1));

[d,ind] = sort(diag(pcaD),1,'descend');

if strcmp(params.estimate,'eigSpectrum')
   nEig = estimatenumberofneurons(d,'linearFit',nnz(params.frames)); 
   params.numOfIC = nEig;
end

%reduce dimensionality to number of neurons
pcaE = pcaE(:,ind(1:nEig));
pcaD = pcaD(ind(1:nEig),ind(1:nEig));

pcsig = pcaE'*X(params.channels, params.frames);
fprintf('Dimensionality reduced to %g out of %g PCs.\n',nEig,...
    nnz(params.channels));


fprintf('Trying to extract %g independent components.\n',params.numOfIC);
t1 = clock;
[A, W] = fastica(pcsig,'g',params.nonlinearity,...
    'approach',params.approach,'numOfIC',params.numOfIC,...
    'verbose',params.verbose);
clear pcsig
%calculate effective As and Ws in original sensor space
A = pcaE*A;
W = W*pcaE';

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
