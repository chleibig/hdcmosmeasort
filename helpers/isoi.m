function [IsoI] = isoi(featP,featQ)
% [IsoI] = isoi(featP,featQ) evaluate isolation information of 1-D feature
% samples from distribution P vs. 1-D feature samples from distribution Q.
% 
% Distributions P and Q are estimated from the given samples featP and
% featQ via Kernel Density Estimation.
% To avoid zero probabilities for evaluation the Kullback-Leibler
% divergence, a small eps is added to all probabilities.
%
% Isolation information measures were proposed and evaluated to assess 
% the quality of spike sorting results in:
% Neymotin et al., 2011: Measuring the quality of neuronal identification
%                        in ensemble recordings, Journal of Neuroscience
%
% NOTE: Neymotin et al., 2011 compute IsoI for features of arbitrary
% dimensionality for which they use nearest neighbour estimates for 
% the Kullback-Leibler divergence
%
% IsoIBg = isoi(featC,featBg)
% IsoINN = isoi(featC,featNN)
%
% Author: Christian Leibig, 08.05.14
%
% Dependencies:
% 
%  kde.m   - http://www.mathworks.com/matlabcentral/fileexchange/
%            14034-kernel-density-estimator/content/kde.m
%
%  kldiv.m - http://www.mathworks.com/matlabcentral/fileexchange/
%            13089-kldiv/content/kldiv.m 
% 

DEBUG = false;

if length(featP) < 2 || length(featQ) < 2
    warning('Not enough samples for either P or Q, IsoI is set to NaN.');
    IsoI = NaN;
    return;
end

maxFeat = max(max(featP),max(featQ));
minFeat = min(min(featP),min(featQ));

%construct a range that captures both distributions and extend it by
%the default fraction from kde.m
RANGE = maxFeat - minFeat;
MIN = minFeat - RANGE/10;
MAX = maxFeat + RANGE/10;
n=2^14; %default from kde
%estimate distributions with kde.
[unused,densityP,xmesh,unused] = kde(featP,n,MIN,MAX);
[unused,densityQ,xmesh,unused] = kde(featQ,n,MIN,MAX);

%normalize to pdfs
pdfP = densityP/sum(densityP);
pdfQ = densityQ/sum(densityQ);
%add eps to avoid zeros for computing kldiv
eps = 1e-10;
pdfP = (pdfP + eps)/(1 + length(pdfP)*eps);
pdfQ = (pdfQ + eps)/(1 + length(pdfQ)*eps);
%compute Kullback-Leibler divergence (a different estimate is used in 
%Neymotin et al., 2011 !)
KLpQ = kldiv(xmesh',pdfP,pdfQ);
KLqP = kldiv(xmesh',pdfQ,pdfP);

IsoI = KLpQ * KLqP / (KLpQ + KLqP);

if DEBUG
    figure;plot(xmesh,pdfP,'b');hold on;plot(xmesh,pdfQ,'r');
end