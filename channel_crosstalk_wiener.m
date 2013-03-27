function [ CM ] = channel_crosstalk_wiener( X, order, SD, d_max )
%CHANNEL_CROSSTALK_WIENER estimate crosstalk between
%channels
%
% Inputs
% ======
%
% X: (channels,samples) matrix
% order: prediction order in samples
% SD: (channels, channels) matrix with the entry SD(i,j) being the distance
% between channels i and j in µm
% d_max: maximal distance of channels in µm for the channel crosstalk to be
% computed
% 
% 
% Outputs
% =======
%
% CM: upper triangular matrix containing the maximum crosstalk between 
%     channel i and j 
%
% author: Christian Leibig, 27.03.13
% 
% references: Mads Dyrholm: "Independent Component Analysis in a Convoluted
%            World", PhD Thesis, 2005, chapter 5.1.1
%
%            J. Granger: Investigating causal relations by econometric 
%            models and cross-spectral methods, "Econometrica",
%            37(3):424-438, 1969

fprintf('Estimating crosstalk...');
[D,T] = size(X);
CM = zeros(D,D);
t1 = clock;
%non-vectorized code (alternatively one could compute all the pairwise
%cross- and autocorrelations for the wiener filter in advance (xcorr 
%supports that)):
for i=1:D
    for j=i+1:D
        if SD(i,j) > d_max; CM(i,j) = 0;
        else
            [ w ] = GetWienerFilterCoeff( X(j,:), X(i,:), order );
            i_from_j = conv(w,X(j,:));
            i_from_j = i_from_j(1:end-order);
            clear w 
            [ w ] = GetWienerFilterCoeff( X(i,:), X(j,:), order );
            j_from_i = conv(w,X(i,:));
            j_from_i = j_from_i(1:end-order);
            clear w
            ci = var(i_from_j)/var(X(i,:));
            cj = var(j_from_i)/var(X(j,:));
            
            CM(i,j) = max(ci,cj);           
        end
    end
end
t2 = clock;
fprintf('done in %g sec.\n',etime(t2,t1));

end