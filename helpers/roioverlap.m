function [ OL ] = roioverlap( pixelIdxList )
%[ OL ] = roioverlap( pixelIdxList )
%roioverlap calculates the overlap between regions of interest
%Normalization is with respect to number of participating sensors of smaller
%region of interest.
%
% Input
% =====
%
% pixelIdxList - (1 x N_ROI) cell array; pixelIdxList{k} is a vector
%                containing all the linear pixel indices of the k-th region
%                of interest
% 
% Output
% ======
%
% OL - (N_ROI x N_ROI) matrix with OL(k,l) being the fraction of
%      pixels being shared by the k-th and l-th region of interest
%
% christian.leibig@g-node.org, 19.07.13 updated 29.08.13
%
%
    
N_ROI = length(pixelIdxList);

OL = eye(N_ROI);
for k = 1:N_ROI
    for l = k+1:N_ROI
        %normalize to smaller region of interest.
        nk = length(pixelIdxList{k});
        nl = length(pixelIdxList{l});
        nNorm = (nk <= nl)*nk + ~(nk <= nl)*nl;
        %OL(k,l) = length(intersect(pixelIdxList{k},pixelIdxList{l}))/nNorm;
        %intersect calls unique twice, union calls unique only once...
        %hence we can optimize the code as follows. (For many ROIs the following
        %step is expensive)
        OL(k,l) = (nk + nl - length(union(pixelIdxList{k},pixelIdxList{l})))/nNorm;
        OL(l,k) = OL(k,l);
    end
end



end

 