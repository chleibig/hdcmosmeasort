function [ OL ] = roioverlap( pixelIdxList )
%roioverlap calculates the overlap between regions of interest
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
% christian.leibig@g-node.org, 19.07.13
%
%
    
N_ROI = length(pixelIdxList);

OL = eye(N_ROI);
for k = 1:N_ROI
    for l = k+1:N_ROI
        OL(k,l) = ...
            length(intersect(pixelIdxList{k},pixelIdxList{l}))/...
            length(union(pixelIdxList{k},pixelIdxList{l}));
        OL(l,k) = OL(k,l);
    end
end



end

 