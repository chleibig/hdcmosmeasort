function [ OL ] = roioverlap( pixelIdxList, varargin )
%[ OL ] = roioverlap( pixelIdxList, varargin )
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
% OPTIONAL call signature: 
%
% [ OL ] = roioverlap( pixelIdxList1, pixelIdxList2 ) if called with two
%          different lists of ROI pixels, comparison is just performed 
%          between, but not within the lists.
% pixelIdxListi - (1 x Ni_ROI) cell arrays as above
%  
% 
% Output
% ======
%           
% [ OL ] = roioverlap( pixelIdxList) ->
% OL - (N_ROI x N_ROI) matrix with OL(k,l) being the fraction of
%      pixels being shared by the k-th and l-th region of interest
%
% [ OL ] = roioverlap( pixelIdxList1, pixelIdxList2) ->
% OL - (N1_ROI x N2_ROI) matrix with OL(k,l) being the fraction of
%      pixels being shared by the k-th region of interest of the first
%      argument and the l-th region of interest of the second argument
%     
%
% christian.leibig@g-node.org, 19.07.13 updated 29.08.13, 11.12.13
%
%

if isempty(varargin)
    %%%%%%%%%%%%%%%%%%%%% All to all comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_ROI = length(pixelIdxList);
    OL = eye(N_ROI);
    %fprintf('Computing %g ROI overlaps...',N_ROI*(N_ROI -1)/2);
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
    %fprintf('done.\n');

else

    %%%%%%%%%%%%% Comparison across two given lists %%%%%%%%%%%%%%%%%%%%%%%
    pixelIdxList1 = pixelIdxList;
    pixelIdxList2 = varargin{1};
    
    N1_ROI = length(pixelIdxList1);
    N2_ROI = length(pixelIdxList2);
    %fprintf('Computing %g ROI overlaps...',N1_ROI*N2_ROI);
    OL = zeros(N1_ROI,N2_ROI);
    for k1 = 1:N1_ROI
        for k2 = 1:N2_ROI
            %normalize to smaller region of interest.
            nk1 = length(pixelIdxList1{k1});
            nk2 = length(pixelIdxList2{k2});
            nNorm = (nk1 <= nk2)*nk1 + ~(nk1 <= nk2)*nk2;
            %OL(k,l) = length(intersect(pixelIdxList{k},pixelIdxList{l}))/nNorm;
            %intersect calls unique twice, union calls unique only once...
            %hence we can optimize the code as follows. (For many ROIs the following
            %step is expensive)
            OL(k1,k2) = ...
                (nk1 + nk2 - length(union(pixelIdxList1{k1},pixelIdxList2{k2})))/nNorm;
        end
    end
    %fprintf('done.\n');

end



 
end

 