function [G] = cluster_similarity_matrix(S, thr)
% Cluster data points if their similarity, given by 
% the upper triangular similarity matrix S exceeds thr
% for at least one other cluster member

%Comment 11.01.13 : this seems to be equivalent to hierarchical clustering!

% Output: G is a cell array with each cell representing
% a cluster whose member indices are given in a numeric array

g = 2:size(S,1);
G = {1}; % the first item is the seed for the first cluster

while ~isempty(g)
    found = false;
    for i = g %remaining items
        for j = G{end}%iterate through current cluster
          if S(j,i) >= thr || S(i,j) >= thr
              G{end} = [G{end} i]; %append item to current cluster
              g(g == i) = [];
              found = true;
              break; %restart loop because newly added item was not yet checked
              %against remaining items before current i.
          end
        end
        if found
            break;
        end    
    end
    if ~found %open a new cluster
        G = {G{:} g(1)};
        g(1) = [];
    end
        
end  

end