function [ filledVector ] = fillvector( seedVector, horizon )
%fillvector sets all elements in seedVector within horizon to the left
%and to the right of true values to true. Seed values that are less than
%horizon from the borders apart are ignored
%
% Input
% =====
% 
% seedVector - boolean vector
% horizon - single sided extension to be filled with truth values
% 
% Output
% ======
%
% filledVector - boolean vector of length(seedVector)
%
% christian.leibig@g-node.org, 16.07.13
%
%

filledVector = seedVector;

seedIdxs = find(seedVector);
for i = 1:length(seedIdxs)
    %seeds at the borders are ignored
    if (seedIdxs(i)-horizon > 0) & (seedIdxs(i)+horizon <= length(seedVector))
        filledVector(seedIdxs(i)-horizon:seedIdxs(i)+horizon) = true;
    end
end

