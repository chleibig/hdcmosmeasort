function [units] = combinerois(ROIs, OL)
%combinerois extracts valid units from all regions of interest (ROIs) by
%checking all pairwisely overlapping combinations as specified in OL.
%Pairs (i,j) are treated sequentially, order by decreasing overlap fraction OL(i,j).
%
% Input
% =====
%
% ROIs - array of structs, with ROIs(i) being a struct with the information 
%        of one region of interest stored in the following fields:
%
% OL -   (N_ROI x R_ROI) - matrix, with OL(i,j) specifying the fraction of
%        overlap between regions of interest i and j.
%
% Output
% ======
%
% units - 
%
% christian.leibig@g-node.org, 22.07.13
%

%CONCEPT: iterate through OL combinations:
% -> duplicate checks
% -> storage units / i.e. separation between removals and valid units...














end