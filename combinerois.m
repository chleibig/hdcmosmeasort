function [roisOut, N_DUPL] = combinerois(roisIn, OL, sr, data,...
                                 sensor_rows, sensor_cols, t_s, t_jitter,...
                                 coin_thr, sim_thr, show, interactive)
%combinerois extracts valid units from all regions of interest (roisIn) by
%checking all pairwisely overlapping combinations as specified in OL.
%Pairs (i,j) are treated sequentially, order by decreasing overlap fraction OL(i,j).
%
% Input
% =====
%
% roisIn - array of structs, with roisIn(i) being a struct with the information 
%        of one region of interest stored in the following fields:
%        
%        *.sensor_rows, *.sensor_cols - vectors       
%        *.units, which in turn contain at least the fields
%        *.units.time
%        
%
% OL -   (N_ROI x N_ROI) - matrix, with OL(i,j) specifying the fraction of
%        overlap between regions of interest i and j.
%
% sr - sampling rate in kHz
% data - (N_ROW x N_COL x N_FRAMES) array
% sensor_rows, sensor_cols - vectors of sensor coordinates of data
% t_s - maximal overall, pairwise shift of spiketrains to detect coincident
%      spikes in ms
% t_jitter - two spikes are considered as coincident if they are maximally
%           t_s + t_jitter ms apart from each other 
% coin_thr - minimum fraction of coincident spikes
% sim_thr - minimum similarity of average waveforms (STA)
%
% show, interactive - flags to control graphical output / interactivity
% 
%
% Output
% ======
%
% roisOut - the same as roisIn, but with the duplicates moved to *_dupl
%           fields in the respective regions of interest
%
% N_DUPL - the total number of duplicates found
%
% christian.leibig@g-node.org, 22.07.13
%

%CONCEPT: iterate through OL combinations:
% -> duplicate checks
% -> storage units / i.e. separation between removals and valid units...

if length(roisIn) == 1
    roisOut = roisIN;
    fprintf('Only one ROI present - nothing to combine!\n');
    return
end


if show
   figure; 
   imagesc(OL);colorbar;
   title('overlap fraction between regions of interest (ROI)');
   xlabel('ROI index');
   ylabel('ROI index');
   axis square;
end


%sort overlapping roi combinations in decreasing order and get ROI indices
%I and J:
OL = triu(OL,1);
[olSorted,ind] = sort(OL(:),1,'descend');
[I,J] = ind2sub(size(OL),ind(olSorted > 0));

N_OL_PAIRS = length(I);
N_DUPL = 0;
for i = 1:N_OL_PAIRS
    
    if ( (length(roisIn(I(i)).units) > 1) &&  (length(roisIn(J(i)).units) > 1) )
        %Select data region on which to calculate the STAs for the
        %duplicate check.
        sensor_rows_IJ = union(roisIn(I(i)).sensor_rows, roisIn(J(i)).sensor_rows);
        sensor_cols_IJ = union(roisIn(I(i)).sensor_cols, roisIn(J(i)).sensor_cols);
                
        dataIJ = data(...
            (sensor_rows_IJ(1) <= sensor_rows) & ...
            (sensor_rows <= sensor_rows_IJ(end)),...
            (sensor_cols_IJ(1) <= sensor_cols) & ...
            (sensor_cols <= sensor_cols_IJ(end)),:);
        
        %Get duplicates.
        [duplicate_pairs, IiDupl, JiDupl] = checkforinterroiduplicates(...
            roisIn(I(i)).units, roisIn(J(i)).units, sr, dataIJ, ...
            t_s, t_jitter, coin_thr, sim_thr, show, interactive);
        N_DUPL = N_DUPL + length(duplicate_pairs);
        %Remove duplicates in ROI I of pair i
        if ~isempty(IiDupl)
            remove = false(length(roisIn(I(i)).units),1);
            remove(IiDupl) = true;
            
            %S_dupl
            if ~isfield(roisIn(I(i)),'S_dupl')
                roisIn(I(i)).S_dupl = roisIn(I(i)).S(remove,:);
            else
                roisIn(I(i)).S_dupl = [roisIn(I(i)).S_dupl;roisIn(I(i)).S(remove,:)];
            end
            roisIn(I(i)).S = roisIn(I(i)).S(~remove,:);
            
            %A_dupl
            if ~isfield(roisIn(I(i)),'A_dupl')
                roisIn(I(i)).A_dupl = [];
            else
                roisIn(I(i)).A_dupl = cat(2,roisIn(I(i)).A_dupl,roisIn(I(i)).A_tau(:,remove,:));
            end
            roisIn(I(i)).A_tau = roisIn(I(i)).A_tau(:,~remove,:);
            
            %units_dupl
            if ~isfield(roisIn(I(i)),'units_dupl')
                roisIn(I(i)).units_dupl = roisIn(I(i)).units(remove);
            else
                roisIn(I(i)).units_dupl = ...
                    [roisIn(I(i)).units_dupl roisIn(I(i)).units(remove)];
            end            
            roisIn(I(i)).units = roisIn(I(i)).units(~remove);            
            clear remove
        end
        
        %Remove duplicates in ROI J of pair i
        if ~isempty(JiDupl)
            remove = false(length(roisIn(J(i)).units),1);
            remove(JiDupl) = true;
            
            %S_dupl
            if ~isfield(roisIn(J(i)),'S_dupl')
                roisIn(J(i)).S_dupl = roisIn(J(i)).S(remove,:);
            else
                roisIn(J(i)).S_dupl = [roisIn(J(i)).S_dupl;roisIn(J(i)).S(remove,:)];                                
            end
            roisIn(J(i)).S = roisIn(J(i)).S(~remove,:);
            
            %A_dupl
            if ~isfield(roisIn(J(i)),'A_dupl')
                roisIn(J(i)).A_dupl = [];
            else
                roisIn(J(i)).A_dupl = cat(2,roisIn(J(i)).A_dupl,roisIn(J(i)).A_tau(:,remove,:));
            end
            roisIn(J(i)).A_tau = roisIn(J(i)).A_tau(:,~remove,:);
            
            %units_dupl
            if ~isfield(roisIn(J(i)),'units_dupl')
                roisIn(J(i)).units_dupl = roisIn(J(i)).units(remove);
            else
                roisIn(J(i)).units_dupl = ...
                    [roisIn(J(i)).units_dupl roisIn(J(i)).units(remove)];
            end            
            roisIn(J(i)).units = roisIn(J(i)).units(~remove);      
            
            clear remove
        end

    end
    
end

roisOut = roisIn;

end