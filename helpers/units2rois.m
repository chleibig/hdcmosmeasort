function [ROIs] = units2rois(units, ROIs)
% [ROIs] = units2rois(units, ROIs) 
%
% Note: stop using variable units after pushing state information back
%       to variable ROIs with this function because indices are not
%       are not consistent anymore. If you want to have a new vector with
%       all units obtain it with the function 'unitsfromrois'
%
% christian.leibig@g-node.org, 2015-04-03
%

for i = 1:length(ROIs);
    %delete <-> state '4'
    toDelete = [units.k] == i & [units.state] == 4;
    
    if ~isfield(ROIs(i),'units_del');
        ROIs(i).units_del = [];
    end
    ROIs(i).units_del = ...
                       [ROIs(i).units_del units(toDelete)];
    
	%map toDelete to index range in ROI because S_del and A_del are taken
	%from ROI
    toDelete = toDelete([units.k] == i);

    if ~isfield(ROIs(i),'S_del');
        ROIs(i).S_del = [];
    end
    ROIs(i).S_del = [ROIs(i).S_del; ROIs(i).S(toDelete,:)];
    
    if ~isfield(ROIs(i),'A_del'); 
        ROIs(i).A_del = [];
    end
    ROIs(i).A_del = cat(2,ROIs(i).A_del,...
                          ROIs(i).A_tau(:,toDelete,:));    
    
    clear toDelete
    
    %single <-> state '3'
    toSave = [units.k] == i & [units.state] <= 3;
    ROIs(i).units = units(toSave);

    %map toSave to index range in ROI because S and A_tau are taken again
    %from ROI
    toSave = toSave([units.k] == i);
    
    ROIs(i).S = ROIs(i).S(toSave,:);
    ROIs(i).A_tau = ROIs(i).A_tau(:,toSave,:);
    
    clear toSave                                                          
end


end