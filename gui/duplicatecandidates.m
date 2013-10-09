function duplicatecandidates(currentUnit, units, ROIs, params, maxDist, data)
%duplicatecandidates(currentUnit, units, ROIs, params, maxDist, data)



%Idea: compare units based on spatial neighbourhood (realised via connected
%components on thresholded ED matrix) with each other. Within in each such
%component compute the STAs on the union of sensors of all participating 
%ROIs.

[ ED ] = euclideandistance([units.boss_row], [units.boss_col], params.pitch, params.pitch);

[NumObjects,sizes,members,unused] = networkComponents(ED <= maxDist);

membersOfCurrentUnit = members{cellfun(@(x) ismember(currentUnit,x),members)};

%skip units marked as 'to be deleted', i.e. with state == 4
membersOfCurrentUnit = ...
    membersOfCurrentUnit([units(membersOfCurrentUnit).state] < 4);

if any([units(membersOfCurrentUnit).state] == 1) % 1 <-> 'unchecked'
    %duplicate test is only worth to do if thresholds of participating
    %units are verified
    msgbox(['Verify thresholds of the following ICs or mark them as '...
        char(39) 'to be saved' char(39) ' or '...
        char(39) 'to be deleted' char(39) ' first: '...
        num2str(membersOfCurrentUnit([units(membersOfCurrentUnit).state] == 1))]);
    return
end

%Compute STAs on set of participating ROIs within each connected
%component.
roiSet = unique([units(membersOfCurrentUnit).k]);
if length(roiSet) > 1
    %combine rois here....
    combinedSensorRows = unique(cat(1,ROIs(roiSet).sensor_rows));
    combinedSensorCols = unique(cat(1,ROIs(roiSet).sensor_cols));
    dataTmp = data(...
        (combinedSensorRows(1) <= params.sensor_rows) & ...
        (params.sensor_rows <= combinedSensorRows(end)),...
        (combinedSensorCols(1) <= params.sensor_cols) & ...
        (params.sensor_cols <= combinedSensorCols(end)),:);
    for i = 1:length(membersOfCurrentUnit)
        units(membersOfCurrentUnit(i)).STA = ...
            GetSTA(dataTmp, units(membersOfCurrentUnit(i)).time, params.sr, 0);
    end
    clear dataTmp
else
    %only one ROI, we do not need to recompute the STA
end



checkforintraroiduplicates(units(membersOfCurrentUnit), params.sr, ...
                           params.t_s, params.t_jitter, params.coin_thr,...
                           params.sim_thr, 1, 0, membersOfCurrentUnit);


end