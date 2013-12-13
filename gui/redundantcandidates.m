function [numDistinct, sizes, members] = redundantcandidates(units, ...
    ROIs, params, data, varargin)
%redundantcandidates(units, ROIs, params, data, varargin)

%Idea: compare units based on spatial neighbourhood (realised via connected
%components on thresholded ED matrix) with each other. Within in each such
%component compute the STAs on the union of sensors of all participating 
%ROIs.
%
% Input
% =====
%
%  varargin{1} <-> currentUnit, if present only associated connected
%                  is evaluated otherwise all are evaluated

% christian.leibig@g-node.org, 10.10.2013

maxDist = params.d_max;

[ ED ] = euclideandistance([units.boss_row], [units.boss_col], params.pitch, params.pitch);

[numSpatialCCs,sizesSpatial,nbrs,unused] = networkComponents(ED <= maxDist);

fprintf(['Found %g spatially connected components with a maximum\n'...
    'distance of %g µm.\n'],numSpatialCCs,maxDist);
fprintf('The largest one contains %g units.\n',max(sizesSpatial));

if ~isempty(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate and visualize only current spatial CC region %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentUnit = varargin{1};
    fprintf(['Analyzing only the spatially connected component,\n',...
        'unit %g is part of.\n'],currentUnit);
    nbrsOfCurrentUnit = nbrs{cellfun(@(x) ismember(currentUnit,x),nbrs)};
    %skip units marked as 'to be deleted', i.e. with state == 4
    nbrsOfCurrentUnit = ...
        nbrsOfCurrentUnit([units(nbrsOfCurrentUnit).state] < 4);
    if any([units(nbrsOfCurrentUnit).state] == 1) % 1 <-> 'unchecked'
        %duplicate test is only worth to do if thresholds of participating
        %units are verified
        msgbox(['Verify thresholds of the following ICs or mark them as '...
            char(39) 'to be saved' char(39) ' or '...
            char(39) 'to be deleted' char(39) ' first: '...
            num2str(nbrsOfCurrentUnit([units(nbrsOfCurrentUnit).state] == 1))]);
        return
    end
    %Compute STAs on set of participating ROIs within each connected
    %component.
    roiSet = unique([units(nbrsOfCurrentUnit).k]);
    if length(roiSet) > 1
        %combine rois here....
        combinedSensorRows = unique(cat(1,ROIs(roiSet).sensor_rows));
        combinedSensorCols = unique(cat(1,ROIs(roiSet).sensor_cols));
        dataTmp = data(...
            (combinedSensorRows(1) <= params.sensor_rows) & ...
            (params.sensor_rows <= combinedSensorRows(end)),...
            (combinedSensorCols(1) <= params.sensor_cols) & ...
            (params.sensor_cols <= combinedSensorCols(end)),:);
        for i = 1:length(nbrsOfCurrentUnit)
            units(nbrsOfCurrentUnit(i)).STA = ...
                GetSTA(dataTmp, units(nbrsOfCurrentUnit(i)).time, params.sr, 0);
        end
        clear dataTmp
    else
        %only one ROI, we do not need to recompute the STA
    end
    %-- Check for redundancy in one spatially connected component region --
    checkforredundancy(units(nbrsOfCurrentUnit), params.sr, ...
                     params.t_s, params.t_jitter, params.coin_thr,...
                     params.sim_thr, 1, 0, 'unitIDs',nbrsOfCurrentUnit);

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterate over all spatial CC regions %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numDistinct = 0;
    sizes = [];
    members = [];
    for j = 1:numSpatialCCs
        %skip units marked as 'to be deleted', i.e. with state == 4
        nbrs{j} = nbrs{j}([units(nbrs{j}).state] < 4);
        if length(nbrs{j}) < 2
            fprintf(['Spatial component %g contains less than two valid units,\n'...
                'no need to chech redundancy, moving on to next one.\n'],j);
            continue;
        end
        if any([units(nbrs{j}).state] == 1) % 1 <-> 'unchecked'
            %duplicate test is only worth to do if thresholds of participating
            %units are verified
            msgbox(['Verify thresholds of the following ICs or mark them as '...
                char(39) 'to be saved' char(39) ' or '...
                char(39) 'to be deleted' char(39) ' first: '...
                num2str(nbrs{j}([units(nbrs{j}).state] == 1))]);
            return;
        end
        %Compute STAs on set of participating ROIs within each connected
        %component.
        roiSet = unique([units(nbrs{j}).k]);
        if length(roiSet) > 1
            %combine rois here....
            combinedSensorRows = unique(cat(1,ROIs(roiSet).sensor_rows));
            combinedSensorCols = unique(cat(1,ROIs(roiSet).sensor_cols));
            dataTmp = data(...
                (combinedSensorRows(1) <= params.sensor_rows) & ...
                (params.sensor_rows <= combinedSensorRows(end)),...
                (combinedSensorCols(1) <= params.sensor_cols) & ...
                (params.sensor_cols <= combinedSensorCols(end)),:);
            for i = 1:length(nbrs{j})
                units(nbrs{j}(i)).STA = ...
                    GetSTA(dataTmp, units(nbrs{j}(i)).time, params.sr, 0);
            end
            clear dataTmp
        else
            %only one ROI, we do not need to recompute the STA
        end
        fprintf('Analyzing redundancy of spatial component %g ...\n',j);
        [numDistinctJ,sizesJ,membersJ] = checkforredundancy(...
                     units(nbrs{j}), params.sr, ...
                     params.t_s, params.t_jitter, params.coin_thr,...
                     params.sim_thr, 0, 0, 'unitIDs',nbrs{j});
        numDistinct = numDistinct + numDistinctJ;
        sizes = [sizes sizesJ];
        members = [members membersJ];
    end
end


end