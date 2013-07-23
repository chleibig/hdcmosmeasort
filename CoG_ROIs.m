function [ROI, varargout] = ...
    CoG_ROIs(filename_events,minAct,N_ROWS, N_COLS)
%CoG_ROIs calculates regions of interest based on the connected components 
%of the center of gravity (CoG) sensors that exhibit at least minAct events
%over the entire recording
%
% Input
% =====
%
% filename_events: *.basic_events file, used for CoG and surrounding
% sensors
% minAct: min activity in # of events, i.e. unitless
% N_ROWS, N_COLS: array extension in sensor coordinates
%
%
% Output
% ======
%   
%    [ROI] = CoG_ROIs(filename_events,minAct,N_ROWS, N_COLS)
%
%            ROI: regions of interest
%
%    [ROI, CC, nEventsPerSensor, minNoEventExistence] = ...
%              CoG_ROIs(filename_events,minAct,N_ROWS, N_COLS)
% 
%            ROI: regions of interest
%            CC: connected components as output from thresholded CoG image
%            nEventsPerSensor: original CoG image
%            minNoEventExistence: thresholded CoG image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get centers of gravity of events and sensors on which they are detected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%from the *.basic_events list:
[time, amplitude, boss_column, boss_row, sensors] = textread(filename_events,...
    '%*s %*u %f %f %u %u %s','delimiter','\t','headerlines',1);
if length(sensors{1}) <= 1
    error('Please guarantee that reasonable sensor values are present in %s',filename_events);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nEvents = length(amplitude);

nEventsPerSensor = zeros(N_ROWS,N_COLS);
for i = 1:nEvents;
    nEventsPerSensor(boss_row(i),boss_column(i)) = ...
        nEventsPerSensor(boss_row(i),boss_column(i)) + 1;
end




minNoEventExistence = (nEventsPerSensor >= minAct);


%Connected components on CoG sensors
CC = bwconncomp(minNoEventExistence);

%For each connected component add all sensors that participate in the
%events which are associated with the respective connected component

%copy instead of appending another field for sensors because LABELMATRIX
%takes CC as it is output from BWCONNCOMP as input

ROI = CC;
tic;
for i = 1:CC.NumObjects
    %Get CoG sensors of CC(i)
    sensorsROI = CC.PixelIdxList{i};
    [cRows,cCols] = ind2sub(CC.ImageSize,CC.PixelIdxList{i});
    %consider computational scalability of following block!!!
    ROI.time{i} = [];
    for j = 1:length(cRows)%for all CoG sensors in ROI
        %Get the participating event indices:
        events = find((boss_row == cRows(j)) & (boss_column == cCols(j)));
        %Append for each event the not yet present eventSensors
        for k = 1:length(events);
            [colsRowsInterleaved] = sscanf(sensors{events(k)},...
                '%*[(] %u %*[,] %u %*[)]');
            eventRows = colsRowsInterleaved(2:2:end);
            eventCols = colsRowsInterleaved(1:2:end);
            eventSensors = sub2ind(CC.ImageSize,eventRows,eventCols);
            sensorsROI = union(sensorsROI,eventSensors);
        end
        ROI.time{i} = [ROI.time{i};time(events)];
        clear events
    end
    ROI.PixelIdxList{i} = sensorsROI;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+

if nargout > 1
    varargout{1} = CC;
    varargout{2} = nEventsPerSensor;
    varargout{3} = minNoEventExistence;
end


end