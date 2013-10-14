function [ROI, varargout] = ...
    CoG_ROIs(filename_events,maxSensorsPerEvent,minAct,N_ROWS, N_COLS, fillColumnFlag,...
    mergeThr)
%[ROI, varargout] = CoG_ROIs(filename_events,maxSensorsPerEvent,minAct,N_ROWS, N_COLS,
%fillColumnFlag, mergeThr)
%calculates regions of interest based on the connected components 
%of the center of gravity (CoG) sensors that exhibit at least minAct events
%over the entire recording
%
% Input
% =====
%
% filename_events: *.basic_events file, used for CoG and surrounding
% sensors
% maxSensorsPerEvent
% minAct: min activity in # of events, i.e. unitless
% N_ROWS, N_COLS: array extension in sensor coordinates
% fillColumnFlag: if true, skipped sensor columns get filled
% mergeThr: regions of interest with at least mergeThr overlap fraction get
%           merged
% maxSensorsPerEvent
%
% Output
% ======
%   
%    [ROI] = CoG_ROIsCoG_ROIs(filename_events,maxSensorsPerEvent,minAct,N_ROWS, N_COLS,
%                             fillColumnFlag, mergeThr)
%
%            ROI: regions of interest
%
%    [ROI, CC, nEventsPerSensor, minNoEventExistence] = ...
%              CoG_ROIs(filename_events,maxSensorsPerEvent,minAct,...
%                        N_ROWS, N_COLS, fillColumnFlag, mergeThr)
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
% Filter events by validity - not necessary if previous basic event
% detection already does this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nEvents = length(time);

%Construct filter according to number of sensors participating in an event
nSensor = zeros(1,nEvents);
for i = 1:nEvents
    nSensor(i) = length(sscanf(sensors{i},'%*[(] %u %*[,] %u %*[)]'))/2;
end
valid = (nSensor <= maxSensorsPerEvent);

time = time(valid);
amplitude = amplitude(valid);
boss_column = boss_column(valid);
boss_row = boss_row(valid);
sensors = sensors(valid);

fprintf('%g out of %g events are accepted for ROI construction\n',...
    nnz(valid),nEvents);
clear nEvents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct regions of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nEvents = length(time);

nEventsPerSensor = zeros(N_ROWS,N_COLS);
for i = 1:nEvents;
    nEventsPerSensor(boss_row(i),boss_column(i)) = ...
        nEventsPerSensor(boss_row(i),boss_column(i)) + 1;
end


minNoEventExistence = (nEventsPerSensor >= minAct);


%Connected components on CoG sensors
%CC = bwconncomp(minNoEventExistence);
%Or skip connected components on CoG sensors:
CC = struct('Connectivity',8,'ImageSize',[N_ROWS N_COLS],...
    'NumObjects',nnz(minNoEventExistence),...
    'PixelIdxList',{num2cell(find(minNoEventExistence)')});

%For each connected component add all sensors that participate in the
%events which are associated with the respective connected component

%copy instead of appending another field for sensors because LABELMATRIX
%takes CC as it is output from BWCONNCOMP as input

ROI = CC;

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
            if fillColumnFlag
                %make use of MATLAB internal vector representation:
                eventRows = [eventRows(:)';eventRows(:)'];
                eventRows = eventRows(:)';
                eventCols = [eventCols(:)';eventCols(:)'+1];
                eventCols = eventCols(:)';
            end
            eventSensors = sub2ind(CC.ImageSize,eventRows,eventCols);
            sensorsROI = [sensorsROI;eventSensors(:)];
        end
        ROI.time{i} = [ROI.time{i};time(events)];
    end
    ROI.PixelIdxList{i} = unique(sensorsROI);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge regions of interest based on overlap fraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following function is espensive ! :
[ OL ] = roioverlap( ROI.PixelIdxList );

%adjacencyMatrix:
A = OL >= mergeThr;

clear OL

%New ROIs are the connected components of that graph, take code from:
%http://www.mathworks.com/matlabcentral/fileexchange/
%42040-find-network-components
oldROI = ROI;
ROI.PixelIdxList = {};
ROI.time = {};
[ROI.NumObjects,sizes,members,unused] = networkComponents(A);

%merge sensors and time stamps for merge cases.
for i = 1:ROI.NumObjects
    if sizes(i) == 1
        ROI.PixelIdxList{i} = oldROI.PixelIdxList{members{i}};
        ROI.time{i} = oldROI.time{members{i}};
    else
        ROI.PixelIdxList{i} = unique(vertcat(oldROI.PixelIdxList{members{i}}));
        ROI.time{i} = unique(vertcat(oldROI.time{members{i}}));
    end
end
    
clear oldROI


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    varargout{1} = CC;
    varargout{2} = nEventsPerSensor;
    varargout{3} = minNoEventExistence;
end


end