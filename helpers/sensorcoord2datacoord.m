function [ dataRows, dataCols ] = sensorcoord2datacoord(...
    sensorRows,sensorCols, SENSOR_ROW_IDS, SENSOR_COL_IDS)
%sensorcoord2datacoord
%
%
% Input
% =====
%
% sensorRows,sensorCols - the sensor coordinates to be converted
%                         to data coordinates (both vectors of the same length)
%
% SENSOR_ROW_IDS, SENSOR_COL_IDS - vectors containing the sensor
%                                  coordinates of the data block
%
% Output
% ======
%
% dataRows, dataCols - indices into a data block of size
%                      (length(SENSOR_ROW_IDS) x length(SENSOR_COL_IDS))
%
% christian.leibig@g-node.org, 18.07.13
%

if (length(sensorRows) ~= length(sensorCols))
    error('sensorRows and sensorCols must have the same length!\n');
end

dataRows = zeros(size(sensorRows));
dataCols = zeros(size(sensorCols));

for i = 1:length(sensorRows)
    try
        dataRows(i) = find(sensorRows(i) == SENSOR_ROW_IDS);
        dataCols(i) = find(sensorCols(i) == SENSOR_COL_IDS);
    catch
        error('sensor [%g %g] not found.\n',sensorRows(i),sensorCols(i));
    end
end


end

