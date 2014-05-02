function [metadata] = read_metadata(filename)
%[metadata] = read_metadata(filename)
%
% reads metadata from HDCMOSMEA data stored according to fileversion
% 'NH5_1.0.0' as specified in NH5_specification.txt
% 
% Input
% =====
%
% filename - string
%
% Output
% ======
%
% metadata - structure
%        metadata.sensorRows - list of sensor row coordinates 
%        metadata.sensorCols - list of sensor column coordinates  
%        metadata.frameStartTimes - in ms
%        metadata.sr - sampling rate in kHz
%        metadata.sensorPitch - in µm
%        metadata.chipType - string
%
%
% christian.leibig@g-node.org, 04.04.14
%


try
    fileVersion = hdf5read(filename,'/Metadata/FileVersion');
    fileVersion = fileVersion.Data;
catch
    error('Fileformat unknown or file not found!');  
end

switch fileVersion
    case 'NH5_1.0.0'      
        metadata.sensorRows = hdf5read(filename,'/Metadata/RowList');
        metadata.sensorCols = hdf5read(filename,'/Metadata/ColumnList');
        metadata.frameStartTimes = hdf5read(filename,'/Metadata/FrameStartTimes');
      
        metadata.sr = length(metadata.frameStartTimes)/...
            (metadata.frameStartTimes(end) - ...
             metadata.frameStartTimes(1));%in kHz
        chipType = h5attget(filename,'/Metadata','ChipType');
        chipType = chipType{1};
        switch chipType
            case 'G1183'
                metadata.sensorPitch = 7.4;
            case 'NCA'
                metadata.sensorPitch = 16.0;
            otherwise
                error(strcat('Pitch of ',chipType,' is unknown'));
        end
        metadata.chipType = chipType;
    otherwise
        error(['Fileversion ' fileVersion ' is unknown']);
end
        


end