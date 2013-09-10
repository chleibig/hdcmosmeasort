function [data] = read_data(filename)
%read_data reads high-density array data
%
% Input
% =====
%
% filename - string
%
% Output
% ======
%
% data - structure
%        data.X - (N_ROW, N_COL, N_FRAMES) - array (double)
%        data.sensorRows - list of sensor row coordinates 
%        data.sensorCols - list of sensor column coordinates  
%        data.sr - sampling rate in kHz
%        data.frameStartTimes - in ms
%        data.sensorPitch - in µm
%
%
% christian.leibig@g-node.org, 12.07.13
%


try
    fileVersion = h5attget(filename,'/Metadata','FileVersion');
    fileVersion = fileVersion{1};
catch
    error(['Format ' fileVersion ' of ' filename ' is unknown']);  
end

switch fileVersion
    case 'NH5_1.0.0'
      
        dataset = hdf5load(filename);
        
        data.X = permute(dataset.Data, [3 2 1]);
        if ~strcmp(class(data.X),'double');data.X = double(data.X);end;
        
        data.sensorRows = dataset.Metadata.RowList;
        data.sensorCols = dataset.Metadata.ColumnList;
        data.sr = length(dataset.Metadata.FrameStartTimes)/...
            (dataset.Metadata.FrameStartTimes(end) - ...
            dataset.Metadata.FrameStartTimes(1));%in kHz
        data.frameStartTimes = dataset.Metadata.FrameStartTimes; %in ms
        chipType = h5attget(filename,'/Metadata','ChipType');
        chipType = chipType{1};
        switch chipType
            case 'G1183'
                data.sensorPitch = 7.4;
            otherwise
                error(strcat('Pitch of ',chipType,' is unknown'));
        end
        
    otherwise
        error(['Fileversion ' fileVersion ' is unknown']);
end
        


end