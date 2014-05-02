function [data] = readDataBlock(filename,...
                    firstRow,lastRow,firstCol,lastCol,firstFrame,lastFrame)
% [data] = readDataBlock(filename,...
%                   firstRow,lastRow,firstCol,lastCol,firstFrame,lastFrame)
% 
% read a block from HDCMOSMEA data stored according to fileversion
% 'NH5_1.0.0' as specified in NH5_specification.txt
% 
% Input
% =====
%
% 2nd to last parameters index into the 3-dim Data array in MATLAB notation
% (starting with '1')
%
% Output
% ======
% 
% data - [NROW, NCOL, NFRAMES] 
% 
% christian.leibig@g-node.org, 04.04.14
%



try
    fileVersion = h5attget(filename,'/Metadata','FileVersion');
    fileVersion = fileVersion{1};
catch
    error('Fileformat unknown or file not found!');  
end

switch fileVersion
    case 'NH5_1.0.0'
      h5PathData = '/Data';
      chipType = h5attget(filename,'/Metadata','ChipType');
    otherwise
        error(['Fileversion ' fileVersion ' is unknown']);
end
        


rmode = 'H5F_ACC_RDONLY';
plist = 'H5P_DEFAULT';
fid = H5F.open(filename, rmode, plist);

dset_id = H5D.open(fid,h5PathData);

NROWS = lastRow - firstRow + 1;
NCOLS = lastCol - firstCol + 1;
NFRAMES = lastFrame - firstFrame + 1;

if ( NROWS == 1 || NCOLS == 1 || NFRAMES == 1)
    %There is an unresolved issue in case of singleton dimensions when
    %reading in data from a hdf5 hyperslab to matlab, probably attributable
    %to the difference of Matlab's column-major order vs hdf5's 
    %row-major order. 07-04-14
    error(['At least one of the dimensions is one, '...
        'data was found to be inconsistent in that case!']);
end

h5_dims = [NROWS, NCOLS, NFRAMES];
h5_maxdims = h5_dims;
h5_start = [firstRow-1, firstCol-1, firstFrame-1];
h5_stride = [1 1 1];
h5_count = [1 1 1];
h5_block = [NROWS, NCOLS, NFRAMES];

mem_space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
file_space_id = H5D.get_space(dset_id);

H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',h5_start,h5_stride,...
                                                    h5_count,h5_block);

%type_id = H5D.get_type(dset_id);
data = H5D.read(dset_id,'H5ML_DEFAULT',mem_space_id,file_space_id,plist);

%H5T.close(type_id);
H5S.close(file_space_id);
H5S.close(mem_space_id);
H5D.close(dset_id);
H5F.close(fid);

data = permute(data, [3 2 1]);
if ~strcmp(class(data),'double');data = double(data);end;


if strcmp(chipType,'NCA')
  data = 1000 * data; %conversion to mV scale
end

end