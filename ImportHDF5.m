function [ fileinfo, data ] = ImportHDF5( filename )
%ImportHDF5 Get info (about Data and Metadata) from HDF5 file

% input_args:
% filename

% output_args:
% fileinfo: struct, that contains all the information about the file (e.g.
% Sampling Rate, Sensor Lists, etc.)
% data: the actual data
fileinfo = hdf5info(filename);
data = hdf5read(filename, '/Data');

end

