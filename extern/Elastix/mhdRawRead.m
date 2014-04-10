function [data] = mhdRawRead(filename, cast)
% [data] = mhdRawRead(filename, cast)
%
% Read an ITK-style MHD and load RAW data from its associated RAW file
% Doesn't take LSB/MSB into account, but has some support for other
% datatypes (add support in switch statement below).
%
% Currently supports: MET_UCHAR, MET_CHAR, MET_USHORT, MET_SHORT, MET_LONG,
% MET_ULONG, MET_UINT, MET_INT, MET_FLOAT, MET_DOUBLE
%
% If cast is specified, will cast result to that type regardless of
% type stored in file.
%
% Based on function RAWfromMHD by Sean
% Revisions:
%   By: Uday Kurkure (06/26/07) [ver 0.9.0]
%       - output [data] is a structure now.
%               data.im, data.Spacing
%   By: Paul Balanca (07/01/09)
%       - use regular expressions to parse mhd file.

% Open MHD file
mhd = fopen(filename, 'rt');
if mhd == -1
    error('Could not open file "%s".\n', filename);
end

% Read lines
ElementNumberOfChannels = 1;
while (~feof(mhd))
    % e.g., 'NDims = 3'
    line = fgetl(mhd);
    
    % Split the line
    matchNames = regexp(line, '^(?<name>\S*)\s*=\s*(?<value>.*)', 'names');
    matchValues = regexp(matchNames.value, '\S*', 'match');
    
    switch matchNames.name
        case 'NDims'
            NDims = str2double(matchValues{1}); %#ok<NASGU>
        case 'ElementType'
            elType = matchValues{1};
        case 'ElementSpacing'
            Spacing = str2double(matchValues(1:end));
        case 'DimSize'
            DimSize = str2double(matchValues(1:end));
        case 'ElementDataFile'
            dataFile = matchNames.value;
        case 'ElementNumberOfChannels'
            ElementNumberOfChannels = str2double(matchValues{1});
    end
end
fclose(mhd);

% One thing not defined ?
if (~exist('NDims','var') || ~exist('elType','var') || ...
        ~exist('DimSize','var') || ~exist('dataFile','var'))
    error('One or more fields undefined in MHD.\n');
end

% Find matlab data type associated with type here
switch elType
    case 'MET_UCHAR'
        elType = 'uint8';
    case 'MET_CHAR'
        elType = 'int8';
    case 'MET_USHORT'
        elType = 'uint16';
    case 'MET_SHORT'
        elType = 'int16';
    case 'MET_ULONG'
        elType = 'uint32';
    case 'MET_LONG'
        elType = 'int32';
    case 'MET_UINT'
        elType = 'uint32';
    case 'MET_INT'
        elType = 'int32';
    case 'MET_FLOAT'
        elType = 'float32';
    case 'MET_DOUBLE'
        elType = 'double';
    otherwise
        error('Unknown data type : %s\nPlease add this type to mhdRawRead.\n', elType);
end

% Forcing cast ?
if (exist('cast','var'))
    elType = sprintf('%s=>%s', elType, cast);
else
    elType = [ '*' elType ];
end

% Open file ...
pathstr = fileparts(filename);
dataFileFull = fullfile(pathstr, dataFile);
rawfile = fopen(dataFileFull, 'rb');
if rawfile == -1
    error('Failure : Could not open "%s".\n', dataFile);
end

% Read data
dataSize = [ElementNumberOfChannels DimSize];
data.im = fread(rawfile, prod(dataSize), elType);
data.im = reshape(data.im, dataSize);
fclose(rawfile);

% Handle different dimensions (permute x / y)
data.im = squeeze(permute(data.im, [3 2 4 1]));
data.Size = size(data.im);

% Save in struct
data.RelFolderPath = pathstr;
data.fname = filename;
data.Spacing = Spacing;

end
