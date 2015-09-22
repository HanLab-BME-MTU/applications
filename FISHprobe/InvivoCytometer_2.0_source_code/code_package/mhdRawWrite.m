function mhdRawWrite(data, spacing, filename)
% Write a RAW file and associated MHD.
% Format: mhdRawWrite(data, spacing, filename)
%
% Revisions:
%   By: Paul (07/01/09)
%       - add spacing input
%

% Already have an edatatension? Remove it...
[pathstr, name] = fileparts(filename);
filename = fullfile(pathstr, name);

mhdFile = [ filename '.mhd' ];
rawFile = [ filename '.raw' ];

% Data type
elType = class(data);
switch elType
    case 'logical'
        elType = 'MET_UCHAR';
    case 'uint8'
        elType = 'MET_UCHAR';
    case 'int8'
        elType = 'MET_CHAR';
    case 'uint16'
        elType = 'MET_USHORT';
    case 'int16'
        elType = 'MET_SHORT';
    case 'uint32'
        elType = 'MET_ULONG';
    case 'int32'
        elType = 'MET_LONG';
    case 'single'
        elType = 'MET_FLOAT';
    case 'double'
        elType = 'MET_DOUBLE';
    otherwise
        error('Unknown data type : %s\nPlease modify mhdRawWrite.', elType);
end

% Write mhd file
mhd = fopen(mhdFile, 'wt');
if mhd == -1
    error('Could not open file "%s".\n', mhdFile);
end
fprintf(mhd, 'ObjectType = Image\n');
fprintf(mhd, 'BinaryData = True\n');
fprintf(mhd, 'BinaryDataByteOrderMSB = False\n');
if size(data,3) == 1
    fprintf(mhd, 'NDims = 2\n');
    fprintf(mhd, 'DimSize = %d %d\n', size(data,2), size(data,1));
    fprintf(mhd, 'ElementSpacing = %.12g %.12g\n', spacing(1), spacing(2));
else
    fprintf(mhd, 'NDims = 3\n');
    fprintf(mhd, 'DimSize = %d %d %d\n', size(data,2), size(data,1), size(data,3));
    fprintf(mhd, 'ElementSpacing = %.12g %.12g %.12g\n', spacing(1), spacing(2), spacing(3));
end
fprintf(mhd, 'ElementType = %s\n', elType);
fprintf(mhd, 'ElementDataFile = %s\n', [name '.raw']);
fclose(mhd);

% Handle different dimensions (permute x / y)
data = permute(data, [2 1 3]);

% Write raw data
raw = fopen(rawFile, 'wb');
if raw == -1
    error('Could not open file "%s".\n', rawFile);
end
fwrite(raw, data, class(data));
fclose(raw);

end
