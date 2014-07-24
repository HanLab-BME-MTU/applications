function obj = read(fullPath)
disp('Data: Reading data ...');

% ---
% --- Read pixel size from *.ini file if it exists
% ---

iniFilePath = fullPath;
iniFilePath = [iniFilePath(1:end-3) 'ini'];
pixelSize = [];
defaultPixelSize = 141;
if exist(iniFilePath,'file')
    disp('Data: Reading pixel size from ini-file!');
    
    fid = fopen(iniFilePath,'r');
    while ~feof(fid)
        line = fgetl(fid);
        if strcmp(line(1:min(11,end)),'pixel size=')
            pixelSize = str2double(line(12:end));
            fprintf('Data: The pixel size is: %.0f nm\n',pixelSize);
            break;
        end
        
    end
    fclose(fid);
    
    if isempty(pixelSize)
        disp('Data: Could not read pixel size: Using default instead!');
        pixelSize = defaultPixelSize;
    end
else
    disp('Data: Using default pixel size: 141 nm!');
    pixelSize = defaultPixelSize;
end


% ---
% --- Read z-scaling factor from *.opt file if it exists
% ---

optFilePath = fullPath;
optFilePath = [optFilePath(1:end-3) 'opt'];
scalingZ = [];
defaultScalingZ = 1;
if exist(optFilePath,'file')
    disp('Data: Reading z-scaling factor from opt-file!');
    
    fid = fopen(optFilePath,'r');
    while ~feof(fid)
        line = fgetl(fid);
        if strcmp(line(1:min(10,end)),'z-scaling=')
            scalingZ = str2double(line(11:end));
            fprintf('Data: The z-scaling is: %.2f\n',scalingZ);
            break;
        end
        
    end
    fclose(fid);
    
    if isempty(scalingZ)
        disp('Data: Could not read z-scaling: Using default instead!');
        scalingZ = defaultScalingZ;
    end
else
    disp('Data: Using default z-scaling: 1!');
    scalingZ = defaultScalingZ;
end


% ---
% --- Open File
% ---

fid = fopen(fullPath,'r');
headerSize = 16; % Bytes; Offset due to the header (version, status, ...)
numParam = 18; % Number of parameter stored for every point


% ---
% --- Read Header
% ---

A = fread(fid,headerSize/4,'int32');
numPoints = A(4); % Number of points in the list


% ---
% --- Read Floats
% ---

frewind(fid);
A = fread(fid,'float32');
combSelector = 0:numParam:(numPoints-1)*numParam; % Vector to select all values of a certain parameter

tmp.x = A(combSelector+7);
tmp.y = A(combSelector+8);
tmp.z = A(combSelector+22);
tmp.h = A(combSelector+9);
tmp.area = A(combSelector+10);
tmp.width = A(combSelector+11);
tmp.phi = A(combSelector+12);
tmp.Ax = A(combSelector+13);
tmp.bg = A(combSelector+14);
tmp.I = A(combSelector+15);

clear A;


% ---
% --- Read Ints
% ---

frewind(fid);
A = int32(fread(fid,'int32'));

tmp.cat = A(combSelector+16);
tmp.frame = A(combSelector+18);
tmp.length = A(combSelector+19);

clear A;
fclose(fid);

validMol = find(tmp.cat);

f.x = tmp.x(validMol);
f.y = tmp.y(validMol);
f.z = tmp.z(validMol);
f.h = tmp.h(validMol);
f.area =  tmp.area(validMol);
f.width = tmp.width(validMol);
f.phi = tmp.phi(validMol);
f.Ax = tmp.Ax(validMol);
f.bg = tmp.bg(validMol);
f.I = tmp.I(validMol);
f.frame = tmp.frame(validMol);
f.length = tmp.length(validMol);
f.cat = tmp.cat(validMol);

f.x = f.x * pixelSize; % Convert unit to nm
f.y = f.y * pixelSize; % Convert unit to nm

f.z = f.z * scalingZ;

obj = Data();
obj.points = [f.x f.y f.z];
obj.intensity = f.I;
obj.frame = f.frame;

disp('Data: Data read!');

end

