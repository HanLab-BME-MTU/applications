% IceImarisConnectorTestUnit
%
% This is a unit test test the ICE Imaris connector
%
% Copyright Aaron Ponti, 2011 - 2012. All rights reserved
function IceImarisConnectorTestUnit

% ImarisConnector version
% =========================================================================
disp(['Testing IceImarisConnector version ', ...
    IceImarisConnector.version()]);

% Create an ImarisConnector object
% =========================================================================
disp('Create an IceImarisConnector object...');
conn = IceImarisConnector;

% Start Imaris
% =========================================================================
disp('Start Imaris...');
assert(conn.startImaris == 1)

% Test that the connection is valid
disp('Get version...');
assert(conn.getImarisVersionAsInteger() > 0)
disp('Test if connection is alive...');
assert(conn.isAlive() == 1)

% Check the starting index
% =========================================================================
disp('Check starting index...');
assert(conn.indexingStart() == 0);

% Open a file
% =========================================================================
disp('Load file...');
filename = fullfile(fileparts(which(mfilename)), 'PyramidalCell.ims');
conn.mImarisApplication.FileOpen(filename, '');

% Check that there is something loaded
% =========================================================================
disp('Test that the file was loaded...');
assert(conn.mImarisApplication.GetDataSet.GetSizeX > 0)

% Check the extends
% =========================================================================
disp('Check the dataset extends...');
EXTENDS = [-0.114 57.841 -0.114 57.8398 -0.151 20.631];
extends = conn.getExtends();
assert(abs(extends(1) - EXTENDS(1)) < 1e-4)
assert(abs(extends(2) - EXTENDS(2)) < 1e-4)
assert(abs(extends(3) - EXTENDS(3)) < 1e-4)
assert(abs(extends(4) - EXTENDS(4)) < 1e-4)
assert(abs(extends(5) - EXTENDS(5)) < 1e-4)
assert(abs(extends(6) - EXTENDS(6)) < 1e-4)

[minX maxX minY maxY minZ maxZ] = conn.getExtends();
assert(abs(minX - EXTENDS(1)) < 1e-4)
assert(abs(maxX - EXTENDS(2)) < 1e-4)
assert(abs(minY - EXTENDS(3)) < 1e-4)
assert(abs(maxY - EXTENDS(4)) < 1e-4)
assert(abs(minZ - EXTENDS(5)) < 1e-4)
assert(abs(maxZ - EXTENDS(6)) < 1e-4)

% Check the voxel size
% =========================================================================
disp('Check the voxel size...');
VOXELSIZE = [0.2273 0.2282 0.3012];
voxelSize = conn.getVoxelSizes();
assert(abs(voxelSize(1) - VOXELSIZE(1)) < 1e-4)
assert(abs(voxelSize(2) - VOXELSIZE(2)) < 1e-4)
assert(abs(voxelSize(3) - VOXELSIZE(3)) < 1e-4)

[vX vY vZ] = conn.getVoxelSizes();
assert(abs(vX - VOXELSIZE(1)) < 1e-4)
assert(abs(vY - VOXELSIZE(2)) < 1e-4)
assert(abs(vZ - VOXELSIZE(3)) < 1e-4)

% Check the dataset size
%
%   X = 255
%   Y = 254
%   Z =  69
%   C =   1
%   T =   1
%
% =========================================================================
disp('Check the dataset size...');
DATASETSIZE = [255 254 69 1 1];
sizes = conn.getSizes();
assert(all(sizes == DATASETSIZE) == 1);

[sizeX sizeY sizeZ sizeC sizeT] = conn.getSizes();
assert(all([sizeX sizeY sizeZ sizeC sizeT] == DATASETSIZE) == 1);

% Get a spot object, its coordinates and check the unit conversions
% =========================================================================
disp('Count all children at root level...');
children = conn.getAllSurpassChildren(0); % No recursion
assert(numel(children) == 4);

% If the casting in getAllSurpassChildren() works, spot is an actual
% spot object, and not an IDataItem. If the casting worked, the object will
% have a method 'GetPositionsXYZ'.
disp('Test autocasting...');
child = conn.getAllSurpassChildren(0, 'Spots');
spot = child{ 1 };
assert(ismethod(spot, 'GetPositionsXYZ'));

% Get the coordinates
pos = spot.GetPositionsXYZ;

% These are the expected spot coordinates
disp('Check spot coordinates and conversions units<->pixels...');
POS = [
    18.5396    1.4178    8.7341
    39.6139   14.8819    9.0352
    35.1155    9.4574    9.0352
    12.3907   21.6221   11.7459 
   ];

assert(all(all(abs(pos - POS) < 1e-4)) == 1)

% Convert
posV = conn.mapPositionsUnitsToVoxels(pos);
posU = conn.mapPositionsVoxelsToUnits(posV);

% Check the conversion
assert(all(all(abs(posU - POS) < 1e-4)) == 1)

% Try also the different synopses
[posVx, posVy, posVz] = ...
    conn.mapPositionsUnitsToVoxels(pos(:, 1), pos(:, 2), pos(:, 3));
[posUx, posUy, posUz] = ...
    conn.mapPositionsVoxelsToUnits(posVx, posVy, posVz);

% Check the conversion
assert(all(all(abs(posUx - POS(:, 1)) < 1e-4)) == 1)
assert(all(all(abs(posUy - POS(:, 2)) < 1e-4)) == 1)
assert(all(all(abs(posUz - POS(:, 3)) < 1e-4)) == 1)

% Test filtering the selection
% =========================================================================
disp('Test filtering the surpass selection by type...');

% "Select" the spots object
conn.mImarisApplication.SetSurpassSelection(children{4});

% Now get it back, first with the right filter, then with the wrong one
assert(isa(conn.getSurpassSelection('Spots'), class(children{4})));
assert(isempty(conn.getSurpassSelection('Surfaces')));

% Test creating and adding new spots
% =========================================================================
disp('Test creation of new spots...');
vSpotsData = spot.Get;
coords = vSpotsData.mPositionsXYZ + 1.00;
timeIndices = vSpotsData.mIndicesT;
radii = vSpotsData.mRadii;
conn.createAndSetSpots(coords, timeIndices, radii, 'Test', rand(1, 4));
spots = conn.getAllSurpassChildren(0, 'Spots');
assert(numel(spots) == 2);

% Check the filtering and recursion of object finding
% =========================================================================
disp('Get all 7 children with recursion (no filtering)...');
children = conn.getAllSurpassChildren(1);
assert(numel(children) == 7);

disp('Check that there is exactly 1 Light Source...');
children = conn.getAllSurpassChildren(1, 'LightSource');
assert(numel(children) == 1);

disp('Check that there is exactly 1 Frame...');
children = conn.getAllSurpassChildren(1, 'Frame');
assert(numel(children) == 1);

disp('Check that there is exactly 1 Volume...');
children = conn.getAllSurpassChildren(1, 'Volume');
assert(numel(children) == 1);

disp('Check that there are exactly 2 Spots...');
children = conn.getAllSurpassChildren(1, 'Spots');
assert(numel(children) == 2);

disp('Check that there is exactly 1 Surface...');
children = conn.getAllSurpassChildren(1, 'Surfaces');
assert(numel(children) == 1);

disp('Check that there is exactly 1 Measurement Point...');
children = conn.getAllSurpassChildren(1, 'MeasurementPoints');
assert(numel(children) == 1);


% Get the type
% =========================================================================
disp('Get and check the datatype...');
type = conn.getMatlabDatatype();
assert(strcmp(type, 'uint8') == 1);

% Get the data volume
% =========================================================================
disp('Get the data volume...');
stack = conn.getDataVolume(0, 0);

disp('Check the data volume type...');
assert(isa(stack, type) == 1);

% Check the sizes
disp('Check the data volume size...');
assert(size(stack, 1) == DATASETSIZE(1) == 1);
assert(size(stack, 2) == DATASETSIZE(2) == 1);
assert(size(stack, 3) == DATASETSIZE(3) == 1);

% Get the data volume by explicitly passing an iDataSet object
% =========================================================================
disp('Get the data volume by explicitly passing an iDataSet object...');
stack = conn.getDataVolume(0, 0, conn.mImarisApplication.GetDataSet);

disp('Check the data volume type...');
assert(isa(stack, type) == 1);

% Check the sizes
disp('Check the data volume size...');
assert(size(stack, 1) == DATASETSIZE(1) == 1);
assert(size(stack, 2) == DATASETSIZE(2) == 1);
assert(size(stack, 3) == DATASETSIZE(3) == 1);

% Check the getDataVolumeRM() method
% =========================================================================
disp('Check getting the volume in row-major order...');
stackRM = conn.getDataVolumeRM(0, 0);
assert(all(all(stack(:, :, 27) == (stackRM(:, :, 27))')));

% Get the rotation matrix from the camera angle
% =========================================================================
disp('Get the rotation matrix from the camera angle...');
R_D = [
    0.8471    0.2345   -0.4769         0
   -0.1484    0.9661    0.2115         0
    0.5103   -0.1084    0.8532         0
         0         0         0    1.0000
        ];
R = conn.getSurpassCameraRotationMatrix();
assert(all(all(abs(R - R_D) < 1e-4)) == 1);

% Check getting/setting colors and transparency
% =========================================================================
disp('Check getting/setting colors and transparency...');
children = conn.getAllSurpassChildren(1, 'Spots');
spots = children{1};

% We prepare some color/transparency commbinations to circle through (to
% check the int32/uint32 type casting we are forced to apply)
clr = [...
    1 0 0 0.00;    % Red, transparency = 0
    0 1 0 0.00;    % Green, transparency = 0
    0 0 1 0.00;    % Blue,  transparency = 0
    1 1 0 0.00;    % Yellow, transparency = 0
    1 0 1 0.00;    % Purple, transparency = 0
    1 0 1 0.25;    % Purple, transparency = 0.25
    1 0 1 0.50;    % Purple, transparency = 0.50
    1 0 1 0.75;    % Purple, transparency = 0.75
    1 0 1 1.00];  % Purple, transparency = 1.00

for i = 1 : size(clr, 1)
    
    % Set the RGBA color
    spots.SetColorRGBA(conn.mapRgbaVectorToScalar(clr(i, :)));
    
    % Get the RGBA color
    current = conn.mapRgbaScalarToVector(spots.GetColorRGBA());
    
    % Compare (rounding erros allowed)
    assert(abs(all(clr(i, :) - current)) < 1e-2);

end
    
% Close imaris
% =========================================================================
disp('Close Imaris...');
assert(conn.closeImaris(1) == 1)

% Create an ImarisConnector object with starting index 1
% =========================================================================
clear 'conn';
disp('Create an IceImarisConnector object with starting index 1...');
conn = IceImarisConnector([], 1);

% Start Imaris
% =========================================================================
disp('Start Imaris...');
assert(conn.startImaris == 1)

% Check the starting index
% =========================================================================
disp('Check starting index...');
assert(conn.indexingStart() == 1);

% Open a file
% =========================================================================
disp('Load file...');
filename = fullfile(fileparts(which(mfilename)), 'PyramidalCell.ims');
conn.mImarisApplication.FileOpen(filename, '');

% Get and compare the data volume
% =========================================================================
disp('Get and compare the data volume...');
stackIndx1 = conn.getDataVolume(1, 1);

cmp = stack == stackIndx1;
assert(all(cmp(:)));

% Close imaris
% =========================================================================
disp('Close Imaris...');
assert(conn.closeImaris == 1)

% Create an ImarisConnector object with starting index 0
% =========================================================================
clear 'conn';
disp('Create an IceImarisConnector object with starting index 0...');
conn = IceImarisConnector([], 0);

% Start Imaris
% =========================================================================
disp('Start Imaris...');
assert(conn.startImaris == 1)

% Create a dataset
% =========================================================================
disp('Create a dataset');
conn.createDataset('uint8', 100, 200, 50, 3, 10, 0.20, 0.25, 0.5, 0.1);

% Check sizes
% =========================================================================
disp('Check sizes...');
sizes = conn.getSizes();
assert(sizes(1) == 100);
assert(sizes(2) == 200);
assert(sizes(3) == 50);
assert(sizes(4) == 3);
assert(sizes(5) == 10);

% Check voxel sizes
% =========================================================================
disp('Check voxel sizes...');
voxelSizes = conn.getVoxelSizes();
assert(voxelSizes(1) == 0.2);
assert(voxelSizes(2) == 0.25);
assert(voxelSizes(3) == 0.5);

% Check the time delta
% =========================================================================
disp('Check time interval...');
assert(conn.mImarisApplication.GetDataSet().GetTimePointsDelta() == 0.1);

% Check transfering volume data
% =========================================================================
disp('Check two-way data volume transfer...');
data(:, :, 1) = [ 1 2 3; 4 5 6 ];
data(:, :, 2) = [ 7 8 9; 10 11 12];
data = uint8(data);
conn.createDataset('uint8', 3, 2, 2, 1, 1);
conn.setDataVolumeRM(data, 0, 0);
dataOut = conn.getDataVolumeRM(0, 0);
r = data == dataOut;
assert(all(r(:)));

% Close imaris
% =========================================================================
disp('Close Imaris...');
assert(conn.closeImaris(1) == 1)

% All done
% =========================================================================
disp('');
disp('All test succesfully run.');
