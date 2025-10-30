% setupOneMDFromMultiNDs.m prepares (specific to Alex Rivera's data) one MD
% file based on multiple ND2 files taken separately for bead and cell

% Run this script per folder
%% Get the directory
curDir = dir(pwd);
%% Get the bead and cell channgels
names = arrayfun(@(x) x.name, curDir,'unif',false);
priorBead = 'bead_t';
priorCell = 'cell_t';

% Filter the names to get only those starting with 'bead_t'
beadFiles = names(startsWith(names, priorBead));
cellFiles = names(startsWith(names, priorCell));

%% Sort the bead files based on the numeric part of the file name
% Extract the numeric part from the file names
beadNumbers = cellfun(@(x) sscanf(x, 'bead_t%d'), beadFiles);
cellNumbers = cellfun(@(x) sscanf(x, 'cell_t%d'), cellFiles);

% Sort the bead files based on the numeric part
[~, sortIdx] = sort(beadNumbers);
sortedBeadFiles = beadFiles(sortIdx);
[~, sortIdxCell] = sort(cellNumbers);
sortedCellFiles = cellFiles(sortIdxCell);

%% Display the sorted bead file names
disp('Sorted Bead files:');
disp(sortedBeadFiles);
disp('Sorted Cell files:');
disp(sortedCellFiles);
%% Creat folders and save tif files per channel
%% Add Bio-Formats package to MATLAB path (adjust the path accordingly)
% Initialize Bio-Formats reader
bfInitLogging('INFO');

%% Create bead and cell folders if they don't exist
if ~exist('Bead', 'dir')
    mkdir('Bead');
end
if ~exist('Cell', 'dir')
    mkdir('Cell');
end

%% Process bead files
for i = 1:length(sortedBeadFiles)
    save_nd2_to_tif(sortedBeadFiles{i}, 'Bead');
end

%% Process cell files
for i = 1:length(sortedCellFiles)
    save_nd2_to_tif(sortedCellFiles{i}, 'Cell');
end

%% Done!
disp('ND2 files have been converted and saved as tif files in the respective folders.');

%% Get the metadata from bead channel
% Get the reader and retrieve the metadata
reader = bfGetReader(sortedBeadFiles{1});

% Get metadata store object
omeMeta = reader.getMetadataStore();

%% Extract some useful metadata
% Get the number of series (image sets) in the file
% Set the current series
reader.setSeries(0);
series = 1;
% Extract relevant metadata for the current series
pixelSizeXObj = omeMeta.getPixelsPhysicalSizeX(series - 1); 
NA = omeMeta.getObjectiveLensNA(series - 1, 0);
bitDepthObj = omeMeta.getPixelsSignificantBits(series - 1); 
if isempty(NA)
    NA = 1.2;
end
reader.close();
%% Registering each folder as channel
beadChan = Channel([pwd filesep 'Bead']);
beadChan.emissionWavelength_= 647;
beadChan.imageType_ = 'Confocal';
cellChan = Channel([pwd filesep 'Cell']);

integChans=[beadChan cellChan];
MD = MovieData(integChans);
MD.timeInterval_ = 5; % assumed
MD.pixelSize_ = round(double(pixelSizeXObj.value)*1000);
MD.numAperture_ = NA;
MD.camBitdepth_ = double(bitDepthObj.getValue());
MD.setPath(pwd)
MD.outputDirectory_ = pwd;
MD.setFilename('movieData.mat')
MD.sanityCheck;
MD.save
%% Save ref file
refFile = 'reference_01.nd2';
save_nd2_to_tif(refFile,pwd)
%% Function to process and save nd2 files to tif
function save_nd2_to_tif(nd2FileName, folderName)
    % Create reader for the nd2 file
    reader = bfGetReader(nd2FileName);
    
    reader.setSeries(0);  % Set the current series
    img = bfGetPlane(reader, 1);   % Read the first image plane
    
    % Create a filename for the tif
    tifFileName = fullfile(folderName, [nd2FileName(1:end-4), '.tif']); %assuming there is only one frame per ND2 file
    
    % Save the image as a tif file
    imwrite(img, tifFileName);

    % Close the reader after processing
    reader.close();
end