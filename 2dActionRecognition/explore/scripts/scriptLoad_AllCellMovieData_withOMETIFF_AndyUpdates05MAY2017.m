loadMatPath = '/work/bioinformatics/shared/dope/export/Andres/cellMovieData_wAnno05May2017_0927.mat'

% My output path
omeTiffDir = '/work/bioinformatics/shared/dope/data/OMETIFF/all/';

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Updating with Andy''s annotations');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

load(loadMatPath);%, 'allCellsMovieData'); 

allCellsSet = cellDataSet;

masterMovieDir = [omeTiffDir filesep 'sdcAll' filesep];
if ~exist(masterMovieDir, 'dir')
    mkdir(masterMovieDir);
end

parfor i = 1:length(allCellsSet)

    iCell = allCellsSet{i};
    oldMD = MovieData.loadMatFile(iCell.cellMD);
    disp(['Creating movieData for cell ID key: ' iCell.key]);
    movieOMETIFF = oldMD.channels_.channelPath_;
    disp(['Making MD for OME-TIFF : ' movieOMETIFF]);

    disp(iCell.key);

    ts = iCell.ts;
    
    startTime = ts;
    ntime = length(ts);
    endTime = ts(1) + ntime;

    iCell.oldkey = iCell.key;
    iCell.key = [iCell.key '_t' num2str(endTime)]

    
    % MovieData Validation
    disp('Loading with MovieData BF reader')
    
    NewMDOutDir = [masterMovieDir filesep iCell.key filesep];
    mkdir(NewMDOutDir);

    MD = MovieData(movieOMETIFF, true, NewMDOutDir);
    MD.sanityCheck();
    MD.save();
    
    extProc = ExternalProcess(MD, 'LBPfeatures');
    MD.addProcess(extProc);
    MD.processes_{1}.setParameters(iCell);

    MD.addProcess(EfficientSubpixelRegistrationProcess(MD));
    SDCindx = MD.getProcessIndex('EfficientSubpixelRegistrationProcess');
    MD.processes_{SDCindx}.run()

    % MDpath = MD.getFullPath();
    % Save individual cell movieData
    allCellsSet{i}.cellMD = MD.getFullPath();
    allCellsSet{i}.key = iCell.key;
    allCellsSet{i}.oldkey = iCell.oldkey;
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Save cell array with cell MD info');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

save([omeTiffDir filesep 'clickFuryCellData.mat'], 'allCellsSet');

% annotationSet = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells (by local master index)
% annotationSet('bleb')=[NaN]
% annotationSet('protrusion')=[NaN]
% annotationSet('small')=[NaN]
% annotationSet('big')=[NaN]
% annotationSet('active')=[NaN]
% annotationSet('inactive')=[NaN]
% annotationSet('weird')=[NaN]
% annotationSet('neat')=[NaN]
annotationSet.keys
cellDataSet = allCellsSet;

save([masterMovieDir filesep 'AndyCellData05MAY2017.mat'], 'cellDataSet', 'annotationSet');



% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Examples to acccess metadata via MovieData/MovieList');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% Access metadata info in the following fashion...
% ML.movies_{3}.reader.formatReader.getMetadataStore().getDatasetID(0)
% ML.movies_{3}.reader.formatReader.getMetadataStore().getDatasetName(0) 
% ML.movies_{3}.reader.formatReader.getMetadataStore().getImageDescription(0) 
% ML.movies_{3}.reader.formatReader.getMetadataStore().getExperimenterGroupDescription(0)
% ML.movies_{3}.reader.formatReader.getMetadataStore().getsetDatasetDescription(0)

% annotationSet = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells (by local master index)
% annotationSet('bleb')=[NaN]
% annotationSet('protrusion')=[NaN]
% annotationSet('small')=[NaN]
% annotationSet('big')=[NaN]
% annotationSet('active')=[NaN]
% annotationSet('inactive')=[NaN]
% annotationSet('weird')=[NaN]
% annotationSet('neat')=[NaN]
% annotationSet.keys
% cellDataSet = allCellsSet;

