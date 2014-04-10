%function trainCellStateClassificationModel_3D

% clc;
% clear;
% close all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data';    
    defaultOutputDir = 'Z:\intravital\data';    

    % there are two ways to specify training data: 
    % (1) a text file with a list of directories containing CellPatternAnnotation.mat files 
    %      (set flagInputRootDirListFile to true for this option)
    % (2) use a GUI to manually select/specify a list of directories
    %      (set flagInputRootDirListFile to false for this option)
    PARAMETERS.flagInputRootDirListFile = false; 
    
    % specify list of models to build
    %PARAMETERS.classificationTask = { 'G1_S_G2_M', 'G1_SG2M', 'S_G2_M', 'G1_S_G2M', 'G2_M', 'G1SG2_M'};
    PARAMETERS.classificationTask = { 'G1_S_G2_M', 'Interphase_Mitotic', 'G1_SG2M' };
    
    % ignore invalid cells?
    PARAMETERS.flagIgnoreBadlyDetectedCells = false;
    
    % ignore badly segmented cells while training?
    PARAMETERS.flagIgnoreBadlySegmentatedCells = false;
    
    % ignore partially visible cells i.e. cells touching image border?
    PARAMETERS.flagIgnoreCellsOnBorder = false;
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PARAMETERS.flagInputRootDirListFile
    
    [fname, pname] = uigetfile( fullfile(defaultDataRootDir, '*.dlist') );
    rootDirListFile = fullfile( pname, fname );
    
    fidDirList = fopen( rootDirListFile );
    rootDirList = {};
    
    curLine = fgetl( fidDirList );
    numLines = 0;
    while ischar(curLine)

        numLines = numLines + 1;
        if isdir(curLine)
           rootDirList{end+1} = curLine; 
        else
           fclose(fidDirList); 
           error( '\nERROR: %s on line %d is not a valid directory. Make sure all the directories in the directory list file exist\n', curLine, numLines);
        end        
        curLine = fgetl(fidDirList);           
        
    end
        
    fclose(fidDirList);
    
else

    % ask the user to provide the list of directories to be processed
    rootDirList = uipickfiles('Prompt', 'Select folders containing the annotation files', ...
                              'REFilter', '(CellPattern).*\.mat$');
    
end

if ~iscell(rootDirList)
    return;
end

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

if ~ischar(outputRootDir)
   return; 
end

hStatusDialog = waitbar(0, 'Building Model ... This might take a while, you may want to get some coffee.');

% make a list of files that need to be processed
annotationFileList = [];

for rid = 1:numel(rootDirList)
    curAnnotationFileList = rdir( fullfile(rootDirList{rid}, '**', 'CellPatternAnnotation.mat') );
    annotationFileList = cat(1, annotationFileList, curAnnotationFileList );
end

if isempty(annotationFileList)        
    msgbox('No annotation files were found in the specified folders', ...
           'Region merging model training', 'error', 'modal' );
    closeStatusDialog(hStatusDialog);
    return;
end

% log
diary off;
diary_file = fullfile( outputRootDir , sprintf( '%s.log' , mfilename ) );
if exist( diary_file, 'file' )
    fclose all;
    delete( diary_file );
end
diary( diary_file );
diary on;

PARAMETERS
rootDirList
outputRootDir

fprintf( '\nList of %d annotation files that will be processed :\n', numel(annotationFileList));
for i = 1:numel(annotationFileList)
    fprintf( '\n%d/%d - %s\n', i, numel(annotationFileList), annotationFileList(i).name );    
end

% process each annotation file
featureData = cell( size(PARAMETERS.classificationTask) );
featureClass = cell( size(PARAMETERS.classificationTask) );
cellInfoData = cell( size(PARAMETERS.classificationTask) );
featureNameList = cell( size(PARAMETERS.classificationTask) );
featureComputationParameters = cell( size(PARAMETERS.classificationTask) );

annotatedCellPatternList = {};
annotatedCellInfoData = {};
stackInfoData = {};

for fid = 1:numel(annotationFileList)
    
   PrettyPrintStepDescription( sprintf( '\n\nComputing Features for file %d/%d ...\n\n', fid, numel(annotationFileList) ) );

   feaCompTimer = tic;

   curAnnotationFilePath = annotationFileList(fid).name   

   % load annotation data
   try
       annotationData = load( curAnnotationFilePath ); 
   catch
      fprintf('\nERROR: could not load annotation file\n' );
      continue;
   end
    
   metadata = annotationData.metadata 
   
   imLabelCellSeg = annotationData.imLabelCellSeg;
   imRegValidMask = annotationData.imRegValidMask;
   
   % pre-processing
   fprintf( '\n>> Pre-processing the image data ... \n' );       
   
   imageDataAdjusted = CellCycleStateClassifier.preprocessImageData( annotationData.imageData );
   
   % make a note of all the cells touching the border
   borderCellIds = setdiff( unique(imLabelCellSeg .* ~imRegValidMask), 0 );
   
   % compute features for each cell
   numCells = numel(annotationData.cellStats);
   
   flagIsValidCell = zeros(numel(annotationData.cellStats),1);
   flagIsCellWellSegmented = zeros(numel(annotationData.cellStats),1);
   flagIsCellOnBorder = zeros(numel(annotationData.cellStats),1);
   
   featureDataLocal = cell( size(PARAMETERS.classificationTask) );
   cellInfoDataLocal = cell( size(PARAMETERS.classificationTask) );
   featureClassLocal = cell( size(PARAMETERS.classificationTask) );
   
   annotatedCellPatternListLocal = cell(numCells, 1);
   
   diary off;
   for cellId = 1:numCells
       
        curCellStats = annotationData.cellStats(cellId);
        
        fprintf( '\n>> Computing Features for cell %d/%d of type - %s ...\n', cellId, numel(annotationData.cellStats), curCellStats.cellPatternType );       
        
        annotatedCellPatternListLocal{cellId} = curCellStats.cellPatternType;

        % compute basic region properties
        cellProps = ComputeRegionProperties( imLabelCellSeg, cellId );
        
        % note down some cell identity information
        curCellInfoStruct.annotationFilePath = curAnnotationFilePath;
        curCellInfoStruct.cellId = cellId;
        curCellInfoStruct.cellPatternType = curCellStats.cellPatternType;
        curCellInfoStruct.Centroid = curCellStats.Centroid;
        curCellInfoStruct.DepthPhysp = curCellStats.Centroid(3) / size(imageDataAdjusted{1},3);
        curCellInfoStruct.VolumePhysp = cellProps.Area * prod(metadata.voxelSpacing);

        [curCellInfoVec , cellInfoLabelList] = ConvertFeatureStructToFeatureVec(curCellInfoStruct);        
        annotatedCellInfoData = cat(1, annotatedCellInfoData, curCellInfoVec);
        
        % ignore badly detected cells 
        flagIsValidCell(cellId) = ~strcmp( curCellStats.cellPatternType, 'Bad_Detection' );
        if PARAMETERS.flagIgnoreBadlyDetectedCells && ~flagIsValidCell(cellId)
            flagIsValidCell(cellId) = false;
            fprintf( '\nCell was erroneously detected and hence will be ignored from training\n' );
            continue;            
        end
        
        % ignore badly segmented cells
        flagIsCellWellSegmented(cellId) = ~ismember( curCellStats.cellPatternType, { 'Under_Segmentation', 'Over_Segmentation'} );
        
        if PARAMETERS.flagIgnoreBadlySegmentatedCells && ~flagIsCellWellSegmented(cellId)
            fprintf( '\nCell was erroneously segmented and hence will be ignored from training\n' );
            continue;            
        end
       
        % check if cell touches the border of valid foreground region
        flagIsCellOnBorder(cellId) = ismember( cellId, borderCellIds );       

        if PARAMETERS.flagIgnoreCellsOnBorder && flagIsCellOnBorder(cellId)
            fprintf( '\nBorder Cell - Part of the cell is outside the valid foreground region. Will be ignored from training.\n' );
            continue;
        end

        % generate features for each classification task
        imCurCellMask = (imLabelCellSeg == cellId);
        
        for tid = 1:numel(PARAMETERS.classificationTask)
            
            fprintf( '\n\tComputing features for task - %s.\n', PARAMETERS.classificationTask{tid} );
            
            % determine feature class
            cellClassFunc = str2func( ['DetermineCellClass_', PARAMETERS.classificationTask{tid}] );
            [curFeatureClass, classNameList] = cellClassFunc( curCellStats.cellPatternType );

            if isempty(curFeatureClass)
                continue;
            end

            featureClassLocal{tid} = cat( 1, featureClassLocal{tid}, {curFeatureClass} ); 

            % compute features
            featurefunc = str2func( ['ComputeCellStateClassificationFeatures_', PARAMETERS.classificationTask{tid}] );

            [featureDataStruct, featureComputationParameters{tid}] = featurefunc(imageDataAdjusted, imRegValidMask, imLabelCellSeg, cellId, metadata.voxelSpacing );          
            [featureVec , featureNameList{tid}] = ConvertFeatureStructToFeatureVec( featureDataStruct );        
            featureDataLocal{tid} = cat(1, featureDataLocal{tid}, featureVec);

            % note down some cell identity information
            cellInfoDataLocal{tid} = cat(1, cellInfoDataLocal{tid}, curCellInfoVec);

        end
        
   end
   diary on;
   
   for tid = 1:numel(PARAMETERS.classificationTask)

       fprintf( '\nclass distribution for classification task - %s: \n', PARAMETERS.classificationTask{tid} );
       tabulate( featureClassLocal{tid} )        
       
       featureData{tid} = cat(1, featureData{tid}, featureDataLocal{tid});       
       featureClass{tid} = cat(1, featureClass{tid}, featureClassLocal{tid});
       cellInfoData{tid} = cat(1, cellInfoData{tid}, cellInfoDataLocal{tid});
       
   end
   
   fprintf( '\nCell_Count: %d\n', numCells );
   curSegmentationAccuracy = 100.0 * mean(flagIsCellWellSegmented & flagIsValidCell);
   fprintf( '\nSegmentation_Accuracy: %.2f%%\n', curSegmentationAccuracy );
   
   fprintf( '\nannotated cell pattern distribution: \n' );
   tabulate( annotatedCellPatternListLocal )   
   annotatedCellPatternList = cat(1, annotatedCellPatternList, annotatedCellPatternListLocal );
   
   computationTime = toc(feaCompTimer)   
    
   waitbar(fid / numel(annotationFileList), hStatusDialog);
   
end

fprintf( '\n\n*****************************************************\n\n' );

if isempty(annotatedCellPatternList)        
    msgbox('No valid annotation files were found in the specified folders', ...
           'Region merging model training', 'error', 'modal' );
    closeStatusDialog(hStatusDialog);
    return;
end

fprintf( '\nannotated cell pattern distribution in all datasets: \n' );
tabulate( annotatedCellPatternList )   

classFeatureName = 'ClassLabel';
for tid = 1:numel(PARAMETERS.classificationTask)
    
    PrettyPrintStepDescription( sprintf( 'Generating model files for classification task - %s', PARAMETERS.classificationTask{tid} ) );    

    fprintf( '\nTotal class distribution: \n' );
    tabulate( featureClass{tid} );
    
    curFeatureComputationParameters = featureComputationParameters{tid}
    
    featureFileNamePrefix = sprintf( 'cellStateFeatures_%s', PARAMETERS.classificationTask{tid} );
    
    % write csv file -- without debug info
    fprintf( '\nWriting csv file without cell info ... \n' );
    
    featureLabels = cat(2, classFeatureName, featureNameList{tid});
    featureMatrix = cat(2, featureClass{tid}, featureData{tid});
    
    WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, [featureFileNamePrefix '.csv']), featureMatrix, featureLabels );

    % write csv file -- with debugging info
    fprintf( '\nWriting csv file with cell info ... \n' );

    featureLabels = cat(2, classFeatureName, cellInfoLabelList, featureNameList{tid});
    featureMatrix = cat(2, featureClass{tid}, cellInfoData{tid}, featureData{tid});
    
    WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, [featureFileNamePrefix '_info.csv']), featureMatrix, featureLabels );

    % write arff file for weka
    fprintf( '\nConverting csv file to arff file ... \n' );

    AddWekaClassesToPath();

    ConvertCSVFileToArffFile( fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );

    % generate classification model
    wekaTrainingDataset = getWekaDatasetFromCSVFile( fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );
    wekaTrainingDataset.setClassIndex(0);
    wekaModel = BuildWekaModelForCellCycleStateIdentification(wekaTrainingDataset, ...
                                                              fullfile(outputRootDir, [PARAMETERS.classificationTask{tid}, '.model']));
    
end

% switch off diary
diary off;
clostStatusDialog(hStatusDialog);

%end
