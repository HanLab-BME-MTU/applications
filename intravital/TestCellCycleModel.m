

clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis\annotations\cell_pattern_annotations';    
    defaultOutputDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis';    
    defaultModelDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis\cell_cycle_model';
    
    flagIgnoreCellsOnBorder = true;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ask the user to provide the list of directories to be processed
rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% ask the user to select the model directory
modelDir = uigetdir( defaultModelDir, 'Select Model Output Directory' ); 

% load cell cycle model
cls = CellCycleClassifier();
cls.LoadModelFromDisk( modelDir );
labelList = { 'G1', 'S', 'G2', 'M' };

for rid = 1:numel(rootDirList)

    fprintf( '\n\n*****************************************************\n\n' );           
    
    % get a directory to process from the list
    dataRootDir = rootDirList{rid};     

    fprintf( '\nGenerating features for root directory - %d/%d\n', rid, numel(rootDirList) );

    % set output directory accordingly
    [pathstr, name, ext] = fileparts( dataRootDir );
    outputDir = fullfile(outputRootDir, name);
    
    if ~isdir(outputDir)
        mkdir(outputDir);
    end

    dataRootDir
    outputDir
    
    % log
    diary off;
    diary_file = fullfile( outputDir , sprintf( '%s.log' , mfilename ) );
    if exist( diary_file, 'file' )
        fclose all;
        delete( diary_file );
    end
    diary( diary_file );
    diary on;

    % get a list of all cell pattern annotation files
    annotationFileList = rdir( fullfile(dataRootDir, '**/CellPatternAnnotation.mat') );

    % run prediction on each annotation file
    globalTrueLabels = {};
    globalPredictedLabels = {};
    globalPerformanceData = [];
    
    for fid = 1:numel(annotationFileList)

       fprintf( '\n\nRunning Cell Cycle Classifier for file %d/%d ...\n\n', fid, numel(annotationFileList) );       

       tic

       curAnnotationFilePath = annotationFileList(fid).name   
       curAnnotationRelativeFilePath = annotationFileList(fid).name((numel(dataRootDir)+2):end);

       [pathstr, name, ext] = fileparts(curAnnotationRelativeFilePath); 
       curOutputDir = fullfile( outputDir, pathstr );
       if ~isdir(curOutputDir)
           mkdir(curOutputDir);
       end

       % load annotation data
       try
            annotationData = load( curAnnotationFilePath ); 
       catch
           fprintf('\nERROR: could not load annotation file\n' );
           continue;
       end

       % run cell cycle classifier
       [ predictedLabels ] = cls.predict( annotationData.imageData, ...
                                          annotationData.imLabelCellSeg, ...
                                          annotationData.cellStats);
       
       
       % make a note of all the cells touching the border
       borderCellIds = setdiff( unique(annotationData.imLabelCellSeg .* ~annotationData.imRegValidMask), 0 );

       % get true labels (reduce to G1, S, G2, M)
       flagIsValidCell = zeros(numel(annotationData.cellStats), 1);
       trueLabels = cell( size(predictedLabels) );
       
       for cid = 1:numel(annotationData.cellStats)
                      
           curCellStats = annotationData.cellStats(cid);
           
           flagIsValidCell(cid) = true;
           
           switch(curCellStats.cellPatternType)

                case { 'Mono_Pre_G1', 'Multi_Pre_G1', 'Mono_G1', 'Multi_G1' } 

                    curTrueLabel = 'G1';

                case { 'Mono_S', 'Multi_S' } 

                    curTrueLabel = 'S';

                case { 'Mono_G2', 'Multi_G2' }     

                    curTrueLabel = 'G2';
                    
                case { 'Mono_Prophase', 'Multi_Prophase' }     

                    curTrueLabel = 'M';                    

                otherwise

                    curTrueLabel = 'Other';
                    flagIsValidCell(cid) = false;
           end             
        
           if ismember( cid, borderCellIds )
                curTrueLabel = 'Border_Cell';
                flagIsValidCell(cid) = false;
           end
           
           trueLabels{cid} = curTrueLabel;
           
       end

       curPerformanceData.annotationFile = curAnnotationFilePath;
       curPerformanceData.flagIsValidCell = flagIsValidCell;
       curPerformanceData.trueLabels = trueLabels;
       curPerformanceData.predictedLabels = predictedLabels;
       curPerformanceData.modelDir = modelDir;
       
       % compute performance stats
       predictedLabels = predictedLabels(flagIsValidCell > 0);
       trueLabels = trueLabels(flagIsValidCell > 0);

       curPerformanceStats = ComputeClassificationPerformance_MultiClass( predictedLabels, trueLabels, labelList )
       
       confusioinMatrix_Display = curPerformanceStats.cMatrix_Display
       
       curPerformanceData.stats = curPerformanceStats;
       curPerformanceData.labelList = labelList;
       
       % write performance data to a mat file
       performanceData = curPerformanceData;
       save( fullfile(curOutputDir, 'cellCycleModelPerformance.mat'), ...
             'performanceData' );
       
       % append to global
       globalPredictedLabels = cat(1, globalPredictedLabels, predictedLabels);
       globalTrueLabels = cat(1, globalTrueLabels, trueLabels);
       globalPerformanceData = cat(1, globalPerformanceData, curPerformanceData);
       
    end
    
    % compute global performance stats
    globalPerformanceStats = ComputeClassificationPerformance_MultiClass( globalPredictedLabels, globalTrueLabels, labelList )
    
    performanceData = [];    
    performanceData.dataRootDir = dataRootDir;
    performanceData.annotationFileList = annotationFileList;
    performanceData.trueLabels = globalTrueLabels;
    performanceData.predictedLabels = globalPredictedLabels;
    performanceData.labelList = labelList;
    performanceData.stats = globalPerformanceStats;
    performanceData.subPerformanceData = globalPerformanceData;
    
   % write performance data to a mat file
   save( fullfile(outputDir, 'cellCycleModelPerformance.mat'), ...
         'performanceData' );
    
   % switch off diary
   diary off;
    
end    
