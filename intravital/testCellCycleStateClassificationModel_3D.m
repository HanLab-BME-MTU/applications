
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data';    
    defaultOutputDir = 'Z:\intravital\data';    
    defaultModelDir = 'Z:\intravital\data';    
    
    PARAMETERS.flagInputRootDirListFile = true;
    
    PARAMETERS.flagIgnoreBadlySegmentatedCells = false;
    PARAMETERS.flagIgnoreCellsOnBorder = false;
    
    PARAMETERS.flagSaveImages = true;
    PARAMETERS.flagSaveErrorImagesOnly = false;
    
    % classification model: options - OneStage, TwoStage_G1_SG2M, TwoStage_G1_S_G2M
    PARAMETERS.classificationModelType = 'OneStage';
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ask the user to provide the list of directories to be processed
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
    rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);
    
end

% ask the user to select the model file
modelPath = uigetdir( defaultModelDir, 'Select Weka Model Dir' ); 

% load classification model
switch PARAMETERS.classificationModelType
    
    case 'OneStage'

        modelFilePath = fullfile(modelPath, 'G1_S_G2_M.model');
        model = CellCycleStateClassifier_OneStage( modelFilePath );
        
    case 'TwoStage_G1_SG2M'
    
        modelFilePath_Stage1 = fullfile(modelPath, 'G1_SG2M.model');
        modelFilePath_Stage2 = fullfile(modelPath, 'S_G2_M.model');
        model = CellCycleStateClassifier_TwoStage_G1_SG2M(modelFilePath_Stage1, modelFilePath_Stage2);
        
    case 'TwoStage_G1_S_G2M'

        modelFilePath_Stage1 = fullfile(modelPath, 'G1_S_G2M.model');
        modelFilePath_Stage2 = fullfile(modelPath, 'G2_M.model');
        model = CellCycleStateClassifier_TwoStage_G1_S_G2M(modelFilePath_Stage1, modelFilePath_Stage2);

    case 'Interphase_Mitotic'
        
        modelFilePath = fullfile(modelPath, 'Interphase_Mitotic.model');
        model = CellCycleStateClassifier_OneStage( modelFilePath );
        
    otherwise 
        
        error( 'Unrecognized Classification Model Type' );
        
end

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% make a list of files that need to be processed
annotationFileList = [];

for rid = 1:numel(rootDirList)
    curAnnotationFileList = rdir( fullfile(rootDirList{rid}, '**', 'CellPatternAnnotation.mat') );
    annotationFileList = cat(1, annotationFileList, curAnnotationFileList );
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
modelPath

if PARAMETERS.flagSaveImages
    if isdir( fullfile(outputRootDir, 'images') )            
        rmdir( fullfile(outputRootDir, 'images'), 's' );
    end
end

fprintf( '\nList of %d annotation files that will be processed :\n', numel(annotationFileList));
for i = 1:numel(annotationFileList)
    fprintf( '\n%d/%d - %s\n', i, numel(annotationFileList), annotationFileList(i).name );    
end

% process each annotation file
cellInfoData = {};
annotatedFeatureClass = {};
predictedFeatureClass = {};
performanceDataGlobal = [];

annotatedCellPatternList = {};
modelClassNameList = model.getClassNameList()

for fid = 1:numel(annotationFileList)
    
   fprintf( '\n\n*****************************************************\n\n' );           
   
   fprintf( '\n\nRunning cell inlier detection model for file %d/%d ...\n\n', fid, numel(annotationFileList) );       
    
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
   
   dataFilePath = annotationData.dataFilePath;
   imLabelCellSeg = annotationData.imLabelCellSeg;
   imRegValidMask = annotationData.imRegValidMask;
   
   [pname, nuclearMarkerFileName, fext] = fileparts( dataFilePath{1} );
   [pname, fucciFileName, fext] = fileparts( dataFilePath{2} );
   
   % pre-processing
   fprintf( '\n>> Pre-processing the image data ... \n' );       

   imageDataAdjusted = CellCycleStateClassifier.preprocessImageData( annotationData.imageData );
   
   % make a note of all the cells touching the border
   borderCellIds = setdiff( unique(imLabelCellSeg .* ~imRegValidMask), 0 );
   
   % compute features for each cell
   numCells = numel(annotationData.cellStats);
   
   flagIsCellWellSegmented = zeros(numel(annotationData.cellStats),1);
   flagIsCellOnBorder = zeros(numel(annotationData.cellStats),1);
   
   cellInfoDataLocal = {};
   annotatedFeatureClassLocal = {};
   predictedFeatureClassLocal = {};
   
   annotatedCellPatternListLocal = cell(numCells, 1);
   
   diary off;
   
   for cellId = 1:numCells
       
        curCellStats = annotationData.cellStats(cellId);
        imCurCellMask = (imLabelCellSeg == cellId);
        
        fprintf( '\n>> Applying classification model on cell %d/%d of type - %s ...\n', cellId, numel(annotationData.cellStats), curCellStats.cellPatternType );       
        
        annotatedCellPatternListLocal{cellId} = curCellStats.cellPatternType;
        
        % ignore badly segmented cells
        flagIsCellWellSegmented(cellId) = ~ismember( curCellStats.cellPatternType, { 'Under_Segmented', 'Over_Segmented'} );
        
        if PARAMETERS.flagIgnoreBadlySegmentatedCells && ~flagIsCellWellSegmented(cellId)
            fprintf( '\n\tCell was erroneously segmented and hence will be ignored from training\n' );
            continue;            
        end
       
        % check if cell touches the image border
        flagIsCellOnBorder(cellId) = ismember( cellId, borderCellIds );       

        if PARAMETERS.flagIgnoreCellsOnBorder && flagIsCellOnBorder(cellId)
            fprintf( '\n\tBorder Cell - Part of the cell is outside the valid foreground region. Will be ignored from training.\n' );
            continue;
        end

        % determine feature class
        [curAnnotatedFeatureClass, classNameList] = DetermineCellClass_G1_S_G2_M( curCellStats.cellPatternType );
        
        if isempty(curAnnotatedFeatureClass)
            fprintf('\n\tInvalid cell pattern for current classification task. Will be ignored.\n');
            continue;
        end
        
        annotatedFeatureClassLocal = cat( 1, annotatedFeatureClassLocal, {curAnnotatedFeatureClass} );
        
        % apply classifier
        [curPredictedFeatureClass, ...
         curClassPredictionProbabilities, ~] = model.predictCell(imageDataAdjusted, annotationData.imRegValidMask, imLabelCellSeg, cellId, metadata.voxelSpacing); 
        
        predictedFeatureClassLocal = cat( 1, predictedFeatureClassLocal, {curPredictedFeatureClass} );
        
        fprintf( '\n\tTrue Label - %s, Predicted Label = %s\n', curAnnotatedFeatureClass, curPredictedFeatureClass);
        
        % note down some auxillary information
        curCellInfoStruct.predictionProbability = max(curClassPredictionProbabilities);
        
        for classId = 1:numel(modelClassNameList)
            curFieldName = ['predictionProbability_', modelClassNameList{classId}];
            curCellInfoStruct.(curFieldName) = curClassPredictionProbabilities(classId);
        end
        
        % note down some cell identity information
        curCellInfoStruct.annotationFilePath = curAnnotationFilePath;
        curCellInfoStruct.nuclearMarkerFileName = nuclearMarkerFileName; 
        curCellInfoStruct.fucciFileName = fucciFileName; 
        curCellInfoStruct.cellId = cellId;
        curCellInfoStruct.cellPatternType = curCellStats.cellPatternType;
        curCellInfoStruct.Centroid = curCellStats.Centroid;
        curCellInfoStruct.Depth = curCellStats.Centroid(3) / size(annotationData.imageData{1},3);
        
        % save images if requested
        if PARAMETERS.flagSaveImages 
            
            if PARAMETERS.flagSaveErrorImagesOnly && strcmp(curAnnotatedFeatureClass, curPredictedFeatureClass)
                continue; % dont save images if the classifier prediction was correct
            end
            
            curImageRelDir = fullfile( 'images', [curAnnotatedFeatureClass, '_' curPredictedFeatureClass], sprintf('fid_%.2d_%s', fid, strtrim(nuclearMarkerFileName)) );
            curImageOutputDir = fullfile(outputRootDir, curImageRelDir );

            if ~isdir( curImageOutputDir )
                mkdir( curImageOutputDir );
            end
            
            curCellStats = annotationData.cellStats(cellId);
            curCellCentroid = curCellStats.Centroid;
            curCellBoundingBox = curCellStats.BoundingBox;
            curCellDisplaySize = max( [curCellBoundingBox(4:5), 70] );
            szOutputImage = [100, 100];
           
            % crop images
            subinds = cell(1,3);
            imsize = size(annotationData.imageData{1});
            for i = 1:2

                xi = round(curCellCentroid(3-i) - 0.5 * curCellDisplaySize);

                xi_low = xi;
                if xi_low < 1 
                    xi_low = 1;
                end

                xi_high = xi + curCellDisplaySize - 1;
                if xi_high > imsize(i)
                    xi_high = imsize(i);
                end

                subinds{i} = xi_low:xi_high;

            end     
            subinds{3} = round(curCellCentroid(3));
            
            imCurCellMidSliceAllChannel = [];
            for j = 1:3
                channelDisplayRange = ComputeImageDynamicRange( annotationData.imageData{j}, 98.0 );
                imCurChannelCropped = mat2gray( annotationData.imageData{j}( subinds{:} ), channelDisplayRange );
                imCurChannelCropped = imresize( imCurChannelCropped, szOutputImage );
                imCurCellMidSliceAllChannel = cat(3, imCurCellMidSliceAllChannel, imCurChannelCropped );
            end
            
            imCurCellSegBndMidSliceCropped = imresize( bwperim( imCurCellMask( subinds{:} ) ), szOutputImage, 'nearest' );
            
            imCurCellMIP = mat2gray( max( annotationData.imageData{1}(subinds{1:2}, :) .* imCurCellMask(subinds{1:2}, :), [], 3) );
            imCurCellMIP = imresize(imCurCellMIP, szOutputImage);
            
            % write images
            histoneImageFileName = sprintf('CellMidSliceHistone_%.3d.png', cellId);
            imwrite( genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 ), ...
                     fullfile(curImageOutputDir, histoneImageFileName), 'png' );   
            
            channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
            
            fucciImageFileName = sprintf('CellMidSliceFUCCI_%.3d.png', cellId);
            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
                     fullfile(curImageOutputDir, fucciImageFileName), 'png' );   

            allChannelImageFileName = sprintf('CellMidSliceAllChannel_%.3d.png', cellId);
            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
                     fullfile(curImageOutputDir, allChannelImageFileName), 'png' );   
                 
            mipImageFileName = sprintf('CellMIP_%.3d.png', cellId);
            imwrite(imCurCellMIP, fullfile(curImageOutputDir, mipImageFileName), 'png' );   
           
            curCellInfoStruct.imageOutputRelDir = curImageRelDir;
            curCellInfoStruct.histoneImageRelPath = fullfile(curImageRelDir, histoneImageFileName);
            curCellInfoStruct.fucciImageRelPath = fullfile(curImageRelDir, fucciImageFileName);
            curCellInfoStruct.allChannelImageRelPath = fullfile(curImageRelDir, allChannelImageFileName);
            curCellInfoStruct.mipImageRelPath = fullfile(curImageRelDir, mipImageFileName);
            
        end
        
        [curCellInfoVec, cellInfoLabelList] = ConvertFeatureStructToFeatureVec(curCellInfoStruct);        
        cellInfoDataLocal = cat(1, cellInfoDataLocal, curCellInfoVec);
        
   end
   
   diary on;
   
   fprintf( '\nannotated cell pattern distribution: \n' );
   tabulate( annotatedCellPatternListLocal )   
   annotatedCellPatternList = cat(1, annotatedCellPatternList, annotatedCellPatternListLocal );
   
   fprintf( '\nannotated class distribution:\n' );
   tabulate( annotatedFeatureClassLocal ) 
   
   fprintf( '\npredicted class distribution:\n' );
   tabulate( predictedFeatureClassLocal ) 
   
   % compute performance stats
   curPerformanceStats = ComputeClassificationPerformance_MultiClass(predictedFeatureClassLocal, annotatedFeatureClassLocal, classNameList)
   
   confusioinMatrix_Display = curPerformanceStats.cMatrix_Display
   
   curPerformanceData.annotationFilePath = curAnnotationFilePath;
   curPerformanceData.trueLabels = annotatedFeatureClassLocal;
   curPerformanceData.predictedLabels = predictedFeatureClassLocal;
   curPerformanceData.stats = curPerformanceStats;
   curPerformanceData.cellInfoData = cellInfoDataLocal;
   
   % append to global
   cellInfoData = cat(1, cellInfoData, cellInfoDataLocal);
   annotatedFeatureClass = cat(1, annotatedFeatureClass, annotatedFeatureClassLocal);
   predictedFeatureClass = cat(1, predictedFeatureClass, predictedFeatureClassLocal);
   performanceDataGlobal = cat(1, performanceDataGlobal, curPerformanceData);
   
   computationTime = toc(feaCompTimer)   
    
end

fprintf( '\n\n*****************************************************\n\n' );

fprintf( '\nAnnotated cell pattern distribution in all datasets: \n' );
tabulate( annotatedCellPatternList )   

fprintf( '\nAnnotated class distribution accross all datasets: \n' );
tabulate( annotatedFeatureClass );

fprintf( '\nPredicted class distribution accross all datasets: \n' );
tabulate( predictedFeatureClass );

% compute global performance stats
fprintf( '\nComputing overall performance stats: \n' );

globalPerformanceStats = ComputeClassificationPerformance_MultiClass( predictedFeatureClass, annotatedFeatureClass, classNameList )
globalConfusioinMatrix_Display = globalPerformanceStats.cMatrix_Display

performanceData = [];
performanceData.annotationFileList = annotationFileList;
performanceData.rootDirList = rootDirList;
performanceData.modelPath = modelPath;
performanceData.classNameList = classNameList;
performanceData.trueLabels = annotatedFeatureClass;
performanceData.predictedLabels = predictedFeatureClass;
performanceData.stats = globalPerformanceStats;
performanceData.cellInfoData = cellInfoData;
performanceData.subPerformanceData = performanceDataGlobal;

save( fullfile(outputRootDir, 'cellCycleStateClassificationPerformance.mat'), 'performanceData' );

% write global files
fprintf( '\n>>Writing csv file for inspection ... \n' );

featureFileNamePrefix = 'cellCycleStateClassificationPerformance';
featureLabels = cat(2, 'TrueClassLabel', 'PredictedClassLabel', cellInfoLabelList);
featureMatrix = cat(2, annotatedFeatureClass, predictedFeatureClass, cellInfoData);

WriteFeatureMatrixToCSVFile( fullfile(pwd, [featureFileNamePrefix '.csv']), featureMatrix, featureLabels );
movefile( fullfile(pwd, [featureFileNamePrefix '.csv']), fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );

% switch off diary
diary off;
