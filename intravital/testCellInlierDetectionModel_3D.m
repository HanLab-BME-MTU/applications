clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data';    
    defaultOutputDir = 'Z:\intravital\data';    
    defaultModelDir = 'Z:\intravital\data';    
    
    PARAMETERS.flagIgnoreBadlySegmentatedCells = false;
    PARAMETERS.flagIgnoreCellsOnBorder = false;
    PARAMETERS.flagIgnoreApoptoticCells = false;
    PARAMETERS.flagConsiderApoptoticCellsAsOutliers = true;
    
    PARAMETERS.flagSaveImages = true;
    PARAMETERS.flagSaveErrorImagesOnly = true;
    
    PARAMETERS.classificationTaskType = 'TwoClass';
    
    strOutlierClassName = 'Outlier';
    strInlierClassName = 'Inlier';
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ask the user to provide the list of directories to be processed
rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);

% ask the user to select the model file
[modelFileName, modelDir] = uigetfile( fullfile(defaultModelDir, '*.model'), 'Select Model Output Directory' ); 
modelFilePath = fullfile(modelDir, modelFileName);

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% load classification model
AddWekaClassesToPath();

switch PARAMETERS.classificationTaskType
    
    case 'OneClass'
        
        model = CellInlierDetectionClassifier_OneClassSVM(modelFilePath);
        
    case 'TwoClass'

        model = CellInlierDetectionClassifier_TwoClass(modelFilePath);
        
end

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
outputRootDir
modelFilePath

if PARAMETERS.flagSaveImages
    if isdir( fullfile(outputRootDir, 'images') )            
        rmdir( fullfile(outputRootDir, 'images'), 's' );
    end
end

% process each annotation file
cellInfoData = {};
annotatedFeatureClass = {};
predictedFeatureClass = {};
performanceDataGlobal = [];

classNameList =  { strInlierClassName, strOutlierClassName };
   

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
   
   imLabelCellSeg = annotationData.imLabelCellSeg;
   imRegValidMask = annotationData.imRegValidMask;
   
   % pre-processing
   fprintf( '\n>> Pre-processing the image data ... \n' );       

   imageDataAdjusted = CellInlierDetectionClassifier.preprocessImageData( annotationData.imageData{1} );
   
   % make a note of all the cells touching the border
   borderCellIds = setdiff( unique(imLabelCellSeg .* ~imRegValidMask), 0 );
   
   % compute features for each cell
   numCells = numel(annotationData.cellStats);
   
   flagIsCellWellSegmented = zeros(numel(annotationData.cellStats),1);
   flagIsCellOnBorder = zeros(numel(annotationData.cellStats),1);
   
   cellInfoDataLocal = {};
   annotatedFeatureClassLocal = {};
   predictedFeatureClassLocal = {};
   
   for cellId = 1:numCells
       
        curCellStats = annotationData.cellStats(cellId);
        imCurCellMask = (imLabelCellSeg == cellId);
        
        fprintf( '\n>> Computing Features for cell %d/%d of type - %s ...\n', cellId, numel(annotationData.cellStats), curCellStats.cellPatternType );       
        
        % ignore badly segmented cells
        flagIsCellWellSegmented(cellId) = ~ismember( curCellStats.cellPatternType, { 'Under_Segmented', 'Over_Segmented'} );
        
        if PARAMETERS.flagIgnoreBadlySegmentatedCells && ~flagIsCellWellSegmented(cellId)
            fprintf( '\nCell was erroneously segmented and hence will be ignored from training\n' );
            continue;            
        end
       
        % check if cell touches the image border
        flagIsCellOnBorder(cellId) = ismember( cellId, borderCellIds );       

        if PARAMETERS.flagIgnoreCellsOnBorder && flagIsCellOnBorder(cellId)
            fprintf( '\nBorder Cell - Part of the cell is outside the valid foreground region. Will be ignored from training.\n' );
            continue;
        end

        % determine feature class
        switch curCellStats.cellPatternType
        
            case { 'Bad_Detection' }
                
                curAnnotatedFeatureClass = strOutlierClassName;
                
            case { 'Apoptotic_Green', 'Apoptotic_Red', 'Apoptotic_Black' }
            
                if PARAMETERS.flagIgnoreApoptoticCells
                    continue;
                else
                    if PARAMETERS.flagConsiderApoptoticCellsAsOutliers
                        curAnnotatedFeatureClass = strOutlierClassName;
                    else
                        curAnnotatedFeatureClass = strInlierClassName;
                    end
                end
               
            otherwise
                
                curAnnotatedFeatureClass = strInlierClassName;
            
        end
        
        annotatedFeatureClassLocal = cat( 1, annotatedFeatureClassLocal, {curAnnotatedFeatureClass} );
        
        % apply classifier
        predictedLabel = model.predictCell(imageDataAdjusted, imLabelCellSeg, cellId, metadata.voxelSpacing); 
        
        if predictedLabel
            curPredictedFeatureClass = strInlierClassName;
        else
            curPredictedFeatureClass = strOutlierClassName;
        end
        
        predictedFeatureClassLocal = cat( 1, predictedFeatureClassLocal, {curPredictedFeatureClass} );
        
        fprintf( '\n\tTrue Label - %s, Predicted Label = %s\n', curAnnotatedFeatureClass, curPredictedFeatureClass);
        
        % note down some cell identity information
        curCellInfoStruct.annotationFilePath = curAnnotationFilePath;
        curCellInfoStruct.cellId = cellId;
        curCellInfoStruct.cellPatternType = curCellStats.cellPatternType;
        curCellInfoStruct.Centroid = curCellStats.Centroid;
        curCellInfoStruct.Depth = curCellStats.Centroid(3) / size(annotationData.imageData{1},3);
        [curCellInfoVec , cellInfoLabelList] = ConvertFeatureStructToFeatureVec(curCellInfoStruct);        
        
        cellInfoDataLocal = cat(1, cellInfoDataLocal, curCellInfoVec');
        
        % save images if requested
        if PARAMETERS.flagSaveImages 
            
            if PARAMETERS.flagSaveErrorImagesOnly 
                
                if strcmp(curAnnotatedFeatureClass, curPredictedFeatureClass) 
                    
                    if strcmp(curAnnotatedFeatureClass, strInlierClassName)
                        continue; % dont save images if the classifier prediction was correct
                    end
                    
                end
                
            end
            
            [afpathstr, afname, ~] = fileparts( annotationData.dataFilePath{1} );
            curImageOutputDir = fullfile(outputRootDir, 'images', [curAnnotatedFeatureClass, '_' curPredictedFeatureClass], sprintf('fid_%.2d_%s', fid, strtrim(afname)) );

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
            imwrite( genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 ), ...
                     fullfile(curImageOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   
            
            channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
                     fullfile(curImageOutputDir, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   
            
            imwrite(imCurCellMIP, fullfile(curImageOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   
                 
        end
        
   end
   
   fprintf( '\nclass distribution:\n' );
   tabulate( annotatedFeatureClassLocal ) 
   
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

fprintf( '\nTotal class distribution: \n' );
tabulate( annotatedFeatureClass );

% compute global performance stats
fprintf( '\nComputing overall performance stats: \n' );

globalPerformanceStats = ComputeClassificationPerformance_MultiClass( predictedFeatureClass, annotatedFeatureClass, classNameList )
globalConfusioinMatrix_Display = globalPerformanceStats.cMatrix_Display

performanceData = [];
performanceData.annotationFileList = annotationFileList;
performanceData.rootDirList = rootDirList;
performanceDatga.modelFilePath = modelFilePath;
performanceData.classNameList = classNameList;
performanceData.trueLabels = annotatedFeatureClass;
performanceData.predictedLabels = predictedFeatureClass;
performanceData.stats = globalPerformanceStats;
performanceData.cellInfoData = cellInfoData;
performanceData.subPerformanceData = performanceDataGlobal;

save( fullfile(outputRootDir, 'cellInlierDetectionModelPerformance.mat'), 'performanceData' );

% write global files
fprintf( '\n>>Writing csv file for inspection ... \n' );

featureFileNamePrefix = 'cellInlierDetectionPerformance';
featureLabels = cat(1, 'TrueClassLabel', 'PredictedClassLabel', cellInfoLabelList);
featureMatrix = cat(2, annotatedFeatureClass, predictedFeatureClass, cellInfoData);

WriteFeatureMatrixToCSVFile( fullfile(pwd, [featureFileNamePrefix '.csv']), featureMatrix, featureLabels );
movefile( fullfile(pwd, [featureFileNamePrefix '.csv']), fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );

% switch off diary
diary off;
