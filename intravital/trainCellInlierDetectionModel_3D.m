
clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data';    
    defaultOutputDir = 'Z:\intravital\data';    
    
    PARAMETERS.flagIgnoreBadlySegmentatedCells = false;
    PARAMETERS.flagIgnoreCellsOnBorder = false;
    PARAMETERS.flagIgnoreApoptoticCells = true;
    PARAMETERS.flagConsiderApoptoticCellsAsOutliers = true;
    
    PARAMETERS.flagSaveOutlierImages = true;
    PARAMETERS.flagSaveInlierImages = false;
    
    strOutlierClassName = 'Outlier';
    strInlierClassName = 'Inlier';
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ask the user to provide the list of directories to be processed
rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);

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
outputRootDir

if PARAMETERS.flagSaveOutlierImages || PARAMETERS.flagSaveInlierImages
    if isdir( fullfile(outputRootDir, 'images') )            
        rmdir( fullfile(outputRootDir, 'images'), 's' );
    end
end

% process each annotation file
featureData = {};
featureClass = {};
cellInfoData = {};
featureNameList = {};

for fid = 1:numel(annotationFileList)
    
   fprintf( '\n\n*****************************************************\n\n' );           
   
   fprintf( '\n\nComputing Features for file %d/%d ...\n\n', fid, numel(annotationFileList) );       
    
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
   
   featureDataLocal = {};
   cellInfoDataLocal = {};
   featureClassLocal = {};
   
   for cellId = 1:numCells
       
        curCellStats = annotationData.cellStats(cellId);
        
        fprintf( '\n>> Computing Features for cell %d/%d of type - %s ...\n', cellId, numel(annotationData.cellStats), curCellStats.cellPatternType );       
        
        % ignore un-annotated cells
        if strcmpi( curCellStats.cellPatternType, 'None' )
            fprintf( '\nCell was not assigned to a pattern. Will be ignored from training.\n' );
            continue;            
        end
        
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
                
                curFeatureClass = strOutlierClassName;
                
            case { 'Apoptotic_Green', 'Apoptotic_Red', 'Apoptotic_Black' }
            
                if PARAMETERS.flagIgnoreApoptoticCells
                    continue;
                else
                    if PARAMETERS.flagConsiderApoptoticCellsAsOutliers
                        curFeatureClass = strOutlierClassName;
                    else
                        curFeatureClass = strInlierClassName;
                    end
                end
               
            otherwise
                
                curFeatureClass = strInlierClassName;
            
        end
        
        featureClassLocal = cat( 1, featureClassLocal, {curFeatureClass} );
        
        % compute features
        imCurCellMask = (imLabelCellSeg == cellId);
        [featureDataStruct, featureComputationParameters] = ComputeCellInlierDetectionFeatures(imageDataAdjusted, imLabelCellSeg, cellId, metadata.voxelSpacing );          
        [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureDataStruct );        
        featureDataLocal = cat(1, featureDataLocal, featureVec');

        % note down some cell identity information
        curCellInfoStruct.annotationFilePath = curAnnotationFilePath;
        curCellInfoStruct.cellId = cellId;
        curCellInfoStruct.cellPatternType = curCellStats.cellPatternType;
        curCellInfoStruct.Centroid = curCellStats.Centroid;
        curCellInfoStruct.Depth = curCellStats.Centroid(3) / size(annotationData.imageData{1},3);
        [curCellInfoVec , cellInfoLabelList] = ConvertFeatureStructToFeatureVec(curCellInfoStruct);        
        
        cellInfoDataLocal = cat(1, cellInfoDataLocal, curCellInfoVec');
        
        % save images if requested
        if (strcmp( curFeatureClass, strOutlierClassName ) && PARAMETERS.flagSaveOutlierImages) || ...
           (strcmp( curFeatureClass, strInlierClassName ) && PARAMETERS.flagSaveInlierImages) 

            [afpathstr, afname, ~] = fileparts( annotationData.dataFilePath{1} );
            curImageOutputDir = fullfile(outputRootDir, 'images', curFeatureClass, sprintf('fid_%.2d_%s', fid, strtrim(afname)) );

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
            
            % write images
            imwrite( genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 ), ...
                     fullfile(curImageOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   
            
            channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
                     fullfile(curImageOutputDir, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   
            
        end
        
   end
   
   fprintf( '\nclass distribution:\n' );
   tabulate( featureClassLocal ) 
   
   featureData = cat(1, featureData, featureDataLocal);       
   featureClass = cat(1, featureClass, featureClassLocal);
   cellInfoData = cat(1, cellInfoData, cellInfoDataLocal);
   
   computationTime = toc(feaCompTimer)   
    
end

fprintf( '\n\n*****************************************************\n\n' );

fprintf( '\nTotal class distribution: \n' );
tabulate( featureClass );

featureComputationParameters

% write global files
fprintf( '\n>>Writing feature files with both inliers and outliers ... \n' );

featureFileNamePrefix = 'cellInlierDetectionFeatures';

    % write csv file -- without debug info
    fprintf( '\nWriting csv files without cell info ... \n' );

    WriteFeatureMatrixToCSVFile( fullfile(pwd, [featureFileNamePrefix '.csv']), featureData, featureNameList, featureClass );
    movefile( fullfile(pwd, [featureFileNamePrefix '.csv']), fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );

    % write global csv file with cell info for debugging
    fprintf( '\nWriting global csv file with cell info ... \n' );

    featureLabels = cat(1, cellInfoLabelList, featureNameList);
    featureMatrix = cat(2, cellInfoData, featureData);

    WriteFeatureMatrixToCSVFile( fullfile(pwd, [featureFileNamePrefix '_info.csv']), featureMatrix, featureLabels, featureClass );
    movefile( fullfile(pwd, [featureFileNamePrefix '_info.csv']), fullfile(outputRootDir, [featureFileNamePrefix '_info.csv']) );
    
    % write arff file for weka
    fprintf( '\nConverting csv files to arff files ... \n' );

    AddWekaClassesToPath();

    ConvertCSVFileToArffFile( fullfile(outputRootDir, [featureFileNamePrefix '.csv']) );
    
% write files separately for inliers and outliers
classNameList = { strInlierClassName, strOutlierClassName }; 

for tid = 1:numel(classNameList)

    fprintf( '\n>>Writing feature files for class - %s ... \n', classNameList{tid} );

    curFeatureFileNamePrefix = sprintf( '%s_%s_only', featureFileNamePrefix, classNameList{tid} );

    flagCurClassFilter = strcmp( featureClass, classNameList{tid} );
    
    curFeatureData = featureData(flagCurClassFilter, :);
    curFeatureClass = featureClass(flagCurClassFilter, :);
    curCellInfoData = cellInfoData(flagCurClassFilter, :);
    
    % write csv file -- without debug info
    fprintf( '\nWriting csv files without cell info ... \n' );

    WriteFeatureMatrixToCSVFile( fullfile(pwd, [curFeatureFileNamePrefix '.csv']), curFeatureData, featureNameList, curFeatureClass);
    movefile( fullfile(pwd, [curFeatureFileNamePrefix '.csv']), fullfile(outputRootDir, [curFeatureFileNamePrefix '.csv']) );

    % write csv file -- with debugging info
    fprintf( '\nWriting csv files with cell info ... \n' );

    featureLabels = cat(1, cellInfoLabelList, featureNameList);
    featureMatrix = cat(2, curCellInfoData, curFeatureData);

    WriteFeatureMatrixToCSVFile( fullfile(pwd, [curFeatureFileNamePrefix '_info.csv']), featureMatrix, featureLabels, curFeatureClass);
    movefile( fullfile(pwd, [curFeatureFileNamePrefix '_info.csv']), fullfile(outputRootDir, [curFeatureFileNamePrefix '_info.csv']) );

    % write arff file for weka
    fprintf( '\nConverting csv files to arff files ... \n' );

    ConvertCSVFileToArffFile( fullfile(outputRootDir, [curFeatureFileNamePrefix '.csv']) );

end

save( fullfile(outputRootDir, [featureFileNamePrefix, '.mat']), ...
      'PARAMETERS', 'rootDirList', 'annotationFileList', ...
      'featureData', 'featureNameList', 'featureClass', 'featureComputationParameters' );

% switch off diary
diary off;
