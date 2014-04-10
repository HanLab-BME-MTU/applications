%function trainRegionMergingModel_3D

%     clc
%     clear 

    %***********************************************************************
    % 
    %                           PARAMETERS 
    %
    %***********************************************************************

        defaultDataRootDir = 'Z:\intravital\data';    
        defaultOutputDir = 'Z:\intravital\data';    

        PARAMETERS.flagParallelize = false;
        
        % less-critical parameters (leave defaults)
        PARAMETERS.threshNeighOverlap = 5;
        PARAMETERS.neighSearchRad = 2;    

    %***********************************************************************

    % ask the user to provide the list of directories to be processed
    rootDirList = uipickfiles('Prompt', 'Select folders containing the annotation files', ...
                              'REFilter', '(Cell).*\.mat$');

    if ~iscell(rootDirList)
        return;
    end
    
    % ask the user to select the output directory
    outputRootDir = uigetdir(defaultOutputDir, 'Select Output Directory' ); 

    if ~ischar(outputRootDir)
       return; 
    end
    
    hStatusDialog = waitbar(0, 'Building Model ...This might take a while, you may want to get some coffee.');
    
    % for parallelization
    if PARAMETERS.flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end    

    % make a list of files that need to be processed
    annotationFileList = [];
    
    for rid = 1:numel(rootDirList)

        curAnnotationFileList = rdir( fullfile(rootDirList{rid}, '**', 'CellSegmentationQualityAnnotation.mat') );
        annotationFileList = cat(1, annotationFileList, curAnnotationFileList );
        
        curAnnotationFileList = rdir( fullfile(rootDirList{rid}, '**', 'CellPatternAnnotation.mat') );
        annotationFileList = cat(1, annotationFileList, curAnnotationFileList );
        
    end
    
    if isempty(annotationFileList)        
        h = msgbox('No annotation files were found in the specified folders', ...
                   'Region merging model training', 'help', 'modal' );
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
    outputRootDir
    
    % process each annotation file
    numMerge = 0;
    numDontMerge = 0;
    
    featureData = [];
    featureClass = [];
    regionInfoData = [];    
        
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

       % get data
       imInput = annotationData.imageData{1};
       imLabelCellSeg = annotationData.imLabelCellSeg;
       imCellSeedPoints = annotationData.imCellSeedPoints;
       numCells = numel(annotationData.cellStats)       
       spacing = annotationData.metadata.voxelSpacing;
       
       % pre-processing

           % standardize the image data before computing any features
           imAdjusted = mat2gray(imInput);

           % denoise it a bit using a median filter
           imAdjusted = matitk('FMEDIAN', [1, 1, 1], imAdjusted);
           
       % compute basic cell properties 
       regPropsArr = cell(numCells, 1);

       if PARAMETERS.flagParallelize
           
           parfor cellId = 1:numCells

                regProps = ComputeRegionProperties(imLabelCellSeg, cellId, spacing);

                % add seed point location
                indCurRegionSeedPoint = regProps.PixelIdxList( imCellSeedPoints(regProps.PixelIdxList) > 0 );
                if isempty( indCurRegionSeedPoint )            
                    regProps.ptCellSeedPoint = regProps.ptCentroid;
                else
                    [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                    indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                    regProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
                end

                regPropsArr{cellId} = regProps;

           end
           
       else
           
           for cellId = 1:numCells

                regProps = ComputeRegionProperties(imLabelCellSeg, cellId, spacing);

                % add seed point location
                indCurRegionSeedPoint = regProps.PixelIdxList( imCellSeedPoints(regProps.PixelIdxList) > 0 );
                if isempty( indCurRegionSeedPoint )            
                    regProps.ptCellSeedPoint = regProps.ptCentroid;
                else
                    [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                    indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                    regProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
                end

                regPropsArr{cellId} = regProps;

           end
           
       end

       cellProps = [];
       for cellId = 1:numCells 
           cellProps = [cellProps; regPropsArr{cellId}];           
       end
       clear regPropsArr;
       
       % build region adjacency graph
       adjMatrix = sparse(numCells, numCells);
       adjMatrixOverseg = sparse(numCells, numCells);

       overseg_pid = find( strcmpi( annotationData.cellPatternTypes, 'Over_Segmentation' ) );
       neighStrel = streldisknd(PARAMETERS.neighSearchRad*ones(1, ndims(imInput)));
       
       for i = 1:numCells

           % get cropped label mask        
           imCurLabelCellSegCropped = imLabelCellSeg( cellProps(i).indCropBBox{:} );
           imCurCellMask = (imCurLabelCellSegCropped == i);

           % check if the cell has significant boundary overlap with any other cells
           imDilatedCellMask = imdilate( imCurCellMask, neighStrel );
           neighPixLabels = imCurLabelCellSegCropped( imDilatedCellMask - imCurCellMask > 0 );
           neighCellIdList = setdiff( unique(neighPixLabels), [0, i] );            
           neighCellIdList( neighCellIdList < i ) = [];

           % add edges to adjacency matrix
           adjMatrix(i, neighCellIdList) = 1;
           adjMatrix(neighCellIdList, i) = 1;

           % get annotated cell pattern type
           curCellPatternType = annotationData.cellStats(i).cellPatternType;
           curCellPatternId = annotationData.cellStats(i).cellPatternId;

           assert( strcmp(curCellPatternType, annotationData.cellPatternTypes(curCellPatternId)) );
           
           % look for overseg-overseg pairs
           neighCellPatternId = [annotationData.cellStats(neighCellIdList).cellPatternId];                
           neighCellPatternTypes = annotationData.cellPatternTypes( neighCellPatternId );
           
           neighOverSegId = neighCellIdList( neighCellPatternId == overseg_pid );
           
           adjMatrixOverseg(i, neighOverSegId) = 1;
           adjMatrixOverseg(neighOverSegId, i) = 1;
           
       end
       
       % compute region merging features
       [v1, v2, ~] = find( adjMatrix );   
       
       flagCellUsedForTraining = zeros(numel(annotationData.cellStats),1);

       ignore_pid = find( ismember( annotationData.cellPatternTypes, { 'Bad_Detection', 'Under_Segmentation', 'None' } ) );
       goodseg_pid = setdiff( 1:numel(annotationData.cellPatternTypes), [ignore_pid ; overseg_pid] );

       featureDataLocal = [];
       regionInfoDataLocal = [];
       featureClassLocal = [];
       
       for k = (find(v1 < v2))'
            
            c1 = v1(k);
            c2 = v2(k);

            patternId = [annotationData.cellStats(c1).cellPatternId, ...
                         annotationData.cellStats(c2).cellPatternId];
            
            if any( ismember( patternId, ignore_pid ) )
                continue;
            end

            % assign to merge/dont-merge class
            if all( patternId == overseg_pid )
                
               oversegNeighId = setdiff( find(adjMatrixOverseg(c1, :) | adjMatrixOverseg(c2, :)), [c1,c2] );
                
               if ~isempty(oversegNeighId)
                   continue;
               end
               
               curFeatureClass = 'Merge';
               
            else
                
               curFeatureClass = 'Dont-Merge';
               
            end
            
            featureClassLocal = cat( 1, featureClassLocal, {curFeatureClass} ); 
                
            % note down some info identifying regions for debugging
            curRegionInfoStruct.annotationFilePath = curAnnotationFilePath;
            curRegionInfoStruct.cellId = [c1, c2];
            curRegionInfoStruct.cellPatternType_1 = annotationData.cellStats(c1).cellPatternType;
            curRegionInfoStruct.cellPatternType_2 = annotationData.cellStats(c2).cellPatternType;
            [curRegionInfoVec , regInfoLabelList] = ConvertFeatureStructToFeatureVec(curRegionInfoStruct);        
            
            regionInfoDataLocal = cat(1, regionInfoDataLocal, curRegionInfoVec);
            
            % compute features
            [featureDataStruct, featureComputationParameters] = ComputeRegionMergingFeatures(imAdjusted, imLabelCellSeg, cellProps, c1, c2, 'spacing', spacing);          
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureDataStruct );        
            featureDataLocal = cat(1, featureDataLocal, featureVec);
            
       end     

       fprintf( '\nclass distribution:\n' );
       tabulate( featureClassLocal ) 
       
       featureData = cat(1, featureData, featureDataLocal);       
       featureClass = cat(1, featureClass, featureClassLocal);
       regionInfoData = cat(1, regionInfoData, regionInfoDataLocal);
       
       computationTime = toc(feaCompTimer)
       
       waitbar(fid / numel(annotationFileList), hStatusDialog);
       
    end

   fprintf( '\n\n*****************************************************\n\n' );
   
    if isempty(featureClass)        
        msgbox('No valid annotation files were found in the specified folders', ...
               'Region merging model training', 'error', 'modal' );
        closeStatusDialog(hStatusDialog);
        return;
    end
   
   fprintf( '\nTotal class distribution: \n' );
   tabulate( featureClass );
        
   featureComputationParameters
   
   % write csv file -- without debug info
    fprintf( '\nWriting csv file without region info ... \n' );
   
   featureLabels = cat(2, 'ClassLabel', featureNameList);
   featureMatrix = cat(2, featureClass, featureData);
   WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, 'regionMergingFeatures.csv'), featureMatrix, featureLabels );

   % write csv file -- with debugging info
    fprintf( '\nWriting csv file with region info ... \n' );
   
   featureLabels = cat(2, 'ClassLabel', regInfoLabelList, featureNameList);
   featureMatrix = cat(2, featureClass, regionInfoData, featureData);
   WriteFeatureMatrixToCSVFile( fullfile(outputRootDir, 'regionMergingFeaturesInfo.csv'), featureMatrix, featureLabels );
    
   % write arff file for weka
   fprintf( '\nConverting csv file to arff file ... \n' );
   
   AddWekaClassesToPath();

   ConvertCSVFileToArffFile( fullfile(outputRootDir, 'regionMergingFeatures.csv') );
   
   % generate classification model
   wekaTrainingDataset = getWekaDatasetFromCSVFile( fullfile(outputRootDir, 'regionMergingFeatures.csv') );
   wekaTrainingDataset.setClassIndex(0);
   wekaModel = BuildWekaModelForRegionMerging(wekaTrainingDataset, fullfile(outputRootDir, 'regionMerging.model'));
   
   %wrapup
   closeStatusDialog(hStatusDialog);
   diary off;
   
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
   
%end
    
    
    
    
    
    
    
    
    