function [flagSuccess] = performCellCycleAnalysis(nuclearMarkerFilePath, fucciFilePath, regionMergingModelFile, cellCycleModelFile, resultsDir, varargin)

    p = inputParser;
    p.addParamValue( 'cellCycleModelClass', @CellCycleStateClassifier_OneStage, @(x) isa(x, 'function_handle') );
    p.addParamValue( 'metaInfoStruct', [], @(x) (isstruct(x) && all(isfield(x, {'header', 'data'}))) );
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagSaveImages', false, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'cellDiameterRange', [8, 20], @(x) (isnumeric(x) && numel(x) == 2) );
    p.addParamValue( 'minCellVolume', 400, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'minCellROIOverlap', 0.5, @(x) (isscalar(x) && x >= 0.0 && x <= 1.0) );
    p.addParamValue( 'minCellBBoxSizeImsp', [7 7 4], @(x) (isnumeric(x) && numel(x) == 3));
    p.addParamValue( 'flagIgnoreCellsOnXYBorder', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'finishStatusReportFile', [], @(x) (ischar(x)) );
    p.addParamValue( 'defaultVoxelSpacing', [], @(x) (isnumeric(x) && numel(x) == 3) );
    p.parse( varargin{:} );
    PARAMETERS = p.Results;      

    if ~isdir( resultsDir ) 
        mkdir(resultsDir);
    end

    % setup diary for logging analysis progress
    diary off;
    diary_file = fullfile( resultsDir , sprintf( '%s.log' , mfilename ) );
    if exist( diary_file, 'file' )
        delete( diary_file );
    end
    diary_file
    diary( diary_file );
    diary on;

    % check if all file provided exist
    if ~exist( nuclearMarkerFilePath, 'file' )
        error( 'unable to find nuclear marker file - %s', nuclearMarkerFilePath );
    end

    if ~exist( fucciFilePath, 'file' )
        error( 'unable to find fucci file - %s', fucciFilePath );
    end

    if ~exist( regionMergingModelFile, 'file' )
        error( 'unable to find region merging model file - %s', regionMergingModelFile );
    end

    if ~exist( cellCycleModelFile, 'file' )
        error( 'unable to find cell cycle classification model file - %s', cellCycleModelFile );
    end

    PARAMETERS
    nuclearMarkerFilePath
    fucciFilePath
    regionMergingModelFile
    cellCycleModelFile
    metaInfoStruct = PARAMETERS.metaInfoStruct
    
    if PARAMETERS.flagSaveImages
        if isdir( fullfile(resultsDir, 'images') )            
            rmdir( fullfile(resultsDir, 'images'), 's' );
        end
        mkdir( fullfile(resultsDir, 'images') );
    end

    AddWekaClassesToPath();

    matlabpath
    javaclasspath('-dynamic')
    version -java

    if ~isempty(PARAMETERS.finishStatusReportFile)
        fidStatus = fopen( PARAMETERS.finishStatusReportFile, 'w' );       
    end
    
    % load input data
    dataLoadTimer = tic;

        % load nuclear marker data    
        PrettyPrintStepDescription( 'Loading Nuclear Marker Data' );

            imageSeriesNucleus = loadIntravitalDataset( nuclearMarkerFilePath );

            metadata_NucleusChannel = imageSeriesNucleus(1).metadata
    
        % load fucci data
        PrettyPrintStepDescription( 'Loading Fucci Data' );

            imageSeriesFUCCI = loadIntravitalDataset( fucciFilePath  );    

            metadata_FUCCI = imageSeriesFUCCI(1).metadata
            
        % run sanity checks on metadata
        metadata_NucleusChannel.channelColors = [0, 0, 1];        
        metadata_FUCCI.channelColors = [ 0 1 0; 1 0 0 ];
        flagMetadataProblem = false;
         
            % make sure fucci data has two channels
            if metadata_FUCCI.numChannels ~= 2
                error( 'Fucci data is expected to contain 2 channels.' );
            end
        
            % make sure that the volume size of histone and fucci data is the same
            if any( metadata_FUCCI.volSize ~= metadata_NucleusChannel.volSize )
                error( 'Volume Size of Nuclear marker data doesnt match with Fucci data' );
            end
            
            % check voxel spacing in both histone and fucci data
            flagBadNuclearMarkerVoxelSpacing = false;
            if (metadata_NucleusChannel.voxelSpacing(1)/ metadata_NucleusChannel.voxelSpacing(3)) >= 1 || ...
               (metadata_NucleusChannel.voxelSpacing(2)/ metadata_NucleusChannel.voxelSpacing(3)) >= 1

                flagBadNuclearMarkerVoxelSpacing = true;

            end
        
            flagBadFucciVoxelSpacing = false;
            if (metadata_FUCCI.voxelSpacing(1)/ metadata_FUCCI.voxelSpacing(3)) >= 1 || ...
               (metadata_FUCCI.voxelSpacing(2)/ metadata_FUCCI.voxelSpacing(3)) >= 1     
           
                flagBadFucciVoxelSpacing = true;
                
            end 
            
            if flagBadNuclearMarkerVoxelSpacing && ~flagBadFucciVoxelSpacing
                
                warning( 'Voxel size in X- or Y- dimension of nuclear marker data is bigger than the Z-dimension. Using the spacing of fucci data.' );
                metadata_NucleusChannel.voxelSpacing = metadata_FUCCI.voxelSpacing;
                flagMetadataProblem = true;
                
            end
                
            if flagBadFucciVoxelSpacing && ~flagBadNuclearMarkerVoxelSpacing
                
                warning( 'Voxel size in X- or Y- dimension of fucci data is bigger than the Z-dimension. Using the spacing of nuclear marker data.' );
                metadata_FUCCI.voxelSpacing = metadata_NucleusChannel.voxelSpacing;       
                flagMetadataProblem = true;
                
            end
            
            if flagBadNuclearMarkerVoxelSpacing && flagBadFucciVoxelSpacing
            
                if isempty(PARAMETERS.defaultVoxelSpacing)
                    error( 'Voxel size in X- or Y- dimension of both nuclear marker and fucci data is bigger than the Z-dimension.' );
                else
                    warning( 'Voxel size in X- or Y- dimension of both nuclear marker and fucci data is bigger than the Z-dimension. Using the specified default value.' );
                    metadata_NucleusChannel.voxelSpacing = PARAMETERS.defaultVoxelSpacing;
                    metadata_FUCCI.voxelSpacing = PARAMETERS.defaultVoxelSpacing;
                    flagMetadataProblem = true;
                end
                
            end
            
            % channel excitation wavelength
            if isempty( metadata_NucleusChannel.channelExcitationWavelength )
                warning( 'metadata of nuclear marker data does not contain channel excitation wavelength. Using a default value of 830' );
                metadata_NucleusChannel.channelExcitationWavelength = 830;
                flagMetadataProblem = true;
            end
            
            if isempty( metadata_FUCCI.channelExcitationWavelength )
                
                warning( 'metadata of fucci data does not contain channel excitation wavelength. Using default values of 473 (GFP) and 559 (RFP), respectively' );
                metadata_FUCCI.channelExcitationWavelength = [473 559];
                flagMetadataProblem = true;
                
            else
                
                if metadata_FUCCI.channelExcitationWavelength(1) > metadata_FUCCI.channelExcitationWavelength(2)

                    % first channel is not green so swap it 
                    warning( 'GFP was not the first channel in confocal data. It will be swapped to the first place');
                    imageSeriesFUCCI(1).imageData = imageSeriesFUCCI(1).imageData(:,[2,1]);
                    metadata_FUCCI.channelExcitationWavelength = metadata_FUCCI.channelExcitationWavelength([2,1]);
                    flagMetadataProblem = true;
                    
                end
                
            end            
            
        % metadata
        metadata = metadata_NucleusChannel;    
        metadata.channelExcitationWavelength = [ metadata_NucleusChannel.channelExcitationWavelength, metadata_FUCCI.channelExcitationWavelength ];
        metadata.channelColors = cat( 1 , metadata_NucleusChannel.channelColors, metadata_FUCCI.channelColors );             

        % store meta info into inspection file
        [pname, nuclearMarkerFileName] = fileparts(nuclearMarkerFilePath);
        stackInfoStruct.nuclearMarkerFilePath = nuclearMarkerFilePath;
        stackInfoStruct.nuclearMarkerFileName = nuclearMarkerFileName;
            
        [pname, fucciFileName, ~] = fileparts(fucciFilePath);
        stackInfoStruct.fucciFilePath = fucciFilePath;
        stackInfoStruct.fucciFileName = fucciFileName;
        
        volVoxel = prod(metadata.voxelSpacing);        
        stackInfoStruct.volSize = metadata.volSize;
        stackInfoStruct.voxelSpacing = metadata.voxelSpacing;
        stackInfoStruct.volSizePhysp = metadata.volSize .* metadata.voxelSpacing;
        stackInfoStruct.stackVolumePhysp = prod(metadata.volSize) * volVoxel;    

        stackInfoStruct.flagMetadataProblem = flagMetadataProblem;
        
        compInfo.dataLoadTime = toc(dataLoadTimer);
        
        diaryFlush();
        
    % analyze data
    analysisTimer = tic;

        % shift correction
        PrettyPrintStepDescription( 'Correcting shift between nuclear marker data and fucci data' );
        
        shiftCorrectionTimer = tic; 
        [imageData, imRegValidMask, ...
         transformParameters] = CorrectTwoPhotonConfocalStageShift( imageSeriesNucleus(1).imageData, metadata_NucleusChannel.voxelSpacing, ...
                                                                    imageSeriesFUCCI(1).imageData, metadata_FUCCI.voxelSpacing, ...
                                                                    'flagParallelize', PARAMETERS.flagParallelize );
        
        stackInfoStruct.shift_correction.gfp.shiftTransform = transformParameters{1};
        stackInfoStruct.shift_correction.gfp.shiftAmount3D = norm(transformParameters{1});
        stackInfoStruct.shift_correction.gfp.shiftAmountXY = norm(transformParameters{1}(1:2));
        
        
        stackInfoStruct.shift_correction.rfp.shiftTransform = transformParameters{2};
        stackInfoStruct.shift_correction.rfp.shiftAmount3D = norm(transformParameters{2});
        stackInfoStruct.shift_correction.rfp.shiftAmountXY = norm(transformParameters{2}(1:2));
        
        stackInfoStruct.shift_correction.validROIVolPercent = 100 * mean(imRegValidMask(:));
        stackInfoStruct.shift_correction.validXPercent = 100 * mean( sum(sum(imRegValidMask, 1), 3) > 0 );
        stackInfoStruct.shift_correction.validYPercent = 100 * mean( sum(sum(imRegValidMask, 2), 3) > 0 );
        stackInfoStruct.shift_correction.validZPercent = 100 * mean( sum(sum(imRegValidMask, 1), 2) > 0 );
        
        compInfo.shiftCorrectionTime = toc(shiftCorrectionTimer);
        
        diaryFlush();
        
        if PARAMETERS.flagSaveImages

            imageDataBefore = [imageSeriesNucleus(1).imageData, imageSeriesFUCCI(1).imageData];
            channelColorMap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
            
            imwrite( generateMultichannelMIPImage(imageData, [], metadata.voxelSpacing), ...
                     fullfile(resultsDir, 'images', 'stackMIP.png'), 'png');
            
            imShiftDisplay = generateMultichannelOrthosliceImage(imageDataBefore, round(0.5 * metadata.volSize), channelColorMap, metadata.voxelSpacing);
            imShiftDisplay(:,end+1:end+20,:) = 1;
            imShiftDisplay = cat(2, imShiftDisplay, generateMultichannelOrthosliceImage(imageData, round(0.5 * metadata.volSize), channelColorMap, metadata.voxelSpacing) );
            
            imwrite( imShiftDisplay, fullfile(resultsDir, 'images', 'shiftCorrection.png'), 'png');
                 
            clear imageDataBefore;
            
        end
        
        clear imageSeriesNucleus;
        clear imageSeriesFUCCI;
        
        % segment nuclei
        PrettyPrintStepDescription( 'Running nuclei segmentation algorithm' );
        
        segTimer = tic;
        
        [imLabelCellSeg, imCellSeedPoints, ...
         segAlgoParameters ] = segmentCellsInIntravitalData( imageData{1}, ...
                                                             metadata.voxelSpacing, ...                                                                      
                                                             'flagParallelize', PARAMETERS.flagParallelize, ...
                                                             'flagDebugMode', false, ...
                                                             'cellDiameterRange', PARAMETERS.cellDiameterRange, ...
                                                             'thresholdingAlgorithm', 'MinErrorPoissonSliceBySliceLocal', ...
                                                             'seedPointDetectionAlgorithm', 'AdaptiveMultiscaleLoG', ...
                                                             'minCellVolume', PARAMETERS.minCellVolume, ...
                                                             'flagIgnoreCellsOnXYBorder', PARAMETERS.flagIgnoreCellsOnXYBorder, ...
                                                             'roiMask', imRegValidMask, ...
                                                             'minCellROIOverlap', PARAMETERS.minCellROIOverlap, ...
                                                             'regionMergingModelFile', regionMergingModelFile);

        numCells = max(imLabelCellSeg(:));       
        cellStats = regionprops( imLabelCellSeg, 'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
        
        fprintf( '\nThe segmentation algorithm found %d cells\n', numCells );
        
        compInfo.segmentationTime = toc(segTimer);        
        stackInfoStruct.CellCount = numCells;
        
        diaryFlush();
        
        % identify cell cycle state of each cell
        PrettyPrintStepDescription( 'Identifying cell cycle state of each cell' );
 
        cellCycleClassificationModel = PARAMETERS.cellCycleModelClass( cellCycleModelFile );
     
        classNameList = cellCycleClassificationModel.getClassNameList()
        compInfo.cellCycleClassificationTime = 0;
        
            % pre-processing
            fprintf('\n>> Pre-processing the image data ... \n');       

            imageDataAdjusted = CellCycleStateClassifier.preprocessImageData( imageData );

            % predict the state of each individual cell
            fprintf('\n>> Applying classification model on each of the %d cells ...\n', numCells);       
            
            predictedClassLabels = cell(numCells, 1);
            classPredictionProbabilities = cell(numCells, 1);  
            cellFeatureStruct = cell( numCells, 1);
            
            classificationTimer = tic;

            if PARAMETERS.flagParallelize                
                
                parfor cellId = 1:numCells

                    AddWekaClassesToPath();
                    curCellCycleModel = PARAMETERS.cellCycleModelClass( cellCycleModelFile );
                    [predictedClassLabels{cellId}, ...
                     classPredictionProbabilities{cellId}, ...
                     cellFeatureStruct{cellId}] = curCellCycleModel.predictCell(imageDataAdjusted, imRegValidMask, imLabelCellSeg, cellId, metadata.voxelSpacing); 

                end
                
            else
                
                fprintf( '\nProgress: \n' );

                last_percent_done = 0;
                numPrint = 0;

                for cellId = 1:numCells

                    [predictedClassLabels{cellId}, ...
                     classPredictionProbabilities{cellId}, ...
                     cellFeatureStruct{cellId}] = cellCycleClassificationModel.predictCell(imageDataAdjusted, imRegValidMask, imLabelCellSeg, cellId, metadata.voxelSpacing); 

                     percent_done = round(100*cellId/numCells);       

                     if percent_done > last_percent_done
                         fprintf( '%.2d%%  ', percent_done );
                         last_percent_done = percent_done;
                         numPrint = numPrint + 1;
                         if mod( numPrint, 10 ) == 0
                            fprintf( '\n' ); 
                         end
                     end        

                     diaryFlush();
                     
                end
                
            end
            
            compInfo.cellCycleClassificationTime = toc(classificationTimer);
            
            tabulate(predictedClassLabels)
            
            diaryFlush();
            
            % store overall class distribution in stack info
            for cid = 1:numel(classNameList)
               
                curClassCount = sum( ismember(predictedClassLabels, classNameList{cid}) );
                curClassPercentage = 100.0 * curClassCount / numel(predictedClassLabels);
                
                classCount.(classNameList{cid}) = curClassCount;
                classPercentage.(classNameList{cid}) = curClassPercentage;
                
            end
            
            stackInfoStruct.classDistribution.Count = classCount;
            stackInfoStruct.classDistribution.Percentage = classPercentage;            
            
            % compute and store some cell info
            fprintf('\n>> Computing some information about each cell ...\n');       
            
            cellInfoStruct = [];
            cellInfoFeatureMatrix = {};
            
            fprintf( '\nProgress: \n' );
            last_percent_done = 0;
            numPrint = 0;
            
            for cellId = 1:numCells

                curCellStruct.cellId = cellId;
                
                % compute region properties
                curCellProps = ComputeRegionProperties( imLabelCellSeg, cellId );

                % class
                curCellStruct.predictedCellCycleState = predictedClassLabels{cellId};

                curCellStruct.predictionProbability = max( classPredictionProbabilities{cellId} );
                
                for classId = 1:numel(classNameList)
                    curFieldName = ['cellCycleStatePredictionProbability_', classNameList{classId}];
                    curCellStruct.(curFieldName) = classPredictionProbabilities{cellId}(classId);
                end

                % note down whether or not the cell touches the border
                ptCellPixelCoord = ind2submat( size(imLabelCellSeg), curCellProps.PixelIdxList);

                curCellStruct.flagIsCellOnBorder = zeros(1,3);
                for dim = 1:3
                    if any( ismember(ptCellPixelCoord(:,dim), [1, size(imLabelCellSeg, dim)]) )
                        curCellStruct.flagIsCellOnBorder(dim) = 1;
                    end
                end

                curCellStruct.flagIsCellOnXYBorder = any(curCellStruct.flagIsCellOnBorder(1:2));
                curCellStruct.flagIsCellOnZBorder = any(curCellStruct.flagIsCellOnBorder(3));
                curCellStruct.volPercentInsideValidROI = 100.0 * mean( imRegValidMask(curCellProps.PixelIdxList) > 0 );

                % some basic properties
                basicProps.CentroidImsp = cellStats(cellId).Centroid;
                basicProps.DepthPhysp = (cellStats(cellId).Centroid(3)-1) .* metadata.voxelSpacing(3);

                basicProps.BBoxSizeImsp = cellStats(cellId).BoundingBox(4:end);
                basicProps.BBoxSizePhysp = cellStats(cellId).BoundingBox(4:end) .* metadata.voxelSpacing; 

                basicProps.Volume = curCellProps.Area * volVoxel;       
                basicProps.VolumeOfConvexHull = curCellProps.ConvexArea * volVoxel;
                basicProps.Convexity = curCellProps.Area / curCellProps.ConvexArea;

                ptCell = ind2submat( size(imLabelCellSeg), curCellProps.PixelIdxList ); 
                ptCell = bsxfun(@times, ptCell, metadata.voxelSpacing);
                ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
                [U, S, V] = svd( (ptCell' * ptCell) / size(ptCell,1) );
                curCellEllipsoidRadiusPhysp = zeros(1,3); 
                for j = 1:3
                    curCellEllipsoidRadiusPhysp(j) = 2 * sqrt(S(j,j)); % eigen-values are a measure of variance
                end
                basicProps.ellipsoidRadiusPhysp = curCellEllipsoidRadiusPhysp;        
                basicProps.fitterEllipsoidVolume = (4/3) * pi * prod(curCellEllipsoidRadiusPhysp);
                basicProps.ellipticVariance = ComputeEllipticVariance( ptCell, metadata.voxelSpacing );

                curCellStruct.properties = basicProps;      

                % intensity stats
                for chid = 1:3
                    cellPixelIntensities = imageDataAdjusted{chid}( curCellProps.PixelIdxList );
                    curIntensityStats.median(chid) = median( cellPixelIntensities );
                    curIntensityStats.mad(chid) = mad( cellPixelIntensities );
                    curIntensityStats.iqr(chid) = iqr( cellPixelIntensities );
                    curIntensityStats.robustSkewness(chid) = skewnessRobustHinkley( cellPixelIntensities );
                    curIntensityStats.robustKurtosis(chid) = kurtosisRobustCrow( cellPixelIntensities );
                end

                curCellStruct.intensityStats = curIntensityStats;

                % store features
                curCellStruct.features = cellFeatureStruct{cellId};

                % append to cellInfoStruct
                cellInfoStruct = [cellInfoStruct ; curCellStruct];                    
                [featureVec, cellInfoFeatureNameList] = ConvertFeatureStructToFeatureVec(curCellStruct);
                cellInfoFeatureMatrix = cat(1, cellInfoFeatureMatrix, featureVec);
                    
                percent_done = round(100*cellId/numCells);       

                if percent_done > last_percent_done
                    fprintf( '%.2d%%  ', percent_done );
                    last_percent_done = percent_done;
                    numPrint = numPrint + 1;
                    if mod( numPrint, 10 ) == 0
                       fprintf( '\n' ); 
                    end
                end        
               
                diaryFlush();
                
            end
        
            % compute and store some cell info
            if PARAMETERS.flagSaveImages 

                fprintf('\n>> Saving snapshot images of cells ...\n');       

                imageSaveTimer = tic;

                if PARAMETERS.flagParallelize
                
                    parfor cellId = 1:numCells
                        
                        imageOutputDir = fullfile(resultsDir, 'images', predictedClassLabels{cellId});

                        if ~isdir( imageOutputDir )
                            mkdir( imageOutputDir );
                        end
                        
                        WriteCellSnapshotImages(imageData, (imLabelCellSeg == cellId), cellStats(cellId), cellId, imageOutputDir);
                        
                    end

                else
                    
                    fprintf( '\nProgress: \n' );

                    last_percent_done = 0;
                    numPrint = 0;
                    
                    for cellId = 1:numCells
                        
                        imageOutputDir = fullfile(resultsDir, 'images', predictedClassLabels{cellId});

                        if ~isdir( imageOutputDir )
                            mkdir( imageOutputDir );
                        end
                        
                        WriteCellSnapshotImages(imageData, (imLabelCellSeg == cellId), cellStats(cellId), cellId, imageOutputDir);
                        
                        percent_done = round(100*cellId/numCells);       

                        if percent_done > last_percent_done
                            fprintf( '%.2d%%  ', percent_done );
                            last_percent_done = percent_done;
                            numPrint = numPrint + 1;
                            if mod( numPrint, 10 ) == 0
                               fprintf( '\n' ); 
                            end
                        end     
                        
                        diaryFlush();
                        
                    end
                    
                end
                
                WriteCellClassificationSummaryImage(imageData, imLabelCellSeg, classNameList, ...
                                                    predictedClassLabels, classPredictionProbabilities, fullfile(resultsDir, 'images'));                
                                                
                compInfo.imageSaveTime = toc(imageSaveTimer);
                
            end        

       % compute and store image quality metrics
       PrettyPrintStepDescription( 'Computing image quality metrics' );

            % cell density
           fprintf( '\n>> Running a thresholding algorithm to get rough estimates of cell density ...\n' );   
           imThresh = segmentCellForegroundUsingLocalMinError( imageDataAdjusted{1}, 30, ...
                                                               'model', 'poisson', ...  
                                                               'flagParallelize', PARAMETERS.flagParallelize );

           stackQuality.foregroundVolumePhysp = sum( imThresh(:) > 0 ) * prod(metadata.voxelSpacing);
           stackQuality.foregroundVolPercent = 100.0 * mean( imThresh(:) > 0 );

           fprintf( '\n\tStack contains %.2f%% foreground ...\n', stackQuality.foregroundVolPercent );   

           diaryFlush();
           
           % noise level
           fprintf( '\n>> Estimating noise level ...\n\n' );   

           stackQuality.noiseAgainstMedianFilteredResponse = zeros(1,3);

           for chid = 1:3
                imSmooth = matitk( 'FMEDIAN', [1,1,1], imageData{chid} );   
                stackQuality.noiseAgainstMedianFilteredResponse(chid) = sqrt(mean((imSmooth(:) - imageData{chid}(:)).^2)); 
           end

           diaryFlush();
           
           % amount of blurring
           stackQuality.radialPowerSpectrumLogLogSlope = zeros(1,3);
           
           for chid = 1:3
                stackQuality.radialPowerSpectrumLogLogSlope(chid) = ComputePowerSpectrumSlope( imageData{chid}, metadata.voxelSpacing );
           end

           for chid = 1:3
               
                curSliceLogLogSlope = zeros(size(imageData{chid},3),1);
                
                for sliceId = 1:size(imageData{chid},3)
                    curSliceLogLogSlope(sliceId) = ComputePowerSpectrumSlope( imageData{chid}(:,:,sliceId), metadata.voxelSpacing(1:2) ); 
                end
                
                curSlopeStats.mean = mean(curSliceLogLogSlope);
                curSlopeStats.stddev = std(curSliceLogLogSlope);
                curSlopeStats.min = min(curSliceLogLogSlope);
                curSlopeStats.max = max(curSliceLogLogSlope);
                
                stackQuality.radialPowerSpectrumLogLogSlopePerSlice(chid) = curSlopeStats;
                
           end
           
       stackInfoStruct.quality = stackQuality;

    compInfo.totalAlgorithmTime = compInfo.shiftCorrectionTime + compInfo.segmentationTime + compInfo.cellCycleClassificationTime;
    compInfo.totalAnalysisTime = toc(analysisTimer);

    stackInfoStruct.computation = compInfo;
    stackInfoStruct.resultsDir = resultsDir;
    
    stackInfoStruct
    
    stackQuality

    compInfo

    metaInfoStruct
    
    % save stuff to files
    
        % write stack info
        fprintf( '\n>>Writing stack analysis csv file ... \n' );

        [ stackFeatureVec , stackFeatureNameList ] = ConvertFeatureStructToFeatureVec( stackInfoStruct );
        
        if ~isempty( PARAMETERS.metaInfoStruct )
            stackFeatureNameList = cat(2, metaInfoStruct.header, stackFeatureNameList);            
            stackFeatureVec = cat(2, metaInfoStruct.data , stackFeatureVec);            
        end
        
        WriteFeatureMatrixToCSVFile( fullfile(resultsDir, 'stackAnalysisInfo.csv' ), ...
                                     stackFeatureVec, stackFeatureNameList );

        diaryFlush();
                                 
        % write cell info
        fprintf( '\n>>Writing cell analysis csv file ... \n' );

        featureNameList = cat(2, stackFeatureNameList, cellInfoFeatureNameList );
        featureMatrix = cat(2, repmat(stackFeatureVec, numCells, 1), cellInfoFeatureMatrix );
        
        WriteFeatureMatrixToCSVFile( fullfile(resultsDir, 'cellAnalysisInfo.csv' ), ...
                                     featureMatrix, featureNameList );

        % save analysis info in a mat file for inspection later
        save( fullfile(resultsDir, 'cellCycleAnalysisInfo.mat'), 'stackInfoStruct', 'cellInfoStruct', 'PARAMETERS', 'segAlgoParameters', 'regionMergingModelFile', 'cellCycleModelFile' );

        % save mat file that can be loaded into the CellStateAnalzer tool
        analysisData.dataFilePath = {nuclearMarkerFilePath, fucciFilePath};
        analysisData.metadata = metadata;
        analysisData.imageData = imageData;
        analysisData.imRegValidMask = imRegValidMask;
        
        analysisData.imLabelCellSeg = imLabelCellSeg;
        analysisData.imCellSeedPoints = imCellSeedPoints;
        [imSegMaskRGB, analysisData.CellSegColorMap] = label2rgbND(imLabelCellSeg);
        
        analysisData.cellPatternTypes = classNameList;
        analysisData.cellCycleModelFile = cellCycleModelFile;
        analysisData.regionMergingModelFile = regionMergingModelFile;
        
        cellStats = regionprops( imLabelCellSeg, 'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
        for cellId = 1:numel(cellStats)
            cellStats(cellId).cellId = cellId;
            cellStats(cellId).cellPatternType = cellInfoStruct(cellId).predictedCellCycleState; 
            cellStats(cellId).classPredictionProbabilities = classPredictionProbabilities{cellId};
        end        
        analysisData.cellStats = cellStats;
        
        save( fullfile(resultsDir, 'CellStateAnalyzer.mat'), '-struct', 'analysisData' );

    % report success
    flagSuccess = true;
    fprintf( '\n\n---- Analysis Succeeded ---- \n\n' );
    
    if ~isempty(PARAMETERS.finishStatusReportFile)
        fprintf( fidStatus, 'Analysis Succeeded' );
        fclose(fidStatus);
    end
    
    diary off;
        
end


function WriteCellClassificationSummaryImage(imageData, imLabelCellSeg, classNameList, ...
                                             predictedClass, classPredictionProbabilities, outDir)

    gapSize = 20;
    szOutputImage = [100, 100];    
    textBoxHeight = 30;
    channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
    
    if size(classPredictionProbabilities{1}, 1) == numel(classNameList)
        classPredictionProbabilities = (reshape(cell2mat(classPredictionProbabilities), [numel(classNameList), numel(predictedClass)]))';
    else
        classPredictionProbabilities = cell2mat(classPredictionProbabilities);
    end
    
    cellStats = regionprops(imLabelCellSeg, {'Centroid', 'BoundingBox'} );
    imOutput = [];
    
    for lid = 1:numel(classNameList)
    
        indCellsInCurClass = find( strcmpi(predictedClass, classNameList{lid}) );
        
        if isempty(indCellsInCurClass)
            continue;
        end
        
        [p, ix] = sortrows( -1 * classPredictionProbabilities(indCellsInCurClass,:), lid );
        indCellsInCurClass = indCellsInCurClass(ix);
        
        if ~isempty(imOutput)
            imOutput = padarray(imOutput, [2*gapSize, 0, 0], 1.0, 'post');
            hzbarpos = size(imOutput,1) - gapSize;
            imOutput(hzbarpos:hzbarpos+2, :, :) = 0;
        end
        
        imCurStrip = repmat(ones( [szOutputImage(1) + textBoxHeight, szOutputImage(2)*48 + 15 * gapSize]), [1, 1, 3]);
        numStripImages = 0;
        nextImageCol = 1;
        
        for clsCellId = 1:numel(indCellsInCurClass)
            
            cellId = indCellsInCurClass(clsCellId);
            curCellStats = cellStats(cellId);
            curCellCentroid = curCellStats.Centroid;
            curCellBoundingBox = curCellStats.BoundingBox;
            curCellDisplaySize = max( [curCellBoundingBox(4:5), 70] );
            imCellMask = (imLabelCellSeg == cellId);
            
            % crop images
            subinds = cell(1,3);
            imsize = size(imageData{1});
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
                channelDisplayRange = ComputeImageDynamicRange( imageData{j}, 98.0 );
                imCurChannelCropped = mat2gray( imageData{j}( subinds{:} ), channelDisplayRange );
                imCurChannelCropped = imresize( imCurChannelCropped, szOutputImage );
                imCurCellMidSliceAllChannel = cat(3, imCurCellMidSliceAllChannel, imCurChannelCropped );
            end

            imCurCellSegBndMidSliceCropped = imresize( bwperim( imCellMask( subinds{:} ) ), szOutputImage, 'nearest' );

            imCurCellMIP = mat2gray( max( imageData{1}(subinds{1:2}, :) .* imCellMask(subinds{1:2}, :), [], 3) );
            imCurCellMIP = repmat(imresize(imCurCellMIP, szOutputImage), [1, 1, 3]);

            imCurCellMidSlice_FUCCI = genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) );
            imCurCellMidSlice_Histone = genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 );
            %imCurCellMidSlice_AllChannel = genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap );
            
            imCurCellMidSliceDisplay = cat(2, imCurCellMIP, imCurCellMidSlice_FUCCI, imCurCellMidSlice_Histone);
            imCurCellMidSliceDisplay = padarray(imCurCellMidSliceDisplay, [textBoxHeight 0 0], 1.0, 'pre');

            % add text: <cell id> <probability> 
            textInserter = vision.TextInserter( sprintf('%d %.2f ', cellId, classPredictionProbabilities(cellId, lid)), 'FontSize', textBoxHeight-4);
            imCurCellMidSliceDisplay = step(textInserter, imCurCellMidSliceDisplay);

            % shown probability distribution
            imProb = [];
            probBoxSize = min( floor(textBoxHeight * 0.75), floor(szOutputImage(2) * 0.2) );
            for j = 1:numel(classNameList)
               imCurProb = repmat( (1-classPredictionProbabilities(cellId,j)) * ones( probBoxSize * [1,1] ), [1, 1, 3]);
               imCurProb = padarray(imCurProb, [1, 1, 0], 0);
               imProb = cat(2, imProb, imCurProb);
            end
            imCurCellMidSliceDisplay(1:size(imProb,1), (end-size(imProb,2)+1):end, :) = imProb;
            
            % add gap to its left
            if numStripImages > 0
                imCurCellMidSliceDisplay = padarray(imCurCellMidSliceDisplay, [0, gapSize, 0], 1.0, 'pre');
            end
            
            % insert into strip    
            imCurStrip(:, nextImageCol:(nextImageCol+size(imCurCellMidSliceDisplay,2)-1), :) = imCurCellMidSliceDisplay;
            numStripImages = numStripImages + 1;
            
            if numStripImages < 16 && clsCellId < numel(indCellsInCurClass)
                nextImageCol = nextImageCol + size(imCurCellMidSliceDisplay, 2);
            else
                imOutput = cat(1, imOutput, imCurStrip);
                imCurStrip = repmat(ones( [szOutputImage(1) + textBoxHeight, szOutputImage(2)*48 + 15 * gapSize]), [1, 1, 3]);
                numStripImages = 0;
                nextImageCol = 1;
            end
            
        end
        
    end

    imwrite( imOutput, fullfile(outDir, 'cellCycleClassificationSummary.png'), 'png' );
    
end

function WriteCellSnapshotImages(imageData, imCellMask, cellStats, cellId, imageOutputDir)

    curCellCentroid = cellStats.Centroid;
    curCellBoundingBox = cellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), 70] );
    szOutputImage = [100, 100];

    % crop images
    subinds = cell(1,3);
    imsize = size(imageData{1});
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
        channelDisplayRange = ComputeImageDynamicRange( imageData{j}, 98.0 );
        imCurChannelCropped = mat2gray( imageData{j}( subinds{:} ), channelDisplayRange );
        imCurChannelCropped = imresize( imCurChannelCropped, szOutputImage );
        imCurCellMidSliceAllChannel = cat(3, imCurCellMidSliceAllChannel, imCurChannelCropped );
    end

    imCurCellSegBndMidSliceCropped = imresize( bwperim( imCellMask( subinds{:} ) ), szOutputImage, 'nearest' );

    imCurCellMIP = mat2gray( max( imageData{1}(subinds{1:2}, :) .* imCellMask(subinds{1:2}, :), [], 3) );
    imCurCellMIP = imresize(imCurCellMIP, szOutputImage);

    % write images
    imwrite( genImageMaskOverlay( imCurCellMidSliceAllChannel(:,:,1), imCurCellSegBndMidSliceCropped, [1, 0, 0], 0.5 ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   

    channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
    imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   

    imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
             fullfile(imageOutputDir, sprintf('CellMidSliceAllChannel_%.3d.png', cellId)), 'png' );   

    imwrite(imCurCellMIP, fullfile(imageOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   

end
