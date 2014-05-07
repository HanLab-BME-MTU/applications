function [flagSuccess] = performDNADamageAnalysis(imageDataFilePath, channelId53BP1, ...
                                                  regionMergingModelFile, resultsDir, varargin)

    p = inputParser;

    p.addParamValue('cellDiameterRange', [8, 20], @(x) (isnumeric(x) && numel(x) == 2) );
    p.addParamValue('minCellVolume', 400, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'flagIgnoreCellsOnXYBorder', true, @(x) (isscalar(x) && islogical(x)) );
    
    p.addParamValue('fociDiameterRange',  0.621 * [3, 7], @(x) (isnumeric(x) && numel(x) == 2) );
    p.addParamValue('minDistanceToROIBoundary', 1.25, @(x) (isnumeric(x) && isscalar(x)) );
    
    p.addParamValue('maxFociCount', 4, @(x) (isnumeric(x) && isscalar(x)));
    p.addParamValue('metaInfoStruct', [], @(x) (isstruct(x) && all(isfield(x, {'header', 'data'}))) );
    p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue('flagSaveImages', false, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'finishStatusReportFile', [], @(x) (ischar(x)) );
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
    if ~exist( imageDataFilePath, 'file' )
        error( 'unable to find the image data file - %s', imageDataFilePath );
    end

    if ~exist( regionMergingModelFile, 'file' )
        error( 'unable to find region merging model file - %s', regionMergingModelFile );
    end

    PARAMETERS
    imageDataFilePath
    regionMergingModelFile
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
    PrettyPrintStepDescription( 'Loading Image Data' );
    
    imageSeries = loadIntravitalDataset( imageDataFilePath );    
    imageData = imageSeries(1).imageData;
    metadata = imageSeries(1).metadata

    compInfo.dataLoadTime = toc(dataLoadTimer);    
    
    if PARAMETERS.flagSaveImages

        imwrite( generateMultichannelMIPImage(imageData, [], metadata.voxelSpacing), ...
                 fullfile(resultsDir, 'images', 'stackMIP.png'), 'png');

    end
    
    % analyze data
    analysisTimer = tic;
    
        % segment nuclei
        PrettyPrintStepDescription( 'Running the nuclei segmentation algorithm' );
        
        segTimer = tic;
        
        [imLabelCellSeg, imCellSeedPoints, ...
         segAlgoParameters ] = segmentCellsInIntravitalData( imageData{channelId53BP1}, ...
                                                             metadata.voxelSpacing, ...                                                                      
                                                             'thresholdingAlgorithm', 'BackgroudRemovalUsingMorphologicalOpening', ...
                                                             'seedPointDetectionAlgorithm', 'AdaptiveMultiscaleLoG', ...
                                                             'cellDiameterRange', PARAMETERS.cellDiameterRange, ...
                                                             'minCellVolume', PARAMETERS.minCellVolume, ...
                                                             'flagIgnoreCellsOnXYBorder', PARAMETERS.flagIgnoreCellsOnXYBorder, ...
                                                             'regionMergingModelFile', regionMergingModelFile, ...
                                                             'flagParallelize', PARAMETERS.flagParallelize, ...
                                                             'flagDebugMode', false);

        
        segAlgoParameters

        numCells = max(imLabelCellSeg(:));       
        cellStats = regionprops( imLabelCellSeg, 'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
        
        fprintf( '\nThe segmentation algorithm found %d cells\n', numCells );
        
        compInfo.segmentationTime = toc(segTimer);        
        stackInfoStruct.cellCount = numCells;
        stackInfoStruct.cellDensity = sum([cellStats.Area]) / numel(imLabelCellSeg);
        
        % detect foci
        PrettyPrintStepDescription( 'Running the foci detection algorithm' );
        
        fociDetectionTimer = tic;
        
        [fociStats, imFociSeedPoints, imLabelFociSeg, ...
         fociDetectionParameters ] = segmentFociInsideNuclei( imageData{channelId53BP1}, ...
                                                              PARAMETERS.fociDiameterRange, ...                      
                                                              'spacing', metadata.voxelSpacing, ...
                                                              'roiMask', imLabelCellSeg, ...
                                                              'minDistanceToROIBoundary', PARAMETERS.minDistanceToROIBoundary);

        fociDetectionParameters
        
        for cid = 1:numel(cellStats)

            curCellPixInd = find(imLabelCellSeg == cid);
            curCellFoci = unique(imFociSeedPoints(curCellPixInd));
            curCellFoci = curCellFoci(curCellFoci > 0);
            cellStats(cid).foci = curCellFoci;
            cellStats(cid).fociCount = numel(curCellFoci);

        end        

        compInfo.fociDetectionTime = toc(fociDetectionTimer);        
        fprintf( '\nThe foci detection algorithm found %d foci\n', numel(fociStats));
        tabulate( [cellStats.fociCount] )
        
        % compute and store stack level info
        stackInfoStruct.totalFociCount = numel(fociStats);

        for i = 0:PARAMETERS.maxFociCount
            curField = (sprintf('cellFociCountDistribution_%.2d', i));
            stackInfoStruct.(curField) = 0;
        end
        
        for i = 1:numel(cellStats)
            
            if cellStats(i).fociCount >= PARAMETERS.maxFociCount
                curField = (sprintf('cellFociCountDistribution_%.2d', PARAMETERS.maxFociCount));
            else
                curField = (sprintf('cellFociCountDistribution_%.2d', cellStats(i).fociCount));
            end
            
            stackInfoStruct.(curField) = stackInfoStruct.(curField) + (1.0 / (eps + numel(fociStats)));
            
        end

        % compute and store cell level info
        fprintf('\n>> Computing some information about each cell ...\n');       
        
        cellInfoStruct = [];
        cellInfoFeatureMatrix = {};
        
        fprintf( '\nProgress: \n' );
        last_percent_done = 0;
        numPrint = 0;
        
        for cellId = 1:numel(cellStats)
            
            curCellStruct.cellId = cellId;
            
            % compute region properties
            curCellProps = ComputeRegionProperties( imLabelCellSeg, cellId );
            
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
            
            % some basic properties
            basicProps.centroidImsp = cellStats(cellId).Centroid;
            basicProps.centroidPhysp = cellStats(cellId).Centroid .* metadata.voxelSpacing;

            basicProps.bboxSizeImsp = cellStats(cellId).BoundingBox(4:end);
            basicProps.bboxSizePhysp = cellStats(cellId).BoundingBox(4:end) .* metadata.voxelSpacing; 

            volVoxel = prod(metadata.voxelSpacing);
            basicProps.volume = curCellProps.Area * volVoxel;       
            basicProps.volumeOfConvexHull = curCellProps.ConvexArea * volVoxel;
            basicProps.convexity = curCellProps.Area / curCellProps.ConvexArea;

            ptCell = ind2submat( size(imLabelCellSeg), curCellProps.PixelIdxList ); 
            ptCell = bsxfun(@times, ptCell, metadata.voxelSpacing);
            ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
            [U, S, V] = svd( (ptCell' * ptCell) / size(ptCell,1) );
            curCellEllipsoidRadiusPhysp = zeros(1,3); 
            for j = 1:3
                curCellEllipsoidRadiusPhysp(j) = 2 * sqrt(S(j,j)); % eigen-values are a measure of variance
            end
            basicProps.ellipsoidRadiusPhysp = curCellEllipsoidRadiusPhysp;        
            basicProps.fittedEllipsoidVolume = (4/3) * pi * prod(curCellEllipsoidRadiusPhysp);
            basicProps.ellipticVariance = ComputeEllipticVariance(ptCell, metadata.voxelSpacing);

            curCellStruct.cellprops = basicProps;      
            
            % foci inside cell
            curCellStruct.fociprops.cellFociCount = cellStats(cellId).fociCount;
            
            curCellStruct.fociprops.fociVolumeFraction = 0;
            
            curCellStruct.fociprops.meanCellFociRadius = 0;
            curCellStruct.fociprops.stdCellFociRadius = 0;
            curCellStruct.fociprops.minCellFociRadius = 0;            
            curCellStruct.fociprops.maxCellFociRadius = 0;           
            
            if cellStats(cellId).fociCount > 0
                curCellFociStats = fociStats(cellStats(cellId).foci);
                curCellStruct.fociprops.fociVolumeFraction = sum((4/3) * pi * [curCellFociStats.Radius].^3) / basicProps.volume;
                
                curCellStruct.fociprops.meanCellFociRadius = min([curCellFociStats.Radius]);
                curCellStruct.fociprops.stdCellFociRadius = max([curCellFociStats.Radius]);
                curCellStruct.fociprops.minCellFociRadius = mean([curCellFociStats.Radius]);
                curCellStruct.fociprops.maxCellFociRadius = std([curCellFociStats.Radius]);
            end
            
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
            
        end
        
        % compute and store foci level info
        fociInfoStruct = [];
        fociInfoFeatureMatrix = {};
        
        for fid = 1:numel(fociStats)                    
            fociInfoStruct(fid).cellId = imLabelCellSeg( fociStats(fid).PixelLocationIndex );
            fociInfoStruct(fid).cellFociCount = cellStats(fociInfoStruct(fid).cellId).fociCount;
            fociInfoStruct(fid).fociprops = fociStats(fid);
            [featureVec, fociInfoFeatureNameList] = ConvertFeatureStructToFeatureVec(fociInfoStruct(fid));
             fociInfoFeatureMatrix = cat(1, fociInfoFeatureMatrix, featureVec);
        end

        % save foci distribution plot
        numCells = numel(cellStats); 
        punctaCounts = [cellStats.fociCount];

        hFociDist = figure();
        [counts, centers] = histogram(punctaCounts, 'discrete');
        bar(centers, counts/sum(counts), ...
            'linewidth', 2.0, ...
            'barwidth', 0.9, ...
            'EdgeColor', 'k', ...
            'FaceColor', 0.7 * ones(1,3));

        [~, fname, ~] = fileparts(imageDataFilePath);
        fname = strrep(fname,'_','\_');
        text(0.5, 0.95, {fname, sprintf('#cells = %d, #puncta= %d', numCells, sum(punctaCounts))}, ...
             'units', 'normalized', 'horizontalalignment', 'center');

        title({'Distribution of Number of Puncta In a Cell'}, ...
              'FontSize', 12.0, 'FontWeight', 'bold');
        xlabel('Number of puncta per cell');
        ylabel('Proportion of cells');
        ylim([0, 1]);
        grid on;
		Maximize(hFociDist);
        SaveFigure(hFociDist , fullfile(resultsDir, 'fociCountDistributionPlot.png') , 'png');
        SaveFigure(hFociDist , fullfile(resultsDir, 'fociCountDistributionPlot.fig') , 'fig');
        close(hFociDist);
            
        % save images
        if PARAMETERS.flagSaveImages 

            fprintf('\n>> Saving snapshot images of cells ...\n');       
            
            imageSaveTimer = tic;

            % save cell images
            if PARAMETERS.flagParallelize

                parfor cellId = 1:numCells

                    if(cellStats(cellId).fociCount == 0)
                        continue;
                    end
                    
                    imageOutputDir = fullfile(resultsDir, 'images', num2str(cellStats(cellId).fociCount));

                    if ~isdir( imageOutputDir )
                        mkdir( imageOutputDir );
                    end

                    WriteCellSnapshotImages(imageData, imLabelCellSeg == cellId, ...
                                            cellStats(cellId), cellId, fociStats(cellStats(cellId)), ...
                                            metadata.voxelSpacing, imageOutputDir)

                end

            else

                fprintf( '\nProgress: \n' );

                last_percent_done = 0;
                numPrint = 0;

                for cellId = 1:numCells

                    if(cellStats(cellId).fociCount == 0)
                        continue;
                    end
                    
                    imageOutputDir = fullfile(resultsDir, 'images', num2str(cellStats(cellId).fociCount));

                    if ~isdir( imageOutputDir )
                        mkdir( imageOutputDir );
                    end

                    WriteCellSnapshotImages(imageData, imLabelCellSeg == cellId, ...
                                            cellStats(cellId), cellId, fociStats(cellStats(cellId)), ...
                                            metadata.voxelSpacing, imageOutputDir)
                    
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

            compInfo.imageSaveTime = toc(imageSaveTimer);

        end        
        
    compInfo.totalAlgorithmTime = compInfo.segmentationTime + compInfo.fociDetectionTime;
    compInfo.totalAnalysisTime = toc(analysisTimer);

    stackInfoStruct.computation = compInfo;
    stackInfoStruct.resultsDir = resultsDir;
    
    metaInfoStruct
    stackInfoStruct
    compInfo
        
    % write stack info
    fprintf( '\n>>Writing stack analysis csv file ... \n' );

    [ stackFeatureVec , stackFeatureNameList ] = ConvertFeatureStructToFeatureVec( stackInfoStruct );

    if ~isempty( PARAMETERS.metaInfoStruct )
        stackFeatureNameList = cat(2, metaInfoStruct.header, stackFeatureNameList);            
        stackFeatureVec = cat(2, metaInfoStruct.data , stackFeatureVec);            
    end

    WriteFeatureMatrixToCSVFile( fullfile(resultsDir, 'stackAnalysisInfo.csv' ), ...
                                 stackFeatureVec, stackFeatureNameList );
                             
    % write cell info
    fprintf( '\n>>Writing cell analysis csv file ... \n' );

    featureNameList = cat(2, stackFeatureNameList, cellInfoFeatureNameList );
    featureMatrix = cat(2, repmat(stackFeatureVec, numCells, 1), cellInfoFeatureMatrix );

    WriteFeatureMatrixToCSVFile( fullfile(resultsDir, 'cellAnalysisInfo.csv' ), ...
                                 featureMatrix, featureNameList );

    % write foci info
    fprintf( '\n>>Writing foci analysis csv file ... \n' );

    featureNameList = cat(2, stackFeatureNameList, fociInfoFeatureNameList );
    featureMatrix = cat(2, repmat(stackFeatureVec, numel(fociInfoStruct), 1), fociInfoFeatureMatrix);

    WriteFeatureMatrixToCSVFile( fullfile(resultsDir, 'fociAnalysisInfo.csv' ), ...
                                 featureMatrix, featureNameList );
                             
    % save mat file that can be loaded into the DNADamageAnalyzer tool
    analysisData.dataFilePath = imageDataFilePath;    
    analysisData.imageData = imageData;
    analysisData.metadata = metadata;
    analysisData.channelId53BP1 = channelId53BP1;

    analysisData.imLabelCellSeg = imLabelCellSeg;
    analysisData.imCellSeedPoints = imCellSeedPoints;
    [imCellSegMaskRGB, analysisData.CellSegColorMap] = label2rgbND(imLabelCellSeg);
    analysisData.segAlgoParameters = segAlgoParameters;
    analysisData.cellStats = cellStats;
    
    analysisData.fociStats = fociStats;
    analysisData.imFociSeedPoints = imFociSeedPoints;
    analysisData.imLabelFociSeg = imLabelFociSeg;
    [imFociSegMaskRGB, analysisData.FociSegColorMap] = label2rgbND(imLabelFociSeg);
    analysisData.fociDetectionParameters = fociDetectionParameters;
    
    analysisData.regionMergingModelFile = regionMergingModelFile;

    save( fullfile(resultsDir, 'DNADamageAnalysis.mat'), '-struct', 'analysisData' );
                             
    % report success
    flagSuccess = true;
    fprintf( '\n\n---- Analysis Succeeded ---- \n\n' );
    
    if ~isempty(PARAMETERS.finishStatusReportFile)
        fprintf( fidStatus, 'Analysis Succeeded' );
        fclose(fidStatus);
    end
    
    diary off;
    
end

function WriteCellSnapshotImages(imageData, imCellMask, cellStats, cellId, fociStats, spacing, imageOutputDir)

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
    subinds{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));

    imCurCellMIP = mat2gray( max(imageData(subinds{:}) .* imCellMask(subinds{:}), [], 3) );

    imLabelFociMask = zeros(size(imCurCellMIP));    
    fociSeedPointLoc = ind2submat(size(imCellMask), [fociStats.PixelLocationIndex]);
    fociSeedPointLoc(:,3) = [];
    for i = 1:2
        fociSeedPointLoc(:,i) = fociSeedPointLoc(:,i) - min(subinds{i}) + 1;
    end        
    seedPos = fociSeedPointLoc * diag(spacing);
    kd = KDTreeSearcher(seedPos);
    pixelPos = ind2submat(size(imCurCellMIP), (1:numel(imCurCellMIP))') * diag(spacing);
    [closestSeedInd, distanceToSeed] = kd.knnsearch(pixelPos);
    flagIsPixelInSeedVicinity = abs(distanceToSeed - blobRadii(closestSeedInd)) <= min(1.5 * spacing);
    imLabelFociMask( flagIsPixelInSeedVicinity ) = closestSeedInd( flagIsPixelInSeedVicinity );
    imFociRGBMask = label2rgbND(imLabelFociMask);
    
    imwrite( genImageRGBMaskOverlay(imCurCellMIP, imFociRGBMask, 0.5 ), ...
             fullfile(imageOutputDir, sprintf('Cell53BP1MIPFociOverlay_%.3d.png', cellId)), 'png' );   

end