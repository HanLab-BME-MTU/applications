%function analyzeApoptosisDataset

clc
clear
close all
fclose all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    apoptosisDataRootDir = 'Z:\intravital\Lee-Data';
    
    elastixParameterFile = fullfile( 'C:\deepak\code\matlab\intravital', 'VideoStabilization_3D.eparam' );
    
    regionMergingModelFile = fullfile( 'Z:\intravital\data\image_analysis', ...
                                       'models', 'region_merging', 'models', ...
                                       'train_M04_M09_M12', ...
                                       'regionMerging.model' );                                      

    flagParallelize = true;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get data file
dataRootDir = apoptosisDataRootDir; 

if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( dataRootDir, '*.oif; *.nucseg' ), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
    [~,~,ext] = fileparts(fileName);
end

% turn on parallelization if requested
if flagParallelize
    flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if ~flagPoolOpenedAlready 
        matlabpool open;
    end            
end

%% Load data
if ~strcmp(ext, '.nucseg')
    
    imageSeries = loadIntravitalDataset( dataFilePath );   

    metadata = imageSeries.metadata;
    imageData = imageSeries.imageData;
    clear imageSeries;

    %% Segment nuclei -- will have to be parallelized on the cluster
    PrettyPrintStepDescription( 'Segmenting Nuclei in all frames' );
    
    cellDiameterRange = [8, 20];
    minCellVolume = 400;
    minCellROIOverlap = 0.5;

    nucleiDetectionResult = cell(metadata.numTimePoints, 2);

    segFunc = @(im) (segmentCellsInIntravitalData( im, metadata.voxelSpacing, ...                                                                      
                                                   'cellDiameterRange', cellDiameterRange, ...
                                                   'minCellVolume', minCellVolume, ...
                                                   'minCellROIOverlap', minCellROIOverlap, ...
                                                   'regionMergingModelFile', regionMergingModelFile, ...                                                           
                                                   'flagParallelize', true, ...
                                                   'flagDebugMode', false) );
                        
    segTimer = tic;

    AddWekaClassesToPath();
    
    for tid = 1:metadata.numTimePoints
        
        curResult = cell(1,2);
        curImageData = imageData(tid,:);
        
        for cid = 1:2

            curSegTimer = tic;
            
            [ cmdOutput, imLabelCellSeg, imCellSeedPoints ] = evalNucleiSegmentation(segFunc, curImageData{cid});
            
            curResult{cid}.imLabelCellSeg = uint16( imLabelCellSeg );
            curResult{cid}.imCellSeedPoints = logical( imCellSeedPoints );
            curResult{cid}.numCells = max(imLabelCellSeg(:));               

            fprintf( '\nSegmented image in frame %d/%d and channel %d/2 ... took %f seconds\n', ...
                      tid, metadata.numTimePoints, cid, toc(curSegTimer) );
        end
        
        nucleiDetectionResult(tid, :) = curResult;
        
    end

    fprintf( '\nSegementation of %d stacks took a total of %f seconds ... \n', ...
             metadata.numTimePoints, toc(segTimer) );
    
else
    
    f = load( dataFilePath, '-mat');
    f = f.dataStruct;
    imageData = f.imageData;
    metadata = f.metadata;
    imageFilePath = f.dataFilePath;
    nucleiDetectionResult = f.nucleiDetectionResult;
    clear f;
    
end

%% Compute global motion by registering consecutive frames of the vascular channel
PrettyPrintStepDescription( 'Registering consecutive frames of vascular channel to estimate global motion' );

elastixParameters{1} = readElastixParameters( elastixParameterFile );
vesselChannelData = imageData(:,3);

frame2frameElastixTransform = cell(metadata.numTimePoints-1, 3);
frame2frameGlobalTransfrom = zeros(metadata.numTimePoints-1, 3);

parfor t = 1:metadata.numTimePoints-1
    
    regTimer = tic;
    
    [imRegistered, regTransform, imValidROI] = registerImagesUsingElastix(vesselChannelData{t+1}, vesselChannelData{t}, elastixParameters, metadata.voxelSpacing);
    frame2frameGlobalTransfrom(t, :) = regTransform{1}.TransformParameters;
    frame2frameElastixTransform{t} = regTransform;
    
    fprintf( '\nRegistered frame %d to %d ... took %d seconds\n', t, t+1, toc(regTimer) );
    
end

clear vesselChannelData;
globalMotionTrajectory = [ zeros(1,3); frame2frameGlobalTransfrom ];

figure;
hold all;
    plot( globalMotionTrajectory(:,1), '.-' );
    plot( globalMotionTrajectory(:,2), '.-' );
    plot( globalMotionTrajectory(:,3), '.-' );
    legend( {'X', 'Y', 'Z'} );
    xlabel('Time-Frame');
    title('1D Projections of the Global Motion Trajectory', 'FontWeight', 'bold' );
    grid on;
hold off;


%% Track nuclei using kalman filters, interacting motion model, and LAP
PrettyPrintStepDescription( 'Tracking nuclei using kalman filters, IMM, and LAP' );

    % specify parameters for tracking
    gapCloseParam.timeWindow = 4;
    gapCloseParam.mergeSplit = 3;
    gapCloseParam.minTrackLen = 3;
    gapCloseParam.diagnostics = 0;
    
    costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
    plinkparams.linearMotion = 1;
    plinkparams.minSearchRadius = 2 * max(metadata.voxelSpacing);
    plinkparams.maxSearchRadius = 4 * max(metadata.voxelSpacing);
    plinkparams.brownStdMult = 3;
    plinkparams.useLocalDensity = 1;
    plinkparams.nnWindow = gapCloseParam.timeWindow;
    costMatrices(1).parameters = plinkparams;

    costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';
    tlinkparams = plinkparams;
    tlinkparams.brownStdMult = plinkparams.brownStdMult * ones(gapCloseParam.timeWindow, 1);
    tlinkparams.brownScaling = [0.25, 0.01];
    tlinkparams.timeReachConfB = gapCloseParam.timeWindow;
    tlinkparams.ampRatioLimit = [0.7, 4];
    tlinkparams.lenForClassify = 4;
    tlinkparams.useLocalDensity = 0;
    tlinkparams.linStdMult = 1*ones(gapCloseParam.timeWindow,1);
    tlinkparams.linScaling = tlinkparams.brownScaling;
    tlinkparams.timeReachConfL = tlinkparams.timeReachConfB;
    tlinkparams.maxAngleVV = 30;
    tlinkparams.gapPenalty = 1;
    tlinkparams.resLimit = [];
    costMatrices(2).parameters = tlinkparams;
    
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

    % get object measurements in each time frame
    movieInfo = cell(1,2);
    
    for cid = 1:2
        
        offsetGlobal = zeros(1,3);
        
        for tid = 1:metadata.numTimePoints

            cellprops = regionprops( nucleiDetectionResult{tid,cid}.imLabelCellSeg, ...
                                     'Centroid', 'Area', 'PixelIdxList' );
        
            for i = 1:numel(cellprops)
                
                movieInfo{cid}(tid).xCoord(i,1) = offsetGlobal(1) + cellprops(i).Centroid(1) * metadata.voxelSpacing(1);
                movieInfo{cid}(tid).yCoord(i,1) = offsetGlobal(2) + cellprops(i).Centroid(2) * metadata.voxelSpacing(2);
                movieInfo{cid}(tid).zCoord(i,1) = offsetGlobal(3) + cellprops(i).Centroid(3) * metadata.voxelSpacing(3);
                movieInfo{cid}(tid).amp(i,1) = median( double(imageData{tid,cid}(cellprops(i).PixelIdxList)) );
                
                movieInfo{cid}(tid).xCoord(:,2) = 0;
                movieInfo{cid}(tid).yCoord(:,2) = 0;
                movieInfo{cid}(tid).zCoord(:,2) = 0;
                movieInfo{cid}(tid).amp(:,2) = 0;
                
            end

            if tid < metadata.numTimePoints
                offsetGlobal = offsetGlobal + frame2frameGlobalTransfrom(tid,:);
            end
            
        end
    end
    
    % track nuclei
    trackInfo = cell(1,2);
    
    for cid = 1:2
        
        [trackInfo{cid}.tracksFinal, ...
         trackInfo{cid}.kalmanInfoLink, ...
         trackInfo{cid}.errFlag] = trackCloseGapsKalmanSparse( ...
                                            movieInfo{cid}, ...
                                            costMatrices, ...
                                            gapCloseParam, ...
                                            kalmanFunctions, ...
                                            3, 0, 1);
    
    end
    
%% Display Tracking Result in Imaris
PrettyPrintStepDescription( 'Displaying tracking result in imaris' );
surfaceQuality = 0.05;

    vis = ImarisDataVisualizer(imageData, 'spacing', metadata.voxelSpacing);    

    objectLocInfo = cell(1,2);
    flagIsApoptotic = cell(1,2);
    
    for cid = 1:2

        % generate isosurfaces
        objectGeometry = cell(1, metadata.numTimePoints);
        
        h = waitbar(0, sprintf('Generating isosurfaces of objects in channel %d for rendering in imaris', cid));
        
        for tid = 1:metadata.numTimePoints
            
            curDetectionResult = nucleiDetectionResult{tid,cid};
            curObjectGeometry = cell(curDetectionResult.numCells, 1);
            
            parfor obid = 1:curDetectionResult.numCells
                imObjectMask = (curDetectionResult.imLabelCellSeg == obid);
                curObjectGeometry{obid} = ImarisDataVisualizer.generateSurfaceFromMask(imObjectMask, 'timepoint', tid, 'surfaceQuality, 'surfaceQuality);
            end
        
            objectGeometry{tid} = curObjectGeometry;
            waitbar(tid / metadata.numTimePoints, h);
            
        end
        
        close(h);
        
        % display tracks
        numTracks = numel( trackInfo{cid}.tracksFinal );
        chgroup = vis.AddDataContainer();

        healthyGroup = vis.AddDataContainer( chgroup );
        sickGroup = vis.AddDataContainer( chgroup );

        h = waitbar(0, sprintf('Rendering tracking result for channel %d', cid));

        flagIsApoptotic{cid} = false(numTracks, 1);

        for trackid = 1:numTracks

            curTrack = trackInfo{cid}.tracksFinal(trackid);

            % get start-frame, end-frame and lifetime of each track segment
            numTrackSegments = size(curTrack.tracksFeatIndxCG, 1);
            segStartEndLifetime = getTrackSEL( curTrack, 1 );

            % determine type of track -- interphase, apoptosis
            flagIsDividing = numTrackSegments > 1;
            flagIsApoptotic{cid}(trackid) = any( segStartEndLifetime(:,2) < metadata.numTimePoints );

            if flagIsApoptotic{cid}(trackid)
                curTrackGroup = sickGroup;
            else
                curTrackGroup = healthyGroup;
            end

            % generate and display surfaces of objects in each segment
            curTrackSurfaceList = cell(numTrackSegments, 1);
            curTrackEdges = [];

            for segid = 1:numTrackSegments

                % get the frames in which each object resides
                curObjectFrames = segStartEndLifetime(segid,1):segStartEndLifetime(segid,2);

                % get local indices i.e. indices of objects in their frame
                curObjectIndices = curTrack.tracksFeatIndxCG(segid, :);                
                validIndices = find(curObjectIndices);
                validIndices = min(validIndices):max(validIndices);
                curObjectIndices = curObjectIndices(validIndices);
                
                % ignore any zero indices
                curObjectFrames(curObjectIndices == 0) = [];
                curObjectIndices(curObjectIndices == 0) = [];

                assert( numel(curObjectIndices) == numel(curObjectFrames) );
                
                % generate isosurfaces for each object
                curSegmentSurfaceList = cell(numel(curObjectIndices), 1);

                for obid = 1:numel(curObjectIndices)
                    curSegmentSurfaceList{obid} = objectGeometry{curObjectFrames(obid)}{ curObjectIndices(obid) };
                end

                curSegmentTrackEdges = size(curTrackEdges, 1) + ([1:numel(curObjectIndices)-1; 2:numel(curObjectIndices)])';
                curTrackEdges = cat(1, curSegmentTrackEdges);

                curTrackSurfaceList{segid} = cat(1, curSegmentSurfaceList{:});

            end

            curTrackSurfaceList = cat(1, curTrackSurfaceList{:});

            vis.AddSurfaces( curTrackSurfaceList, curTrackGroup, ...
                            'tracks', curTrackEdges, ...    
                            'name', sprintf('Track_%d', trackid) );

            waitbar(trackid / numTracks, h);

        end

        close(h);

        fprintf( '\n%.2d (%d %%) apoptotic candidates were found in channel %d\n', ...
                 sum(flagIsApoptotic{cid}), mean(flagIsApoptotic{cid}) * 100, cid );

    end
    
