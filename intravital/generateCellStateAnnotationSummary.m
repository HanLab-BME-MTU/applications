clc;
clear;
fclose all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultDataRootDir = 'Z:\intravital\data';    
    defaultOutputDir = 'Z:\intravital\data';    

    PARAMETERS.flagInputRootDirListFile = true; 
    
    PARAMETERS.flagIgnoreBadlyDetectedCells = false;
    PARAMETERS.flagIgnoreBadlySegmentatedCells = false;
    PARAMETERS.flagSaveImages = false;
    
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
    rootDirList = uipickfiles('FilterSpec', defaultDataRootDir);
    
end

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% for parallelization
if matlabpool( 'size' ) == 0
    matlabpool open;
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
rootDirList
outputRootDir

fprintf( '\nList of %d annotation files that will be processed :\n', numel(annotationFileList));
for i = 1:numel(annotationFileList)
    fprintf( '\n%d/%d - %s\n', i, numel(annotationFileList), annotationFileList(i).name );    
end

if PARAMETERS.flagSaveImages
    if isdir( fullfile(outputRootDir, 'images') )            
        rmdir( fullfile(outputRootDir, 'images'), 's' );
    end
end

% process each annotation file
annotatedCellPatternList = {};
annotatedCellInfoData = {};
stackInfoData = {};

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
   
   imageDataAdjusted = CellCycleStateClassifier.preprocessImageData( annotationData.imageData );
   
   % compute features for each cell
   numCells = numel(annotationData.cellStats);
   
   annotatedCellInfoDataLocal = {};
   annotatedCellPatternListLocal = cell(numCells, 1);
   
   for cellId = 1:numCells
       
        curCellStats = annotationData.cellStats(cellId);
        
        if strcmpi( curCellStats.cellPatternType, 'None' )
            error( 'Incomplete annotation. File %s contains unannotated cells', annotationFilePath );
        end
        
        annotatedCellPatternListLocal{cellId} = curCellStats.cellPatternType;
        
        % compute basic region properties
        cellProps = ComputeRegionProperties( imLabelCellSeg, cellId );
        
        % note down some cell identity information
        curCellInfoStruct.annotationFilePath = curAnnotationFilePath;
        curCellInfoStruct.cellId = cellId;
        curCellInfoStruct.cellPatternType = curCellStats.cellPatternType;
        curCellInfoStruct.Centroid = curCellStats.Centroid;
        curCellInfoStruct.DepthPhysp = (curCellStats.Centroid(3)-1) .* metadata.voxelSpacing(3);
        
        volVoxel = prod(metadata.voxelSpacing);
        curCellInfoStruct.Volume = cellProps.Area * volVoxel;
        curCellInfoStruct.BBoxSize = curCellStats.BoundingBox(4:end) .* metadata.voxelSpacing;
        curCellInfoStruct.VolumeOfConvexHull = cellProps.ConvexArea * volVoxel;
        curCellInfoStruct.Convexity = cellProps.Area / cellProps.ConvexArea;
                
        ptCell = ind2submat( size(imLabelCellSeg), cellProps.PixelIdxList ); 
        ptCell = bsxfun(@times, ptCell, metadata.voxelSpacing);
        ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
        [U, S, V] = svd( (ptCell' * ptCell) / size(ptCell,1) );
        curCellEllipsoidRadiusPhysp = zeros(1,3); 
        for j = 1:3
            curCellEllipsoidRadiusPhysp = 2 * sqrt(S(j,j)); % eigen-values are a measure of variance
        end
        curCellInfoStruct.ellipsoidRadiusPhysp = curCellEllipsoidRadiusPhysp;        
        curCellInfoStruct.ellipticVariance = ComputeEllipticVariance( ptCell, metadata.voxelSpacing );
        
        [curCellInfoVec , cellInfoLabelList] = ConvertFeatureStructToFeatureVec(curCellInfoStruct);        

        annotatedCellInfoDataLocal = cat(1, annotatedCellInfoDataLocal, curCellInfoVec);
        
        % save images if requested
        if PARAMETERS.flagSaveImages 
            
            [afpathstr, afname, ~] = fileparts( annotationData.dataFilePath{1} );
            curImageOutputDir = fullfile(outputRootDir, 'images', curCellStats.cellPatternType, sprintf('fid_%.2d_%s', fid, strtrim(afname)) );

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

            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
                     fullfile(curImageOutputDir, sprintf('CellMidSliceAllChannel_%.3d.png', cellId)), 'png' );   
                 
            imwrite(imCurCellMIP, fullfile(curImageOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   
                 
        end        
        
   end

   fprintf( '\nannotated cell pattern distribution: \n' );
   tabulate( annotatedCellPatternListLocal )   
   annotatedCellPatternList = cat(1, annotatedCellPatternList, annotatedCellPatternListLocal );
   
   fprintf( '\nTotal Cell Count: %d\n', numCells );
   
   % basic stack info
   curStackInfoStruct.annotationFilePath = curAnnotationFilePath; 
   curStackInfoStruct.volSize = metadata.volSize;
   curStackInfoStruct.voxelSize = metadata.voxelSpacing;
   curStackInfoStruct.stackVolumePhysp = prod(metadata.volSize) * prod(metadata.voxelSpacing);
   
   curStackInfoStruct.CellCount = numCells;
   
   assert( numel(annotatedCellPatternListLocal) == numCells );
   
   % segmentation quality
   curSegPerformance.BadDetectionCount = sum( strcmp( annotatedCellPatternListLocal, 'Bad_Detection' ) );   
   curSegPerformance.UnderSegmentationCount = sum( strcmp( annotatedCellPatternListLocal, 'Under_Segmentation' ) );
   curSegPerformance.OverSegmentationCount = sum( strcmp( annotatedCellPatternListLocal, 'Over_Segmentation' ) );   
   curSegPerformance.TotalErrorCount = sum(ismember(annotatedCellPatternListLocal, {'Bad_Detection', 'Under_Segmentation', 'Over_Segmentation'}));
   curSegPerformance.WellSegmentedCellCount = sum(~ismember(annotatedCellPatternListLocal, {'Bad_Detection', 'Under_Segmentation', 'Over_Segmentation'}));
   
   curSegPerformance.BadDetectionPercentage = 100.0 * curSegPerformance.BadDetectionCount / numCells;
   curSegPerformance.OverSegmentationPercentage = 100.0 * curSegPerformance.OverSegmentationCount / numCells;
   curSegPerformance.UnderSegmentationPercentage = 100.0 * curSegPerformance.UnderSegmentationCount / numCells;
   curSegPerformance.TotalErrorPercentage = 100.0 * curSegPerformance.TotalErrorCount / numCells;
   
   curSegPerformance.segmentationAccuracy = 100.0 * curSegPerformance.WellSegmentedCellCount / numCells;

   curSegPerformance
   
   curStackInfoStruct.segQuality = curSegPerformance;
   
   fprintf( '\nSegmentation_Accuracy: %.2f%%\n', curSegPerformance.segmentationAccuracy );
   
   % estimate cell density
   fprintf( '\nRunning a thresholding algorithm to get rough estimates of cell density ...\n' );   
   localThresholdWindowRadiusPhysp = 30;
   localWindowRadius = round(localThresholdWindowRadiusPhysp ./ metadata.voxelSpacing(1));
   localWindowPace = round(localWindowRadius / 3);
   minLocalGlobalThresholdRatio = 0.6;   
   imThresh = segmentCellForegroundUsingLocalMinError( imageDataAdjusted{1}, localWindowRadius, ...
                                                       'model', 'poisson', ...  
                                                       'localWindowPace', localWindowPace, ...
                                                       'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                       'flagParallelize', true );
   
   curCellDensity.foregroundVolumePhysp = sum( imThresh(:) > 0 ) * prod(metadata.voxelSpacing);
   curCellDensity.foregroundVolPercent = 100.0 * curCellDensity.foregroundVolumePhysp / curStackInfoStruct.stackVolumePhysp;
   
   fprintf( '\n\tStack contains %.2f%% foreground ...\n', curCellDensity.foregroundVolPercent );   
   
   avgCellVolume = (4/3) * pi * 6^3;
   curCellDensity.estimatedCellCount = curCellDensity.foregroundVolumePhysp / avgCellVolume;
   
   curCellDensity
   
   curStackInfoStruct.cellDensity = curCellDensity;
   
%    % estimate noise level of the image
%    fprintf( '\nEstimating the noise level of the image ...\n' );
%    
%    imStandardized = mat2gray( annotationData.imageData{1} ) * 4096;
%    
%    imNoiseFree = matitk( 'FMEDIAN', [1,1,1], imStandardized );    
%    imNoiseEstimate = imStandardized - imNoiseFree;
%    funcGlobalNoise = @(g,a) (g^2 * imNoiseFree(:) + a^2 - var(imimNoiseEstimate(:)))
%    
%    noiseWindowRad = 15;
%    noiseStrel = ones((2*noiseWindowRad+1)*[1,1]);
%    imNoiseVar = (eps + stdfilt( imStandardized, noiseStrel )).^2;
%    imBgnd = imerode( ~imThresh, noiseStrel );      
%    curImageNoise.signalIndependentNoise = median( imNoiseVar(imBgnd > 0) );
%    
%    
%    curImageNoise.medianSNR = median( imAvg(imBgnd > 0) ./ imNoiseVar(imBgnd > 0) );
%    
%    curImageNoise
%    
%    curStackInfoStruct.noiseLevel = curImageNoise;
   
   % store stack info
   [featureVec , stackInfoLabelList] = ConvertFeatureStructToFeatureVec( curStackInfoStruct );        
   stackInfoData = cat(1, stackInfoData, featureVec);   
   
   curStackInfoStruct
   
   computationTime = toc(feaCompTimer)   
    
end

fprintf( '\n\n*****************************************************\n\n' );

fprintf( '\nannotated cell pattern distribution in all datasets: \n' );
tabulate( annotatedCellPatternList )   

stackInfoFileName = 'stackAnnotationInfo.csv';
WriteFeatureMatrixToCSVFile( fullfile(pwd, stackInfoFileName), stackInfoData, stackInfoLabelList);
movefile( fullfile(pwd, stackInfoFileName), fullfile(outputRootDir, stackInfoFileName) );

cellInfoFileName = 'cellAnnotationInfo.csv';
WriteFeatureMatrixToCSVFile( fullfile(pwd, cellInfoFileName), annotatedCellInfoData, cellInfoLabelList);
movefile( fullfile(pwd, cellInfoFileName), fullfile(outputRootDir, cellInfoFileName) );

% switch off diary
diary off;
