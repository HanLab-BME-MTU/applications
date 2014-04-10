
    clc
    clear 
    close all

    % load data
    stefanOldDataDir = 'C:\deepak\data\Stefan_June2012';
    stefanNewDataDir = 'C:\deepak\data\Stefan_July2012';
    ralphDataDir = 'C:\deepak\data\Ralph_June2012';
    adhesionDir = 'C:\deepak\data\adhesions';

    dataFileDir = stefanNewDataDir;

    if ~exist( 'dataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.mat' ), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    dataFileContents = load( dataFilePath );
    imInput = dataFileContents.im;

    % adjust the intensity range of the image 
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = medfilt2( imAdjusted, [3,3] );    

    % compute threshold
    fprintf( '\nThresholding the Image using a Global Thresholding Algorithm ...\n' );
    numHistogramBins = 256; 
    thresholdingAlgorithm = 'FOtsuThreshold';
    imLOG = mat2gray( ComputeImageLogTransform( imAdjusted ) ); 
    imThresh = callMatitkFilter( thresholdingAlgorithm, { numHistogramBins } , imLOG );                        
    
    imseriesmaskshow( imAdjusted, imThresh );
    set( gcf, 'Name', 'Result of Thresholding' );
    
    % compare different methods for cells seed point detection   
    meanCellDiameter = 10;    
    cellDiameterRange = [8, 12];
    spacing = [0.5, 0.5];
    imCellSeedDetectionResult = {};
    strSeedDetectionAlgorithm = {};
    
        % Local Intensity Maxima in Blurred Image
        fprintf( '\n\nDetecting seed points as local maxima in Blurred Image ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'IntensityMaxima';
        imCellSeedDetectionResult{end+1} = detect_cell_seeds_IntensityMaxima( imAdjusted, ...
                                                                              meanCellDiameter, ...
                                                                              'spacing', spacing, ...
                                                                              'debugMode', true );
                                                        
        % Local Maxima in LoG Filtered Image
        fprintf( '\n\nDetecting seed points as local maxima in LoG Filtered Image ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'FixedScaleLoG';
        imCellSeedDetectionResult{end+1} = detectBlobsUsingLoG( imAdjusted, ...
                                                                meanCellDiameter, ...
                                                                'spacing', spacing, ...
                                                                'debugMode', true );

        % Local Maxima in a Multiscale LoG Filtered Image
        fprintf( '\n\nDetecting seed points as local maxima in Multiscale LoG Filtered Image ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'MultiScaleLoG';        
        imCellSeedDetectionResult{end+1} = detectBlobsUsingMultiscaleLoG( imAdjusted, ...
                                                                          cellDiameterRange, ...
                                                                          'spacing', spacing, ...
                                                                          'debugMode', true );

        % Local Maxima in a Multiscale LoG Filtered Image
        fprintf( '\n\nDetecting seed points as local maxima in Adaptive Multiscale LoG Filtered Image ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'AdaptiveMultiScaleLoG';
        imCellSeedDetectionResult{end+1} = detect_cell_seeds_adaptive_multiscale_LoG( imAdjusted, ...
                                                                                      imThresh, ...
                                                                                      'cellDiameterRange', cellDiameterRange, ...
                                                                                      'spacing', spacing, ...
                                                                                      'debugMode', true );

        % Local Maxima in a Multiscale LoBG Filtered Image
        fprintf( '\n\nDetecting seed points as local maxima in Adaptive Multiscale LoG Filtered Image ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'MultiScaleLoBG';
        imCellSeedDetectionResult{end+1} = detectBlobsUsingMultiscaleLoBG( imAdjusted, ...
                                                                           cellDiameterRange, ...
                                                                           'spacing', spacing, ...
                                                                           'debugMode', true );

        % Local Maxima in a radial symmetry transform
        fprintf( '\n\nDetecting seed points as local maxima in multiscale radial symmetry transform ...\n\n' );
        
        strSeedDetectionAlgorithm{end+1} = 'RadialSymmetryTransform';
        roiMask = imdilate( imThresh, ones(3,3) );
        imCellSeedDetectionResult{end+1} = detect_cell_seeds_radial_symmetry( imAdjusted, ...
                                                                              0.5 * cellDiameterRange, ...
                                                                              'roiMask', roiMask, ...
                                                                              'spacing', spacing, ...
                                                                              'flagParallelize', false, ...
                                                                              'debugMode', true );
                                                                          
        % suppress seed points outside the thresholded foreground region
        for i = 1:numel( imCellSeedDetectionResult )
            imCellSeedDetectionResult{i}( ~imThresh ) = 0;
        end
        
        % display result
        imSeedMaskForDisplay = imCellSeedDetectionResult;
        for i = 1:numel( imSeedMaskForDisplay )
            imSeedMaskForDisplay{i} = imdilate( imSeedMaskForDisplay{i}, strel('disk', 3) );
        end
        imseriesmaskshow( imAdjusted, imSeedMaskForDisplay );
        set( gcf, 'Name', [ 'Cell Seed Point Detection Results', ...
                            sprintf( ' -- %s', strSeedDetectionAlgorithm{:} ) ] );
