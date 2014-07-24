classdef segmentCellsInIntravitalDataUsingRandomWalks
    
    properties (SetAccess = private)
        
        imInput;     
        spacing;
        
        cellDiameterRange = [12, 20];
        minCellVolume = 400;

        flagParallelize = false;
        flagUseGPU = false;        
        
        flagDebugMode = false;
        
    end

    methods
        
        function obj = segmentCellsInIntravitalDataUsingRandomWalks( imInput, spacing )
           
            obj.imInput = imInput;
            obj.spacing = spacing;
            
            obj.cellDiameterRange = [12, 20];
            obj.minCellVolume = 800;
            
            obj.flagParallelize = false;
            obj.flagUseGPU = false;
            
            obj.flagDebugMode = false;
            
        end        
        
    end
end

function [ imLabelCellSeg, varargout ] = segmentCellsInIntravitalDataUsingRandomWalks( imInput, spacing, varargin )

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
    p.parse( imInput );

    imdims = ndims(imInput);   
    
    p.addRequired( 'spacing', @(x) (numel(x) == imdims) );    
    p.addParamValue( 'thresholdingAlgorithm', ...
                     'MinErrorPoissonSliceBySliceLocal', ...
                     @(x) (ischar(x) && ismember(x, {'OtsuSliceBySliceLocal', 'MinErrorPoissonSliceBySliceLocal'}) ));
    p.addParamValue( 'seedPointDetectionAlgorithm', ...
                     'MultiscaleLoG', ...
                     @(x) (ischar(x) && ismember(x,{'IntensityMaxima', 'MultiscaleLoG', 'AdaptiveMultiscaleLoG'})) );
    p.addParamValue( 'cellDiameterRange', [12, 20], @(x) (isnumeric(x) && numel(x) == 2) );
    p.addParamValue( 'minCellVolume', 800, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, spacing, varargin{:} );    
    
    thresholdingAlgorithm = p.Results.thresholdingAlgorithm;
    seedPointDetectionAlgorithm = p.Results.seedPointDetectionAlgorithm;
    
    minCellVolume = p.Results.minCellVolume;
    cellDiameterRange = p.Results.cellDiameterRange;
    
    flagParallelize = p.Results.flagParallelize;
    flagDebugMode = p.Results.flagDebugMode;
    
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    totalSegTimer = tic;

    %% standardize the image
    imAdjusted = mat2gray( imInput ) * 4096;    

    %% apply median filter to remove any spotty noise
    fprintf( '\n\n>> Appling median filter to remove noise ...\n\n' );
    
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );    
    
    %% threshold the image using some algorithm
    fprintf( '\n\n>> Running a thresholding algorithm ...\n\n' );    
    
    localWindowRadius = round(30 ./ spacing(1));
    localWindowPace = round(localWindowRadius / 3);
    minLocalGlobalThresholdRatio = 0.6;
    
    switch thresholdingAlgorithm
        
        case 'OtsuSliceBySliceLocal'
            
            imThresh = segmentCellForegroundUsingLocalOtsu( imAdjusted, localWindowRadius, ...
                                                            'localWindowPace', localWindowPace, ...
                                                            'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                            'flagParallelize', flagParallelize );
                                                               

        case 'MinErrorPoissonSliceBySliceLocal'
            
                                                        
            imThresh = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                                'model', 'poisson', ...  
                                                                'localWindowPace', localWindowPace, ...
                                                                'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                'flagParallelize', flagParallelize );

                                                            
        otherwise
            
            error( 'ERROR: invalid thresholding method' );
        
    end

    imThresh = imfill( imThresh );        
    
    if flagDebugMode
        imseriesmaskshow( imInput, {imThresh}, 'spacing', spacing );
        set( gcf, 'Name', 'Thresholding Result' );
    end
    
    imAdjusted( ~imThresh ) = min( imAdjusted(:) );   
    
    %% detect cell seed points using multiscale LoG   
    fprintf( '\n\n>> Detecting seed points ...\n\n' );
    
    switch seedPointDetectionAlgorithm

        case 'IntensityMaxima'
            
            [imCellSeedPoints, imResponse] = detect_cell_seeds_IntensityMaxima( imAdjusted, ...
                                                                               mean( cellDiameterRange ), ...
                                                                               'spacing', spacing, ...
                                                                               'debugMode', true );
                                                                           
        case 'MultiscaleLoG'            
            
            numLoGScales = min(10, cellDiameterRange(2)-cellDiameterRange(1)+1);
            [imCellSeedPoints, imResponse] = detectBlobsUsingMultiscaleLoG( imAdjusted, ...
                                                                            cellDiameterRange, ...
                                                                            'numLoGScales', numLoGScales, ...
                                                                            'spacing', spacing, ...
                                                                            'debugMode', true );

        case 'AdaptiveMultiscaleLoG'
            
            fprintf( '\n\n\tComputing Distance Map of Binary Foreground Mask ...\n\n' );

            imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap',{} , imThresh, [], [], spacing  );
            imDistMap( imDistMap > 0 ) = 0;
            imDistMap = -1 * imDistMap;

            imseriesmaskshow( imDistMap, imThresh );
            set( gcf, 'Name', 'Foreground Distance Map' );    
            
            numLoGScales = min(10, cellDiameterRange(2)-cellDiameterRange(1)+1);
            [imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                                    imDistMap, ... 
                                                                                    'blobDiameterRange', cellDiameterRange, ...
                                                                                    'numLoGScales', numLoGScales, ...
                                                                                    'spacing', spacing, ...
                                                                                    'debugMode', true );
                                                                           
    end
    
    imCellSeedPoints = double( imCellSeedPoints > 0 );
    imCellSeedPoints(~imThresh) = 0;    
    
    if flagDebugMode
        
        imseriesmaskshow( max(imResponse,[],3), imdilate( max(imCellSeedPoints,[],3), strel('disk', 3) ) );
        set( gcf, 'Name', 'Seed Points Overlayed on MIP Response' );

        imseriesmaskshow( imResponse, imdilate(imCellSeedPoints, ones(3,3,3)) );
        set( gcf, 'Name', 'Seed Points Overlayed on Response' );

    end
    
    %% segment cells using random walks
    
    %% postprocess segmentation result    
    
    
end
    