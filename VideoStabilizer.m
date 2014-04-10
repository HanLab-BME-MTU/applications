function [videoStabilized] = VideoStabilizer( videoInput, elastixParameterFile, varargin )

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'videoInput', @(x) (ismember(ndims(x), [3,4])) );
    p.addRequired( 'elastixParameterFile', @(x) (ischar(x) && exist(x, 'file')) );
    p.parse( videoInput, elastixParameterFile );
    
    numDims = ndims(videoInput) - 1; % last dimension assumed to be time
    imsize = size(videoInput);    
    numTimePoints = imsize(end);
    imsize(end) = [];
    
    p.addParamValue( 'spacing', ones(1,numDims), @(x) (isnumeric(x) && numel(x) == numDims) );
    p.addParamValue( 'wOldRef', 0.9, @(x) (isscalar(x) && isnumeric(x) && x > 0 && x < 1.0) );
    p.addParamValue( 'refTimePoint', round(0.5 * numTimePoints), @(x) (isscalar(x) && isnumeric(x) && x >= 1 && x <= numTimePoints) );
    p.parse( videoInput, elastixParameterFile, varargin{:} );
    
    PARAMETERS = p.Results;
    flagDebugMode = false;
    frameSlicer = repmat( {':'}, [1, ndims(videoInput)] );
    
    videoStabilized = videoInput;
    
    % read elastix parameter file
    elastixParameters{1} = readElastixParameters( elastixParameterFile );
    
    % stabilize forward in time starting from the reference time point
    frameSlicer{end} = PARAMETERS.refTimePoint;
    imRef = videoInput( frameSlicer{:} );
    
    h = waitbar(0, 'stabilizing video in forward direction' );
    
    for t = PARAMETERS.refTimePoint+1:numTimePoints
        
        frameSlicer{end} = t;
        
        imCurFrame = videoInput(frameSlicer{:});
        [imRegistered, regTransform, imValidROI] = registerImagesUsingElastix(imRef, imCurFrame, elastixParameters, PARAMETERS.spacing, flagDebugMode);
        imRegistered( ~imValidROI ) = imRef( ~imValidROI );
        
        imRef = PARAMETERS.wOldRef * imRef + (1 - PARAMETERS.wOldRef) * imRegistered;
        videoStabilized( frameSlicer{:} ) = imRegistered;
        
        waitbar( (t - PARAMETERS.refTimePoint) / (numTimePoints - PARAMETERS.refTimePoint), h );
        
    end
    
    close(h);
    
    % stabilize backward in time starting from the reference time point
    frameSlicer{end} = PARAMETERS.refTimePoint;
    imRef = videoInput( frameSlicer{:} );

    h = waitbar(0, 'stabilizing video in backward direction' );
    
    for t = PARAMETERS.refTimePoint-1:-1:1

        frameSlicer{end} = t;

        imCurFrame = videoInput(frameSlicer{:});
        [imRegistered, regTransform, imValidROI] = registerImagesUsingElastix(imRef, imCurFrame, elastixParameters, PARAMETERS.spacing, flagDebugMode);
        imRegistered( ~imValidROI ) = imRef( ~imValidROI );
        
        imRef = PARAMETERS.wOldRef * imRef + (1 - PARAMETERS.wOldRef) * imRegistered;
        videoStabilized( frameSlicer{:} ) = imRegistered;
        
        waitbar( (PARAMETERS.refTimePoint - t) / (PARAMETERS.refTimePoint - 1), h );
        
    end

    close(h)
    
end

function [imRegistered, regTransform, imValidROI, itersInfo] = registerImagesUsingElastix(imFixed, imMoving, elastixParameters, spacing, flagDebugMode)

    if ~exist( 'flagDebugMode', 'var' )
        flagDebugMode = false;
    end
    
    % setup fixed and moving image data
    dataFixed.im = imFixed;
    dataFixed.Spacing = spacing;   

    dataMoving.im = imMoving;
    dataMoving.Spacing = spacing;

    curElastixParameters{1}.DefaultPixelValue = min(dataMoving.im(:)) - 1;
    
    if flagDebugMode
        [dataMovingR, regTransform, itersInfo] = elastixRegister(dataFixed, dataMoving, elastixParameters);                         
    else
        [cmdlOutput, dataMovingR, regTransform, itersInfo] = evalc( 'elastixRegister(dataFixed, dataMoving, elastixParameters)' );                         
        if nargout > 3
           varargout{1} = cmdlOutput;
        end
    end    
    
    imRegistered = dataMovingR.im;
    imValidROI = (imRegistered > curElastixParameters{1}.DefaultPixelValue);
    
end

function [imRegistered, varargout] = applyElastixTransform(imInput, spacing, regTransform, flagDebugMode)

    if ~exist( 'flagDebugMode', 'var' )
        flagDebugMode = false;
    end

    if flagDebugMode
        imRegistered = elastixTransform( imInput, spacing, regTransform );            
    else
        [cmdlOutput, imRegistered] = evalc( 'elastixTransform( imInput, spacing, regTransform )' );                
        if nargout > 1
           varargout{1} = cmdlOutput;
        end
    end

end
