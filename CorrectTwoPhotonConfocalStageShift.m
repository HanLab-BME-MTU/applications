function [mergedImageData, imRegValidMask, varargout] = CorrectTwoPhotonConfocalStageShift( twoPhotonImageData, twoPhotonPixelSpacing, confocalImageData, confocalPixelSpacing, varargin )

    p = inputParser;
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( varargin{:} );
    PARAMETERS = p.Results;
    
    if PARAMETERS.flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end

    channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
    
    % display result before registration
    if PARAMETERS.flagDebugMode
    
        imseriesshow_multichannel( {twoPhotonImageData{1,1}, confocalImageData{1,1:2}}, ...
                                   'displayColors', channelcolormap, ...
                                   'spacing', twoPhotonPixelSpacing );
                           
        set( gcf, 'Name', 'Overlay of Moving and Fixed Volume Before Registration' );
        
    end
    
    % correct-shift using registration
    [mfilepath, name, ext] = fileparts( mfilename('fullpath') );    
    elastixParameters{1} = readElastixParameters( fullfile( mfilepath, 'elastixParamFile_rigid3D.txt' ) );
    mergedImageData{1} = twoPhotonImageData{1,1};    
    transformParameters = cell(1,2);
    
    imTwoPhoton = PreprocessData( twoPhotonImageData{1,1} );
    
    regTimer = tic;
    
    if PARAMETERS.flagParallelize

        imRegValidMask = cell(1,2);        

        parfor i = 1:2

            curRegTimer = tic;
            
            % setup fixed and moving image data
            dataFixed = struct;
            dataFixed.im = confocalImageData{1,i};
            dataFixed.Spacing = confocalPixelSpacing;   

            dataMoving = struct;
            dataMoving.im = imTwoPhoton;
            dataMoving.Spacing = twoPhotonPixelSpacing;
            
            % compute fixed mask
            imAdjusted = PreprocessData( dataFixed.im );
            fixedMask = segmentCellForeground( imAdjusted );
            fixedMask = imdilate(fixedMask, streldisknd([1,1,1]));

            % perform registration
            curElastixParameters = elastixParameters;
            curElastixParameters{1}.DefaultPixelValue = min(dataMoving.im(:)) - 1;
            [dataMovingR, regTransform, itersInfo] = registerImagesUsingElastix(dataFixed, dataMoving, curElastixParameters, fixedMask);

            % apply transform
            regTransform{1}.DefaultPixelValue = min(dataFixed.im(:)) - 1;
            regTransform{1}.TransformParameters = -1 * regTransform{1}.TransformParameters;
            dataFixedR = applyElastixTransform( dataFixed.im, dataFixed.Spacing, regTransform);            
            
            % mark out region valid for analysis where all channels overlap after registration
            imCurRegValidMask = (dataFixedR > regTransform{1}.DefaultPixelValue);            
            dataFixedR( ~imCurRegValidMask ) = min(dataFixed.im(:));
            
            % store output
            mergedImageData{i+1} = dataFixedR;                        
            transformParameters{i} = regTransform{1}.TransformParameters;
            imRegValidMask{i} = imCurRegValidMask;
            
            fprintf( '\nRegistration for confocal channel-%d took %f seconds', i, toc(curRegTimer) );            
            fprintf( '\nTransform Parameters for confocal channel-%d: [ %s ] \n', i, sprintf( ' %f ', transformParameters{i} ) );
            
        end        
        
        imRegValidMask = and( imRegValidMask{:} );        
        
    else
        
        imRegValidMask = cell(1,2);
        
        for i = 1:2
        
            curRegTimer = tic;
            
            % setup fixed and moving image data
            dataFixed.im = confocalImageData{1,i};
            dataFixed.Spacing = confocalPixelSpacing;   

            dataMoving.im = imTwoPhoton;
            dataMoving.Spacing = twoPhotonPixelSpacing;

            % compute fixed mask
            imAdjusted = PreprocessData( dataFixed.im );            
            fixedMask = segmentCellForeground( imAdjusted );
            fixedMask = imdilate(fixedMask, streldisknd([1,1,1]));

            % perform registration
            curElastixParameters = elastixParameters;
            curElastixParameters{1}.DefaultPixelValue = min(dataMoving.im(:)) - 1;
            [dataMovingR, regTransform, itersInfo] = registerImagesUsingElastix(dataFixed, dataMoving, curElastixParameters, fixedMask, PARAMETERS.flagDebugMode);

            % apply transform
            regTransform{1}.DefaultPixelValue = min(dataFixed.im(:)) - 1;
            regTransform{1}.TransformParameters = -1 * regTransform{1}.TransformParameters;
            dataFixedR = applyElastixTransform( dataFixed.im, dataFixed.Spacing, regTransform, PARAMETERS.flagDebugMode);            
            
            % mark out region valid for analysis where all channels overlap after registration
            imCurRegValidMask = (dataFixedR > regTransform{1}.DefaultPixelValue);            
            dataFixedR( ~imCurRegValidMask ) = min(dataFixed.im(:));
            
            % store output
            mergedImageData{i+1} = dataFixedR;                        
            transformParameters{i} = regTransform{1}.TransformParameters;
            imRegValidMask{i} = imCurRegValidMask;
            
            fprintf( '\nRegistration for confocal channel-%d took %f seconds', i, toc(curRegTimer) );            
            fprintf( '\nTransform Parameters for confocal channel-%d: [ %s ] \n', i, sprintf( ' %f ', transformParameters{i} ) );
            
        end                
        
        imRegValidMask = and( imRegValidMask{:} );
        
    end
    
    totalRegistrationTime = toc(regTimer);
    
    fprintf( '\nInter-channel registation took a total time of %f seconds\n', totalRegistrationTime );
    
    % display result after registration    
    if PARAMETERS.flagDebugMode
        
        imseriesshow_multichannel( mergedImageData, 'displayColors', channelcolormap, 'spacing', twoPhotonPixelSpacing );
        set( gcf, 'Name', 'Overlay of Moving and Fixed Volume After Registration' );        

    end
    
    if PARAMETERS.flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
    if nargout > 2
        varargout{1} = transformParameters;
    end
    
end

function [dataMovingR, regTransform, itersInfo, varargout] = registerImagesUsingElastix(dataFixed, dataMoving, elastixParameters, fixedMask, flagDebugMode)

    if ~exist( 'flagDebugMode', 'var' )
        flagDebugMode = false;
    end
    if flagDebugMode
        [dataMovingR, regTransform, itersInfo] = elastixRegister(dataFixed, dataMoving, elastixParameters, 'fixedMask', fixedMask);                         
    else
        [cmdlOutput, dataMovingR, regTransform, itersInfo] = evalc( 'elastixRegister(dataFixed, dataMoving, elastixParameters, ''fixedMask'', fixedMask)' );                         
        if nargout > 3
           varargout{1} = cmdlOutput;
        end
    end    
    
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

function [ imPreprocessed ] = PreprocessData( imInput )

    % standardize the intensity range
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 98.0 );
    imPreprocessed = mat2gray(imInput, ImageIntensityRange) * 4096;
    
%     % denoise it
%     imPreprocessed = ordfilt3(imPreprocessed, 'med', 3);    
    
end

function [imForegroundMask] = segmentCellForeground( imInput )

    %[level, imForegroundMask] = thresholdOtsu( ComputeImageLogTransform(imInput) );
    [level, imForegroundMask] = thresholdOtsu( imInput );
    
end
