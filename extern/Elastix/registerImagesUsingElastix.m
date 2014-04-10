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
