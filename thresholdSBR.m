function [imMask] = thresholdSBR(imInput, maxObjectRadius, sbrCutOff, varargin)

    p = inputParser;
    p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
    p.addRequired('maxObjectRadius', @(x) (isscalar(x) && isnumeric(x) && x > 0) );
    p.addRequired('sbrCutOff', @(x) (isscalar(x) && isnumeric(x) && x >= 1.0) );
    p.addParamValue('spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue('kernelDimensions', ndims(imInput), @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x <= ndims(imInput)));
    p.addParamValue('minObjectRadius', [], @(x) (isempty(x) || (isnumeric(x) && isscalar(x))));
    p.addParamValue('downsamplingFactor', 1.0, @(x) (isnumeric(x) && isscalar(x) && x > 0 && x <= 1));
    p.addParamValue('flagDebugMode', false, @(x) (islogical(x) && isscalar(x)));
    p.parse(imInput, maxObjectRadius, sbrCutOff, varargin{:});
    
    PARAMETERS = p.Results;
    
    % estimate local background using morphological opening
    if PARAMETERS.downsamplingFactor < 1
        imResized = imresizend(imInput, [PARAMETERS.downsamplingFactor * ones(1, PARAMETERS.kernelDimensions), ones(1, ndims(imInput) - PARAMETERS.kernelDimensions)] );
        %krnlMax = streldisknd( round(maxObjectRadius * PARAMETERS.downsamplingFactor ./ PARAMETERS.spacing(1:PARAMETERS.kernelDimensions)) );
        krnlMax = ones( 2 * round(maxObjectRadius * PARAMETERS.downsamplingFactor ./ PARAMETERS.spacing(1:PARAMETERS.kernelDimensions)) + 1 ); % flat kernel is much more faster
        imLocalBackground = imopen(imResized, krnlMax);
        imLocalBackground = imresizetotarget(imLocalBackground, size(imInput));
    else
        % krnlMax = streldisknd( round(maxObjectRadius ./ PARAMETERS.spacing(1:PARAMETERS.kernelDimensions)) );
        krnlMax = ones( 2 * round(maxObjectRadius ./ PARAMETERS.spacing(1:PARAMETERS.kernelDimensions)) + 1 ); % flat kernel is much more faster
        imLocalBackground = imopen(imInput, krnlMax);
    end
    
    % compute signal to background ration
    imSignalToBackgroundRatio = imInput ./ (eps + imLocalBackground);
    
    % apply cutoff
    imMask = imSignalToBackgroundRatio > PARAMETERS.sbrCutOff;

    % post-processing (optional)
    if ~isempty(PARAMETERS.minObjectRadius)
        krnlMin = streldisknd( round(PARAMETERS.minObjectRadius ./ PARAMETERS.spacing) );
        imMask = imopen(imMask, krnlMin);
    end
    
    % display stuff    
    if PARAMETERS.flagDebugMode
       imseriesmaskshow(imSignalToBackgroundRatio, imMask);
       set(gcf, 'Name', 'Foreground mask overlayed with signal to background ratio image');
    end
end


