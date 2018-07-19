function [imThresh] = thresholdBlobness(imInput, blobDiameter, varargin)

    p = inputParser;
    p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
    p.addRequired( 'blobDiameter', @(x) isscalar(x) );

    p.addParamValue('choice_of_threshold', 'MinimumError', ...
                    @(x) ((ischar(x) && ismember(x, {'Otsu', 'Rosin', 'MinimumError'})) || ...
                          isa(x, 'function_handle')));
    p.addParamValue('spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue('flagDebugMode', false, @(x) (islogical(x) && isscalar(x)));
    p.parse(imInput, blobDiameter, varargin{:});
    
    PARAMETERS = p.Results;

    sigmaLoG = 0.5 * blobDiameter / sqrt(ndims(imInput));
    imLoG = -1 * filterLoGND(imInput, sigmaLoG, 'spacing', PARAMETERS.spacing, 'UseNormalizedDerivatives', true);
    
    if isa(PARAMETERS.choice_of_threshold, 'function_handle')
        level = PARAMETERS.choice_of_threshold(imLoG);                
    else        
        switch PARAMETERS.choice_of_threshold
            case 'Otsu'
                level = thresholdOtsu(imLoG);
            case 'Rosin'
                level = thresholdRosin(imLoG);
            case 'MinimumError'
                level = thresholdMinimumError(imLoG);
        end
    end    
    
    imThresh = imLoG > level;
    
end