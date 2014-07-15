function [imMacrophageSeg] = segmentMacrophages(imMacrophage, varargin)

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('imInput', @(x) (ismember( ndims(x), [2,3] )));
    p.addParamValue('spacing', ones( 1, ndims(imMacrophage) ), @(x) (isnumeric(x) && numel(x) == ndims(imMacrophage)));
    p.addParamValue('maxObjectRadius', 20, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('minSignalToBackgroundRatio', 2.0, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('minObjectRadius', 2.0, @(x) (isscalar(x) && isnumeric(x)));
    p.parse(imMacrophage, varargin{:});
    
    PARAMETERS = p.Results;
    
    imMacrophageAdjusted = matitk('FMEDIAN', round(2 * min(PARAMETERS.spacing) ./ PARAMETERS.spacing), imMacrophage);
    
    imMacrophageSeg = thresholdSBR(imMacrophageAdjusted, ...
                                   PARAMETERS.maxObjectRadius, PARAMETERS.minSignalToBackgroundRatio, ...
                                   'spacing', PARAMETERS.spacing, ...
                                   'minObjectRadius', PARAMETERS.minObjectRadius, ...
                                   'kernelDimensions', 2, ...
                                   'downsamplingFactor', 0.5);

end