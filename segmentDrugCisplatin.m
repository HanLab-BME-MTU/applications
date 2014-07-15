function [imDrugSeg] = segmentDrugCisplatin(imDrug, varargin)

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('imInput', @(x) (ismember( ndims(x), [2,3] )));
    p.addParamValue('spacing', ones( 1, ndims(imDrug) ), @(x) (isnumeric(x) && numel(x) == ndims(imDrug)));
    p.addParamValue('maxObjectRadius', 20, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('minSignalToBackgroundRatio', 2.0, @(x) (isscalar(x) && isnumeric(x)));
    p.addParamValue('minObjectRadius', [], @(x) (isscalar(x) && isnumeric(x)));
    p.parse(imDrug, varargin{:});
    
    PARAMETERS = p.Results;
    
    imDrugAdjusted = matitk('FMEDIAN', round(2 * min(PARAMETERS.spacing) ./ PARAMETERS.spacing), double(imDrug));
    
    imDrugSeg = thresholdSBR(imDrugAdjusted, ...
                             PARAMETERS.maxObjectRadius, PARAMETERS.minSignalToBackgroundRatio, ...
                             'spacing', PARAMETERS.spacing, ...
                             'minObjectRadius', PARAMETERS.minObjectRadius, ...
                             'kernelDimensions', 2, ...
                             'downsamplingFactor', 0.5);
    
end