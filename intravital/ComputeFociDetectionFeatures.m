function [featureStruct] = ComputeFociDetectionFeatures(imInput, fociStats, fid, varargin)

    p = inputParser;
    p.addRequired( 'imInput', @(x) (ismember( ndims(x), [2,3])) ); 
    p.addRequired( 'fociStats', @(x) ( isstruct(x) && all(isfield(x, {'PixelLocationIndex', 'Radius'})) ));
    p.addRequired( 'fid', @(x) ( isnumeric(x) && isscalar(x) && ~(x - floor(x) > 0) && x <= numel(fociStats)));
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.parse(imInput, fociStats, fid, varargin{:});

    curFociStats = fociStats(fid);
    
    featureStruct.SignalToBackgroundRatio = curFociStats.SignalToBackgroundRatio;
    featureStruct.DoGResponse = curFociStats.DoGResponse;

    featureStruct.HessianEigenValues = curFociStats.HessianEigenValues;
    featureStruct.HessianEigenRatio21 = curFociStats.HessianEigenValues(2) / curFociStats.HessianEigenValues(1);
    featureStruct.HessianEigenRatio31 = curFociStats.HessianEigenValues(3) / curFociStats.HessianEigenValues(1);
    featureStruct.HessianDeterminant = prod(curFociStats.HessianEigenValues);
    
end