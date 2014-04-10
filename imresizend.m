function [ imResized ] = imresizend( imInput, resizeFactor, varargin )

    p = inputParser;
    p.addRequired( 'imInput', @(x) (max(size(x)) ~= numel(x) && ndims(x) >= 2) );
    p.parse( imInput );    
    p.addRequired( 'resizeFactor', @(x) (isscalar(x) || (isnumeric(x) && numel(x) == ndims(imInput))) );
    p.addParamValue( 'interpolationMethod', 'linear', @(x) ( ismember(x, {'cubic', 'linear', 'nearest'}) ) );
    p.parse( imInput, resizeFactor, varargin{:} );
    
    interpolationMethod = p.Results.interpolationMethod;
    
    if isscalar( resizeFactor )
        resizeFactor = resizeFactor * ones(1,  ndims(imInput));
    end
    
    scalingTransformMatrix = diag([resizeFactor,1]);
    scalingTransform = maketform( 'Affine' , scalingTransformMatrix );
    linearImageResampler = makeresampler( interpolationMethod , 'replicate' );
    imResized = tformarray( imInput, scalingTransform, linearImageResampler, ...
                            1:ndims(imInput), 1:ndims(imInput), ...
                            round( size(imInput) .* resizeFactor ) , [] , [] );
    
end