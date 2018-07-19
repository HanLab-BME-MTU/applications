function [ imResized ] = imresizetotarget( imInput, targetSize, varargin )

    p = inputParser;
    p.addRequired( 'imInput', @(x) (max(size(x)) ~= numel(x) && ndims(x) >= 2) );
    p.parse( imInput );    
    p.addRequired( 'targetSize', @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.addParamValue( 'interpolationMethod', 'linear', @(x) ( ismember(x, {'cubic', 'linear', 'nearest'}) ) );
    p.parse( imInput, targetSize, varargin{:} );
    
    interpolationMethod = p.Results.interpolationMethod;
    
    % setup scaling transform
    resizeFactor = targetSize ./ size(imInput);
    scalingTransformMatrix = diag([resizeFactor,1]);
    scalingTransform = maketform( 'Affine' , scalingTransformMatrix );

    % image sampler
    linearImageResampler = makeresampler( interpolationMethod , 'replicate' );
    
    % apply transform
    imResized = tformarray( imInput, scalingTransform, linearImageResampler, ...
                            1:ndims(imInput), 1:ndims(imInput), ...
                            targetSize, [] , [] );
    
end