function [imFiltered] = medfilt3(imInput, windowRadius, varargin)
    
    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ndims(x) == 3));
    p.addRequired('windowRadius', @(x) (isnumeric(x) && (isscalar(x) || numel(x) == 3) && ~any(x - floor(x) > 0)) );
    p.addOptional( 'padval', 'replicate', @(x) (isscalar(x) || (ischar(x) && ismember( lower(x), {'zeros', 'symmetric', 'replicate', 'circular' } ))) );
    p.parse(imInput, windowRadius, varargin{:});
    
    padval = p.Results.padval;
    
    if isscalar(windowRadius)
        windowRadius = windowRadius + zeros(1,3);
    end
    
    % pad the input image with windowRadius 
    imInputPadded = padarray(imInput, windowRadius, padval); 
    
    % call mex code for median filtering
    imFiltered = medfilt3mex(imInputPadded, windowRadius);
    
    % unpad the result
    indCropBox = cell(1,3);
    for i = 1:3
        indCropBox{i} = windowRadius(i) + (1:size(imInput,i));
    end
    imFiltered = imFiltered( indCropBox{:} );

end