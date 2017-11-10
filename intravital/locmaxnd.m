function [lm] = locmax3d(img, windowRadius)

    ip = inputParser;
    ip.CaseSensitive = false;
    ip.addRequired('img', @(x) isnumeric(x));
    ip.addRequired('windowRadius', @(x) isnumeric(x) &&  all(x >= 1) && (isscalar(x) || numel(x) == ndims(img)));
    ip.parse(img, windowRadius);

    if any(windowRadius - floor(windowRadius) > 0)
        error('windowRadius must be integer');
    end
    
    if numel(windowRadius) == 1
        windowRadius = windowRadius * ones(1, ndims(img));
    end
    
    mask = true(2 * windowRadius + 1);
    lm = imdilate(img, mask);
    lm(lm~=img) = 0;
    
end