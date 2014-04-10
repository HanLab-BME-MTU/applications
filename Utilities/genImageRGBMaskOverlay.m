function [ imMaskOverlay ] = genImageRGBMaskOverlay( im, rgbMask, maskAlpha )

    if ~iscell( rgbMask )
        rgbMask = { rgbMask };
    end

    if ~exist( 'maskAlpha', 'var' )
        maskAlpha = 0.5 * ones(1, numel(rgbMask));
    end
    
    if isscalar(maskAlpha)
        maskAlpha = maskAlpha * ones(1, numel(rgbMask));
    end
    
    imMaskOverlay = repmat( mat2gray(im), [1,1,3] );
    
    for i = 1:numel(rgbMask)
        
        blnMask = repmat( max( rgbMask{i}, [], 3 ) > 0, [1, 1, 3] );
        imMaskOverlay(blnMask) = (1 - maskAlpha(i)) * imMaskOverlay(blnMask) + maskAlpha(i) * rgbMask{i}(blnMask);
        imMaskOverlay( imMaskOverlay > 1 ) = 1;
        
    end

end