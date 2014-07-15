function [ imMaskOverlay ] = genMultiChannelImageRGBMaskOverlay( im, channelDisplayColors, rgbMask, maskAlpha )

    if ~iscell( rgbMask )
        rgbMask = { rgbMask };
    end

    if ~exist( 'maskAlpha', 'var' )
        maskAlpha = 0.5 * ones(1, numel(rgbMask));
    end
    
    if isscalar(maskAlpha)
        maskAlpha = maskAlpha * ones(1, numel(rgbMask));
    end
    
    imMaskOverlay = genMultiChannelOverlay(im, channelDisplayColors);
    
    for i = 1:numel(rgbMask)
        
        blnMask = repmat( max( rgbMask{i}, [], 3 ) > 0, [ones(1,ndims(im)), 3] );
        imMaskOverlay(blnMask) = (1 - maskAlpha(i)) * imMaskOverlay(blnMask) + maskAlpha(i) * rgbMask{i}(blnMask);
        imMaskOverlay( imMaskOverlay > 1 ) = 1;
        
    end

end