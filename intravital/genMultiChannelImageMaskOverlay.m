function [ imMaskOverlay ] = genMultiChannelImageMaskOverlay( im, channelDisplayColors, masks, maskColors, maskAlphas )

    imChannelOverlay = genMultiChannelOverlay(im, channelDisplayColors);

    imr = imChannelOverlay(:,:,1);
    img = imChannelOverlay(:,:,2);
    imb = imChannelOverlay(:,:,3);
    
    if ~iscell( masks )
        masks = { masks };
    end
       
    for i = 1:numel(masks)
        
        curMask = logical(masks{i});
        curMaskColor = maskColors(i,:);
        curMaskAlpha = maskAlphas(i);
                
        imr(curMask) = double( (1 - curMaskAlpha) * imr(curMask) + curMaskAlpha * curMaskColor(1) );
        img(curMask) = double( (1 - curMaskAlpha) * img(curMask) + curMaskAlpha * curMaskColor(2) );
        imb(curMask) = double( (1 - curMaskAlpha) * imb(curMask) + curMaskAlpha * curMaskColor(3) );
        
    end
    
    imMaskOverlay = squeeze(cat(ndims(im)+1, imr, img, imb));
    imMaskOverlay( imMaskOverlay > 1 ) = 1;
    
end