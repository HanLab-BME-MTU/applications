function [ movie ] = MakeRGBMaskOverlayMovie( imInput, imRGBMask, alphaValue )

    if ~exist( 'alphaValue', 'var' )
        alphaValue = 0.2;
    end
    
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjustedLog = ComputeImageLogTransform( AdjustImageIntensityRange( imInput, ImageIntensityRange ) );        
    displayrange = [min(imAdjustedLog(:)), max(imAdjustedLog(:))];
    for i = 1:size(imInput,3)
       imSlice = im2uint8( mat2gray(imAdjustedLog(:,:,i), displayrange) );
       rgbSlice = cat(3, imSlice, imSlice, imSlice );
       rgbMaskSlice = squeeze(imRGBMask(:,:,i,:));
       rgbSlice = uint8( double((1 - alphaValue) * rgbSlice) + 255 * alphaValue * rgbMaskSlice );
       movie(i) = im2frame(rgbSlice); 
    end

end