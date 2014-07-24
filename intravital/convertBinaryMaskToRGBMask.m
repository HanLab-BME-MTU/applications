function [imMaskRGB] = convertBinaryMaskToRGBMask(imMask, maskColor)

    imMaskRGB = [];
    
    for i = 1:3
        imMaskRGB = cat(ndims(imMask)+1, imMaskRGB, imMask * maskColor(i));
    end

end