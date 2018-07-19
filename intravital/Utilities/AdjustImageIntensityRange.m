function [ imCorrected ] = AdjustImageIntensityRange( im, ImageIntensityRange )
    
    imCorrected = mat2gray( im, ImageIntensityRange ) * range(ImageIntensityRange);
    
end