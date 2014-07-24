function [ intensityRange ] = ComputeImageIntensityRange( im )

    intensityRange = [min(im(:)) max(im(:))];
    
end