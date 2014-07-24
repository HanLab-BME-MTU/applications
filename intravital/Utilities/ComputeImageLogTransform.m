function [ imLog ] = ComputeImageLogTransform( im )

    ImageIntensityRange = ComputeImageDynamicRange( im, 99.0 );    
    imLog = 1 + AdjustImageIntensityRange( im, ImageIntensityRange );
    imLog = log( imLog );
    
end