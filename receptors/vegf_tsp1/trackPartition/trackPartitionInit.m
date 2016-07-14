function [mask,tracksOut] = trackPartitionInit(tracksIn,movieInfo,MD,threshold,minSize,upscale)
% Create upscaled mask and track struct
fprintf('Generating mask, please wait... \n')
if threshold == 0
    mask = constantSizeMasking(movieInfo,MD,minSize,upscale);
else
    mask = gaussianMasking(movieInfo,MD,threshold,minSize,upscale);
end

tracksOut = trackScalar(tracksIn,upscale);
end