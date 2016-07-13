% Create upscaled mask and track struct
function [mask,tracksOut] = trackPartitionInit(tracksIn,movieInfo,MD,threshold,minSize,upscale)

mask = gaussianMaskingInner(movieInfo,MD,threshold,minSize,upscale);
tracksOut = tracksIn;
tracksOut.tracksCoordAmpCG = arrayfun(@(x) x.tracksCoordAmpCG*upscale,tracksIn,'UniformOutput',false);


end