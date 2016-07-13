function [mask,tracksOut] = trackPartitionInit(tracksIn,movieInfo,MD,threshold,minSize,upscale)
% Create upscaled mask and track struct
fprintf('Generating mask, please wait... \n')
if threshold == 0
    mask = constantSizeMasking(movieInfo,MD,minSize,upscale);
else
    mask = gaussianMasking(movieInfo,MD,threshold,minSize,upscale);
end

% Multiply upscale factor to every compound track 
newCoord = arrayfun(@(x) x.tracksCoordAmpCG*upscale,tracksIn,'UniformOutput',false);
% newCoordStruct = cell2struct(newCoord','tracksCoordAmpCG',1);

tracksOut = tracksIn;
[tracksOut.tracksCoordAmpCG] = newCoord{:};
% tracksOut = setfield(tracksOut,{1:size(tracksOut,1)},'tracksCoordAmpCG',{1},newCoordStruct);

% tracksOut(:).tracksFeatIndxCG = tracksIn.tracksFeatIndxCG;
% tracksOut(:).seqOfEvents = tracksIn.seqOfEvents;

% Make the output track struct have the same field name order as the input
tracksOut = orderfields(tracksOut,tracksIn);
end