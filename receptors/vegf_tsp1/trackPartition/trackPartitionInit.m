function [mask,tracksOut] = trackPartitionInit(tracksIn,movieInfo,MD,threshold,minSize,upscale)
%TRACKPARTITIONINIT Prepares mask and track struct for use by
%trackPartition.
%   [mask,tracksOut] = trackPartitionInit(tracksIn,movieInfo,MD,threshold,minSize,upscale)
%
%   Inputs:
%       tracksIn:       track struct, such as one produced by
%                       TrackingProcess
%
%       movieInfo:      particle detection result
%
%       MD:             Movie Data with particles to be masked
%
%       threshold:      threshold for producing masks from the Gaussian
%                       detection information. Masks will have value 'true'
%                       where the fitted Gaussian exceeds the threshold.
%                       Set to 0 (default) to use constant size masks.
%
%       minSize:        minimum mask diameter (in nm); if threshold is set
%                       to 0 for constant size masking, this parameter 
%                       specifies the diameter of the masks; default is 100 nm
%
%       upscale:        factor by which to upscale the mask track coordinates
%                       from the size and resolution of the movie. This 
%                       allows a more circular mask when dealing with mask 
%                       diameters that are only a couple of pixels in the   
%                       original movie image size. (default 1)
%
%   Outputs:
%       mask:           binary mask with value 'true' at particle locations
%
%       tracksOut:      original track struct with every coordinate
%                       multiplied by upscale
%
%Kevin Nguyen, July 2016
fprintf('Generating mask, please wait... \n')
if threshold == 0
    mask = constantSizeMasking(movieInfo,MD,minSize,upscale);
else
    mask = gaussianMasking(movieInfo,MD,threshold,minSize,upscale);
end

tracksOut = trackScalar(tracksIn,upscale);
end