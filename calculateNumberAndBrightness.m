function [number,brightness] = calculateNumberAndBrightness(imageStack,background,backgroundVariance,varargin)

%calculateNumberAndBrightness takes an image stack and calculates the
%number and brightness of flourescent particles in each pixel

%estimate background and backgroundVariance from dark image
%as far as I can tell these are the average intensity and variance of a
%dark image....the average variance can also be calculated from the
%intercept of a variance vs intencity plot of a camera illuminated using
%differnt intensities (see most of Digman and Gratton's work)
%determine average and variance for each pixel

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageStack', @isnumeric);
ip.addRequired('background', @isnumeric);
ip.addRequired('backgroundVariance', @isnumeric);
ip.addParamValue('segmentLength', [], @isscalar);
ip.addParamValue('segmentSpacing', [], @isscalar);
ip.parse(imageStack,background,backgroundVariance,varargin);

segmentLength = ip.Results.segmentLength;
if isempty(segmentLength)
    segmentLength = size(imageStack,3);
end

number = nan(size(imageStack,1),size(imageStack,2),floor(size(imageStack,3)/segmentLength));
brightness = number;
for isegment = 1:floor(size(imageStack,3)/segmentLength)
%calculate number
number(:,:,isegment) = (nanmean(imageStack,3) - background).^2./(nanstd(imageStack,0,3).^2-backgroundVariance.^2);
%calculate brightness
brightness(:,:,isegment)  = (nanstd(imageStack,0,3).^2-backgroundVariance.^2)./(nanmean(imageStack,3) - background);
end
