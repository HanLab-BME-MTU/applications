function [number,brightness] = calculateNumberAndBrightness(imageStack,background,backgroundStd,varargin)

%calculateNumberAndBrightness takes an image stack and calculates the
%number and brightness of flourescent particles in each pixel

%estimate background and backgroundStd from dark image
%as far as I can tell these are the average intensity and variance of a
%dark image....the average variance can also be calculated from the
%intercept of a variance vs intencity plot of a camera illuminated using
%differnt intensities (see most of Digman and Gratton's work)
%determine average and variance for each pixel

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageStack', @isnumeric);
ip.addRequired('background', @isnumeric);
ip.addRequired('backgroundStd', @isnumeric);
ip.addParamValue('segmentLength', 0, @isscalar);
ip.addParamValue('segmentSpacing', 1, @isscalar);
ip.addParamValue('S', 0, @isscalar);
ip.parse(imageStack,background,backgroundStd,varargin{:});
segmentSpacing = ip.Results.segmentSpacing;
S = ip.Results.S;

segmentLength = ip.Results.segmentLength;
if segmentLength == 0
    segmentLength = size(imageStack,3);
end
%make vector of segment starts
%the last segent must fit fully within the available number of frames
segmentStarts = 1:segmentLength+segmentSpacing-1:size(imageStack,3);
segmentStarts = segmentStarts(find(segmentStarts <= size(imageStack,3)-segmentLength+1));
display(num2str(segmentStarts))

number = nan(size(imageStack,1),size(imageStack,2),length(segmentStarts));
brightness = number;
for isegment = 1:length(segmentStarts)
    segmentImages = imageStack(:,:,segmentStarts(isegment):segmentStarts(isegment)+segmentLength-1);
    %calculate brightness
    brightness(:,:,isegment)  = (nanstd(segmentImages,0,3).^2-backgroundStd.^2 - ...
        S*(nanmean(segmentImages,3) - background))./(nanmean(segmentImages,3) - background);
    
    %calculate number
    number(:,:,isegment) = (nanmean(segmentImages,3) - background)./brightness(:,:,isegment);
    
%     brightness(:,:,isegment)  = (nanstd(segmentImages,0,3).^2-backgroundStd.^2)./...
%         (nanmean(segmentImages,3) - min(segmentImages,[],3)); % - background);

end