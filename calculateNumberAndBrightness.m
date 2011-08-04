function [number,brightness] = calculateNumberAndBrightness(imageStack,background,backgroundVariance)

%calculateNumberAndBrightness takes an image stack and calculates the
%number and brightness of flourescent particles in each pixel

%estimate background and backgroundVariance from dark image
%as far as I can tell these are the average intensity and variance of a
%dark image....the average variance can also be calculated from the
%intercept of a variance vs intencity plot of a camera illuminated using
%differnt intensities (see most of Digman and Gratton's work)
%determine average and variance for each pixel


%calculate number
number = (nanmean(imageStack,3) - background).^2./(nanstd(imageStack,0,3).^2-backgroundVariance.^2);
%calculate brightness
brightness = (nanstd(imageStack,0,3).^2-backgroundVariance.^2)./(nanmean(imageStack,3) - background);
