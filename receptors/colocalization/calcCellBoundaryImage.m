function [threshold] = calcCellBoundaryImage(image,fixFrameUp, fixFrameDown)
%CALCCELLBOUNDARY applies mutli-otsu threshold to determine cell boundary
% Mask is best suited for images in which cell is large relative to image
% dimensions
% Synopsis: [maskList,mask] = calcCellBoundary(image)
% Input:
%       image - image that will be thresholded
%
%       fixFrameUp - since the function separates an image into multiple
%       tiers of intensity to find a threshold and the edge of the 
%       foreground maybe poorly defined, it may occur that the
%       threshold chosen is off by a tier. Enter any value to increase the
%       threshold by one tier.
%
%       fixFrameDown - similar to fixFrameUp, enter any value to decrease 
%       the chosen threshold by one tier
%
% Output:
%   threshold - value chosen to separate foreground and background of image 

%% Masking Process

 
     % Compute the thresholds
    Nvals = 1:20;
    metric = zeros(length(Nvals),1);
    for i = 1:length(Nvals)
        [~, metric(i)] = multithresh(image, Nvals(i) );
    end
 
     
    %Apply multi-Otsu threshold on image
    thresh = multithresh(image,Nvals(find(metric == (max(metric)))));
    
    %Attempt to find largest gap in thresholds
    diffThresh = zeros(length(thresh)-1,1);
    for i = 1:(length(thresh)-1)
       diffThresh(i) = thresh(i+1)/thresh(i);
    end
    threshLevel =1 + (find(diffThresh == max(diffThresh)));
    
    %Change threshold on tier higher or lower if selected
     if isempty(fixFrameUp) && isempty(fixFrameDown)
         threshold = thresh(threshLevel);
     elseif isempty(fixFrameDown)
         threshold = thresh(threshLevel+1);
     else
         threshold =  thresh(threshLevel-1);
     end
    
end