function BW = shadowSeg(I)
% This function segments the image by pulling out all of the darkest
% pixels
%
%   I = uint16 intensity matrix
%   BW = binary segmentation
  
    D = im2double(I);
    % First threshold using Rosin
    [N, X] = optimalHistogram(D(:));
    Ni = flipud(N); % Invert the histogram
    Xi = flipud(1-X);
    [cutoffIndex, cutoffValue] = cutFirstHistMode(Ni, Xi, 0);
    
    % Second threshold by finding biggest decrease in bin membership
    [dummy, mindiffIdx] = min(diff(Ni));
    t1 = Xi(mindiffIdx);
    t2 = cutoffValue;
    BW = hysteresisThreshold(1-D, t1, t2);
end