function BW = shadowSeg(I)
% This function segments the image by pulling out all of the darkest
% pixels
%
%   I = uint16 intensity matrix
%   BW = binary segmentation
    
    % scale pixel intensities to avoid loss of precision
    Isc = imadjust(I, stretchlim(I,0), []);
    Dsc = im2double(Isc);
    
    % First threshold using Rosin
    [N, X] = histogram(Dsc(:));
    Ni = flipud(N); % Invert the histogram
    Xi = flipud(1-X);
    [cutoffIndex, cutoffValue] = cutFirstHistMode(Ni, Xi, 0);
    
    % Second threshold by finding biggest decrease in bin membership
    [dummy, mindiffIdx] = min(diff(Ni));
    t1 = Xi(mindiffIdx);
    t2 = cutoffValue;
    BW = hysteresisThreshold(1-Dsc, t1, t2);
end