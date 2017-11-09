function [BW2, td] = clusterSeg(BW1)
% This function finds the largest cluster of pixels from a provided binary
% image by thresholding on the distance transform
%
%   BW2 = binary segmentation
%   rad = radius of an individual pixel

    R = bwdist(BW1); % look for ridges in distance between mask points
    td = thresholdRosin(R);  
    BW2 = R > td;
end