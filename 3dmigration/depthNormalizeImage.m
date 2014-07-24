function normIm = depthNormalizeImage(im,maskOrDx,nucMask)
%DEPTHNORMALIZEIMAGE normalizes the masked image intensities based on their depth within the mask
%
% normIm = depthNormalizeImage(im,mask)
% normIm = depthNormalizeImage(im,distX)
% normIm = depthNormalizeImage(...,nucMask)
%
%   Normalizes the masked intensities within the input image based on
%   distance from the mask edge, so that at a given depth within the mask
%   the average intensity is always 1. This can optionally exclude values
%   within the input nuclear mask.
%
% ADD MORE DOCUMENTATION YOU LAZY BASTARD!!
%
% Hunter Elliott
% 2/2012

if islogical(maskOrDx)
    distX = bwdist(~maskOrDx);
else
    distX = maskOrDx;
end

%Avoid rounding error.
im = double(im);

%Get all unique distance transform values. Should I convert this to a
%lookup table, since it will be mostly the same for a given maximum depth?
distVals = unique(distX);
distVals(distVals==0) = [];

normIm = zeros(size(im));

for j = 1:numel(distVals)
   
   currPix = distX==distVals(j);
    
   normIm(currPix) = im(currPix) ./ mean(im(currPix(:) & ~nucMask(:)));
    
end


    




