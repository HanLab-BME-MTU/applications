function [ x] = getLargestCC( mask )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mask = logical(mask); 
CC = bwconncomp(mask);
csize = cellfun(@(x) length(x), CC.PixelIdxList); 
CC.PixelIdxList(csize<max(csize))=[]; 
CC.NumObjects  = CC.NumObjects -sum(csize<max(csize)); 
x = labelmatrix(CC); 
x(x>0) = 1; 

end

