function [n, binMap] =speckleDetector(img,parameters)
%SPECKLEDETECTOR selects significant local intensity maxima based on a noise model
%
% SYNOPSIS [n, binMap] = speckleDetector(img,parameters)
%
% INPUT      img         : image
%            parameters  : [k sDN sP]
%                                k: parameter
%                                sDN : dark noise
%                                sP  : Poisson noise (dependent on I)
%
% OUTPUT     n           : number of speckles
%            binMap      : binarized map of the speckles

% Select local maxima as speckle candidates
lMax = locmax2d(img,[5 5]);

% Local minima will be used for Delaunay triangulation
lMin = locmin2d(img,[3 3]);   % Mask set to [3 3] to get more points


% Perform Delaunay triangulation
info=analyzeSpeckles(size(img),lMax,lMin);

% Loc max are selected based on a noise model % parameters=[0.4 2.73e-4 2e-4];
% but check first for pathalogical cases where info is empty
if isempty(info) 
   n = 0;
   binMap = zeros(size(img));
   return;
end;
lMax=validateSpeckles(img,lMax,info,parameters);

% Count speckles
binMap = lMax>0;
n = length(find(binMap));
