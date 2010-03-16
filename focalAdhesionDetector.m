function [Fit, Im] = focalAdhesionDetector(I, BW, sigma)

% The main goal of this function is to detect focal adhesions (FAs). A FA
% is assumed to be a diffraction-limited rod-like segment defined by:
% A: average amplitude along the segment
% B: average baseline intensity
% (x0,y0): the location of the center
% l: the length
% theta: the orientation
%
% We try to fit N segments to the image using a rod-like segment model
% (see dlSegment2D.m). The initial number of segment is defined by the
% number of connected components (CCs) provided by a initial coarse
% segmentation BW.
%
% STEP 1: Data fitting
%
% 1.1 Initial estimation
% 
% To get an initial candidate for the segment orientation, we filter the
% image using an anisotropic filter (streerableFiltering.m) and get the
% list of pixels with strictly positive value (after non-maximal
% suppression) and average their orientation locally for each connected
% component.
%
% The latter filtering step allows-us to remove outlier CCs when
% - there is no local maxima (or <= 2)
% or
% - the orientations are not normally distributed (kstest)
% or
% - the average response on NMS is too small
% According to these rules, we defined a boolean vector validCC which tells
% us which CC can be analyzed. When nnz(validCC) == 0, STEP 1 stops.
%
% The initial position (x0,y0) is defined using a Radon-like approach:
% Every pixel belonging to a CC are projected along the line oriented along
% theta0 and passing through the center of mass of CC. (x0,y0) is equal to
% the center of mass shifted by the mean value of the 1-dimensional
% projection of pixels.
%
% 1.2 Data Fitting
%
% We fit N rod-like segments to the image, one per CC, using the estimated
% initial model parameters (lsqnonlin).
%
% 1.3 Loop
%
% We iterate STEP 1 until non local maxima of the filtered image belongs to
% the CCs. We generate an image model from 1.2 and substrat it to the
% current image. I = I - Im;
%
% STEP 2: Merging (TO BE DEFINED)

I = double(I);

[R,T,NMS] = steerableFiltering(I, 2, sigma * sqrt(2)); % use M = 4
ind = find(NMS .* BW);
[y x] = ind2sub(size(I), ind);

edgesBW = double(edge(BW));
Is = (I - min(I(:))) / (max(I(:)) - min(I(:)));
Is(edgesBW == 1) = 0;
IC = repmat(Is,[1 1 3]);
IC(:,:,1) = IC(:,:,1) + edgesBW;

imshow(IC); hold on;
quiver(x, y, cos(T(ind)), sin(T(ind)), 0); hold off;