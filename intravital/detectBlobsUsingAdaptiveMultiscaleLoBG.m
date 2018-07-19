function [ imBlobSeedPoints, varargout ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% Detects blobs as local maxima of the Multiscale LoBG filtered image
%
% [ imBlobSeedPoints ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% [ imBlobSeedPoints, imMultiscaleLoGResponse ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% [ imBlobSeedPoints, imMultiscaleLoGResponse, imBlobSize ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% 
%  The maximum scale at each pixel is set using a distance map of the
%  binary mask of cell foreground obtained using a thresholding algorithm.
% 
%  This is claimed to be better than detecting blobs using multiscale LoG
%  because setting maximum sigma/scale for each pixel is based on its distance
%  to the nearest background pixel avoids overblurring.
% 
%  The input intensity image is first filtered with a Laplacian of Bi-Gaussian
%  (LoBG) Filter accross multiple scales/sigmas. The blobs are then detected as 
%  local maxima in the scale space.
% 
%  The maximum scale at each pixel is set using a distance map of the
%  binary foreground mask obtained using a thresholding algorithm.
%  
%   References:
% 
%   Xiao, C., M. Staring, et al. (2012). 
%   "A multiscale bi-Gaussian filter for adjacent curvilinear structures 
%   detection with application to vasculature images." 
%   IEEE Transactions on Image Processing, PP(99): 1-1.
% 
%   See: filterMultiscaleLoBGND.m, filterLoBGND.m
% 
%   Author: Deepak Roy Chittajallu
% 
% 
