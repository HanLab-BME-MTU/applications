function [I0,sDN,GaussRatio]=calculateStackNoiseParam(stack,bitDepth,sigma)
% calculateStackNoiseparam calculates three parameters for the noise model applied to speckle selection
%
% Run this funtion on a background region (non-speckled) cropped from an image stack 
% (use cropStack for this purpose). 
% calculateStackNoiseparam returns the mean backgroung intensity (I0), the mean standard deviation
% (sDN) and GaussRatio, which indicates how well the dark noise of the camera approximates 
% a normal distribution. More precisely, GaussRatio corresponds to the ratio:
%
%                                std(background image)
%                           -------------------------------,
%                            std(low-pass filtered image)
%
% where 'low-pass filtered image' is the result of the convolution of the image with a 
% Gaussian kernel with standard deviation = sigma.
%
% SYNOPSIS [I0,sDN,GaussRatio]=calculateStackNoiseparam(stack,bitDepth,sigma)
%
% INPUT    stack        : stack of images - should be a n x m ... x z matrix
%                         where z is the dimension of the stack
%          bitDepth     : bit depth of the camera for normalization
%          sigma        : sigma for the gaussian kernel
%
% OUTPUT   I0           : average intensity
%          sDN          : average standard deviation (sigmaDarkNoise)
%          GaussRatio   : ratio std(image)/std(filtered_image)
%
% Copied from fsmCalcNoiseParam

% Input check
ip =inputParser;
ip.addRequired('stack',@isnumeric);
ip.addRequired('bitDepth',@isscalar);
ip.addRequired('sigma',@isscalar);
ip.parse(stack,bitDepth,sigma);

% Normalization boundaries
mx=2^bitDepth-1;
mn=0;

% Compute Ascombes transform (SB - 4/8/10)
% DON'T DO THAT UNTIL IT IS THOROUGHLY CHECKED
%stack = 2 * sqrt(stack + 3/8);

% Normalize stack
stack=(double(stack)-mn)/(mx-mn);

% Mean
I0=mean(stack(:));

% Standard deviation

%(Aaron calculates the sDN by taking the std over time for each pixel, and
%then averaging over all pixels in the chosen part of the background. This
%does not work if you give the software just one image (not a movie), and is 
%probably quite inaccurate for movies that have a small number of time points. 
%So I changed it such that, if the number of images in a movie is less than 5, 
%the sDN is calculated as the std of all pixels in the chosen part at all
%time points. -KJ)

n=size(stack);
if n(end) < 5
    sDN = std(stack(:),1);
else
    % Standard deviation over time
    S=std(stack,1,3);

    % Get a mean value for the standard deviation over time
    sDN=mean(S(:));
end

% Calculate GaussRatio
%
% Calculate how many pixels have to be cropped from the borders
border=3*fix(sigma)+1;

% Initialize vector
nFrames=size(stack,3);
GaussRatios=zeros(1,nFrames);

% Calculate all ratios
for i=1:nFrames
    % Get current image
    rImg=stack(:,:,i);
    % Filter it with user-input sigma
    fImg=filterGauss2D(rImg, sigma);
    % Crop border
    rImg=rImg(border:end-border+1,border:end-border+1);
    fImg=fImg(border:end-border+1,border:end-border+1);
    % Calculate current ratio
    GaussRatios(i)=std(rImg(:))/std(fImg(:));
end

% Calcaulate GaussRatio as the mean of the vector GaussRatios
GaussRatio=mean(GaussRatios);