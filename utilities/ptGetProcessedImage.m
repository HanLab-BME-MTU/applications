function [outImage, background] = ptGetProcessedImage (varargin)
% ptGetProcessedImage uses a number of image filtering functions to improve the
% image quality of the phase-contrast cell images so that an optimal segmentation 
% can take place. In particular the functions imSubtractBackground, imtophat,
% imbothat, imsubtract and imadd are used.
%
% SYNOPSIS       outImage = ptGetProcessedImage (varargin)
%
% INPUT          inputImage: the image that has to be processed
%                diskSize  : the size of the disk used in top- and bottom-hat
%
% OUTPUT         outImage  : the processed image which can be used for segmentation
%                background: the subtracted background
%
% DEPENDENCIES   ptGetProcessedImage uses { imSubtractBackground
%                                           imtophat
%                                           imbothat
%                                           imsubtract
%                                           imadd
%                                           strel }
%                                  
%                ptGetProcessedImage is used by { ptTrackCells }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Apr 04          Initial release

% Test the number of input variables
if nargin < 2
   error('The input image and disk size has to be provided. See help ptGetProcessedImage.');
end

% Get the input variables
inputImage = varargin{1};
diskSize = varargin{2};

% Generate the disk structuring element
se = strel ('disk', diskSize);

% Subtract the background from the input image
[backImage, background] = imSubtractBackground (inputImage);

% Get the tophat filtered image
tophatImage = imtophat (backImage, se);

% Get the bottomhat filtered image
bottomhatImage = imbothat (backImage, se);

% The output image has the nuclei and halos amplified to make segmenting easier
% later on in the process
outImage = imsubtract (imadd (tophatImage, backImage), bottomhatImage);

% Let's normalize the image back to [0..1] again
imageMinimum = min (min (outImage));
imageMaximum = max (max (outImage));
outImage = (outImage - imageMinimum) / (imageMaximum - imageMinimum);

% Do some cleaning up
clear se; clear backImage; 
clear tophatImage; clear bottomhatImage;
