function [localMax,imgLM]=findLocalMaxima(img,sigma,varargin)
%FINDLOCALMAXIMA -> local maxima detection in image array
%   This function detects local maxima in an image using first a Laplacian
%   of Gaussian filter and then an amplitude criterion to determine regions
%   where local maxima are significant. It also excludes maxima that are 
%   closer to the image boundary than four times the expected width of the 
%   point pread function. 
%   The code is in part adjusted from 'pointSourceDetection'.
%   
%   Input:
%       img: input image
%       sigma: standard deviation of Gaussian PSF
%       {'alpha',alpha}: alpha value for significance test (default: 0.01)
%       {'pMask',pMask}: alpha value for generating mask (default: 0.05)
%         {'mask',mask}: cell mask (binary array)
%   Output:
%       localMaxima: structure with fields
%           .x: x coordinates of local maxima
%           .y: y coordinates of local maxima
%           .amp: amplitude of local maxima
%       imgLM: image array with local maxima positions
%
%   2011/07/27, US
%

ip=inputParser;
ip.CaseSensitive=false;
ip.addRequired('img',@isnumeric);
ip.addRequired('sigma',@isscalar);
ip.addParamValue('alpha',0.01,@isscalar);
ip.addParamValue('pMask',0.05,@isscalar);
ip.addParamValue('mask',true(size(img)),@islogical);

ip.parse(img,sigma,varargin{:});
alpha=ip.Results.alpha;
pMask=ip.Results.pMask;

cellMask=ip.Results.mask;


sizeX=size(img,1);
sizeY=size(img,2);

% Gaussian Kernel
w=ceil(4*ceil(sigma));
x=-w:w;
g=exp(-x.^2/(2*sigma^2));
u=ones(1,length(x));

% convolutions
imgXT=padarrayXT(img,[w w],'symmetric');
fg=conv2(g',g,imgXT,'valid');
fu=conv2(u',u,imgXT,'valid');
fu2=conv2(u',u,imgXT.^2,'valid');

% Laplace of Gaussian
gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - ...
    (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

% 2-D kernel
g = g'*g;
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% solution to linear system
% A_est: estimated amplitude, c_est: estimated background
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

J = [g(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
sigma_e2 = RSS/(n-3);

% standard error of amplitude A
sigma_A = sqrt(sigma_e2*C(1,1));

% standard deviation of residuals
sigma_res = sqrt((RSS - (A_est*gsum+n*c_est - fu)/n)/(n-1));

kLevel = norminv(1-alpha/2.0, 0, 1);

SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
T = (A_est - sigma_res*kLevel) ./ scomb;
pval = tcdf(real(T), df2);

% mask of admissible positions for local maxima
mask = pval > 1.0-pMask;

% all local max, keep also flat maxima
allMax = locmax2d(imgLoG, [3 3],1);

% local maxima above threshold in image domain
imgLM = allMax .* cellMask;
imgLM=imgLM.*mask;

if sum(imgLM(:)) ~= 0
    
    % -> set threshold in LoG domain
    logThreshold = min(imgLoG(imgLM~=0));
    logMask = imgLoG >= logThreshold;
    
    % combine masks
    mask = mask | logMask;
    
    % re-select local maxima
    imgLM=allMax.*cellMask;
    imgLM = imgLM .* mask;
    
    % find maxima, discard maxima too close to rim of image
    [lmx,lmy] = find(imgLM~=0);
    indX=lmx > w & lmx <= sizeX-w;
    indY=lmy > w & lmy <= sizeY-w;
    ind=indX & indY;
    lmx=lmx(ind);
    lmy=lmy(ind);
    
    nMax=numel(lmx);
    amp=NaN(nMax,1);
    da=NaN(nMax,1);
    bg=NaN(nMax,1);
    
    for iMax=1:nMax
        i=lmx(iMax);
        j=lmy(iMax);
        amp(iMax)=A_est(i,j);
        da(iMax)=sigma_A(i,j);
        bg(iMax)=c_est(i,j);
       
        %cropROI=img(i-w:i+w,j-w:j+w);
        %cropROI=cropROI(gMaskBin);
        %meanBG=mean(cropROI(:));
        %stdBG=std(cropROI(:));
        %bg(iMax)=meanBG;
        %pValue=1-normcdf(amp(iMax),meanBG,stdBG);
        %f pValue < 0.001
        %   keepMax(iMax)=true;
        %nd
    end
        
    localMax.x=lmx;
    localMax.y=lmy;
    localMax.amp=amp;
    localmax.da=da;
    localMax.bg=bg;
else
    localMax.x=[];
    localMax.y=[];
    localMax.amp=[];
    localMax.da=[];
    localMax.bg=[];
end