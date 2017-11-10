function [backMask,backMu,backSig] = estimateBackgroundArea(ims,varargin)
%ESTIMATEBACKGROUNDAREA estimates background mask, std and mean from input image
%
% [backMask,backMu,backSig] = estimateBackgroundArea(ims)
% [backMask,backMu,backSig] = estimateBackgroundArea(ims,ParamName',paramValue,...)
%
%  The input image can be a single 2D matrix or a 3D matrix of images.
%
%

%Hunter Elliott 10/2012


ip = inputParser;
ip.addParamValue('PostProcess',false,@(x)(islogical(x) && numel(x) == 1));%Mask post-processing
ip.addParamValue('nSTD',2,@(x)(numel(x) == 1));% # of STD above estimated background to threshold mask at
ip.addParamValue('TimeFiltSigma',3,@(x)(numel(x)==1 && x >=0)); %Sigma for gradient filtering in 3D if multiple images input
ip.addParamValue('ShowPlots',false,@(x)(numel(x)==1 && islogical(x))); %Display histogram fitting plots
%Post-processing parameters
ip.addParamValue('CloseRad',3,@(x)(numel(x)==1 && x >=0)); %Foreground closure radius
ip.addParamValue('AreaCutoff',[],@(x)(numel(x)==1 && x >=0)); %Foreground area/volume opening size

ip.parse(varargin{:});
p = ip.Results;


tSig = p.TimeFiltSigma;
nMax = 1e5;
subSamp = max(round(numel(ims)/nMax),1);
nIm = size(ims,3);

if isempty(p.AreaCutoff)
    %Adapt area to maximum-adjacency neighborhood given dimension
    if nIm > 1;
        p.AreaCutoff = 27;
    else
        p.AreaCutoff = 9;
    end
end

%Estimate background statistics
%[backMu,backSig] = estimateBackgroundStat(combIm(:)+eps)
[backMu,backSig] = fitGaussianModeToPDF(double(ims(1:subSamp:end)),'Display',p.ShowPlots);

%Get spatial and temporal gradient via filtering and add back to image.
%This helps in very-low SNR movies with motion. In other cases it
%contributes little to the final mask.
if nIm >= tSig && tSig > 0
    [dX,dY,dZ] = gradientFilterGauss3D(ims,tSig);
    tVar = sqrt(dX .^2 + dY .^2 + dZ .^2);
    combIm = double(ims) + (tVar*tSig);
else
    combIm = double(ims);
end

backMask = combIm < (backMu+backSig*p.nSTD);

%Do post-processing if requested.
if p.PostProcess    
    backMask = ~bwareaopen(~backMask,p.AreaCutoff);
    if nIm > 1
        backMask = ~imclose(~backMask,binarySphere(p.CloseRad));
    else
        backMask = ~imclose(~backMask,strel('disk',p.CloseRad));
    end
end
