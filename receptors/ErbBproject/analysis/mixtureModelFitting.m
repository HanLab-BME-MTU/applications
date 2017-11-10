function [ features ]=mixtureModelFitting(img,localMaxima,sigma,varargin)
%MIXTUREMODELFITTING Gaussian mixture model fitting to image data
%   Detailed explanation goes here
%   Input:
%       img: image array
%       localMaxima: structure with fields as output from 'findLocalMaxima'
%           .x: x coordinate of detected maxima
%           .y: y coordinate of detected maxima
%           .amp: amplitude of maxima
%           .da: standard error of amplitude 
%           .bg: estimated background around maxima
%       {'alphaA',alphaA}: alpha value for amplitude test (default: 0.05)
%       {'alphaD',alphaD}: alpha value for distance test (default: 0.05)
%       {'alphaF',alphaF}: alpha value for model selection (default: 0.05)
%       {'alphaR',alphaR}: alpha value for residual test (default: 0.05)
%
%   Output:
%       features: array containing information about PSF centroids
%           1st, 2nd column: x, y coordinate
%           3rd column: signal amplitude
%           4th column: signal width (PSF sigma)
%           5th to 8th column: standard deviation of above quantities
%           9th column: valency, i.e. size of cluster the signal lies in
%
%   2011/07/27 US           
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;
ip.addRequired('img',@isnumeric);
ip.addRequired('localMaxima',@isstruct);
ip.addRequired('sigma',@isscalar);
ip.addParamValue('alphaA',0.05,@isscalar);
ip.addParamValue('alphaD',0.05,@isscalar);
ip.addParamValue('alphaF',0.05,@isscalar);
ip.addParamValue('alphaR',0.05,@isscalar);

ip.addOptional('doMMF',false,@islogical);

ip.parse(img,localMaxima,sigma,varargin{:});
alphaA=ip.Results.alphaA;
alphaD=ip.Results.alphaD;
alphaF=ip.Results.alphaF;
alphaR=ip.Results.alphaR;

doMMF=ip.Results.doMMF;

Q=struct('alphaA',alphaA,'alphaD',alphaD','alphaF',alphaF,'alphaR',alphaR);

% clear return variable
features=[];

[sizeX,sizeY]=size(img);

s0=max(ceil(sigma),1.0);

% re-order localMaxima structure to be compatible with already existing
% functions like 'findOverlapPSFs2D'
nMax=size(localMaxima.x,1);
cands = repmat(...
    struct('status',1,'IBkg',[],'Lmax',[],'amp',[],'pValue',[]),nMax,1);
for kMax=1:nMax
    cands(kMax).IBkg=localMaxima.bg(kMax);
    cands(kMax).Lmax=[localMaxima.x(kMax) localMaxima.y(kMax)];
    cands(kMax).amp=localMaxima.amp(kMax);
    %cands(kMax).pValue=localMaxima.pVal(kMax);
end

% find signals that are less than a certain distance apart
% one should be able to give this distance as an argument to findOverlapPSFs

[clusters,~]=findOverlapPSFs2D(cands,sizeX,sizeY,sigma,0,img);

% how much signal clusters do exist?
nClust=size(clusters,2);

% sort cluster according to size
tmp=vertcat(clusters.numMaxima);
[tmp,index]=sort(tmp);
clusters=clusters(index);
% number of clusters of size 1
numSingles=sum(tmp == 1);

% options for fitting routine
%maxIter=500;
% epsAbs=1e-6;
% epsRel=1e-6;
% fitOptions=[maxIter epsAbs epsRel];
    
% fit isolated signals first
for kClust=1:numSingles
    x0=clusters(kClust).maximaPos(1);
    y0=clusters(kClust).maximaPos(2);
    A0=clusters(kClust).maximaAmp;
    
    clusterRegion=img(y0-4*s0:y0+4*s0,x0-4*s0:x0+4*s0);
    bg=min(clusterRegion(:));
    
%     [prmVec prmStd C res J]=...
%        fitGaussian2D(clusterRegion,[0 0 A0 s0 bg],'xyAsc',fitOptions);
%     
%     
%     %collect sub-pixel positions, signal amplitudes & widths
%     %and their corresponding errors in array 'features'
%     tmp=zeros(1,9);
%     if res.pval > alphaR
%        tmp(1)=x0;
%        tmp(2)=y0;
%        tmp=tmp+[prmVec(1:4) prmStd(1:4) 1];
%     end
%     features=[features; tmp];
    
    prmVec=[0 0 A0 sigma bg];
    [prmVec prmStd C res J status]=...
        doRepetitiveMMF(clusterRegion,prmVec,Q,'doMMF',doMMF);
    
    if status
        prmVec(:,1)=prmVec(:,1)+x0;
        prmVec(:,2)=prmVec(:,2)+y0;
        features=[features;[prmVec prmStd]];
    end
%         s=prmVec(end-1);
%         ds=prmStd(end-1);
%         nSig=(numel(prmVec)-2)/3;
%         for k=0:nSig-1% numSignal-1
%             x=prmVec(k*3+1)+x0;
%             y=prmVec(k*3+2)+y0;
%             a=prmVec(k*3+3);
%             
%             dx=prmStd(k*3+1);
%             dy=prmStd(k*3+2);
%             da=prmStd(k*3+3);
%             
%             if dx > 1.0 || dy > 1.0
%                 continue;
%             end
%             
%             tmp=[x y a s dx dy da ds nSig prmVec(end)];
%             
%             features=[features;tmp];
%        end
%    end
    
end

% continue with higher-order signal clusters using mixture-model fitting
for kClust=numSingles+1:nClust
    % axis convention:
    % x-axis: from left to right
    % y-axis: from top to bottom
    % left part of image corresponds to x = -inf, right to x = +inf
    % top part of image corresponds to y = -inf, bottom part to y = +inf
    
    % number, positions and amplitudes of signals
    numSignal=clusters(kClust).numMaxima;
    posSignal=clusters(kClust).maximaPos(:,1:2);
    ampSignal=clusters(kClust).maximaAmp;
    
    % min/max coordinates to crop region of cluster
    % column: x position, rows: y position
    xCropMin=min(posSignal(:,2))-3*s0;
    xCropMax=max(posSignal(:,2))+3*s0;
    yCropMin=min(posSignal(:,1))-3*s0;
    yCropMax=max(posSignal(:,1))+3*s0;
    
    clusterRegion=img(xCropMin:xCropMax,yCropMin:yCropMax);
    
    [data,oX,oY]=makeSquareImage(clusterRegion);
    
    allPrm=NaN(1,numSignal*3);
    bg=min(clusterRegion(:));
    sigmaAndBG=[sigma bg];
    ss=floor(size(data,1)/2)+1;
    
    for k=0:numSignal-1
        allPrm(1,k*3+1)=posSignal(k+1,1)-yCropMin+oX+1-ss;
        allPrm(1,k*3+2)=posSignal(k+1,2)-xCropMin+oY+1-ss;
        allPrm(1,k*3+3)=ampSignal(k+1);
    end
    
    prmVec=[allPrm sigmaAndBG];
    [prmVec prmStd C res J status]=doRepetitiveMMF(data,prmVec,'doMMF',doMMF);
        
    % transform fitted positions back to real image coordinates
    % save sub-pixel position amplitude and width in array 'features'
    
    %if res.pval > alphaR
    if status
        prmVec(:,1)=prmVec(:,1)+yCropMin-oX-1+ss;
        prmVec(:,2)=prmVec(:,2)+xCropMin-oY-1+ss;
        features=[features;[prmVec prmStd]];
    end
%     if status
%         s=prmVec(end-1);
%         ds=prmStd(end);
%         nSig=(numel(prmVec)-2)/3;
%         sub=0;
%         tmp=[];
%         for k=0:nSig-1%umSignal-1
%             x=prmVec(k*3+1)+yCropMin-oX-1+ss;
%             y=prmVec(k*3+2)+xCropMin-oY-1+ss;
%             a=prmVec(k*3+3);
%             
%             dx=prmStd(k*3+1);
%             dy=prmStd(k*3+2);
%             da=prmStd(k*3+3);
%             % discard signals with too large uncertainty
%             if double(dx) > 1.0 || double(dy) > 1.0
%                 sub=sub+1;
%                 continue;
%             end                
%             
%             tmp=[tmp;x y a s dx dy da ds nSig bg];
%         end
%         
%         if nSig > sub
%             tmp(:,end)=nSig-sub;
%             features=[features;tmp];
%         end
%         
%     end
    %else
    %    features=[features; zeros(1,9)];
    %end
    
end

% filter features for unsuitable maxima
%imMax=max(img(:));
%imN=img/imMax;
%threshold=graythresh(imN)*imMax;
%indA=features(:,3) > threshold;
%indS=features(:,4) > 0.9*sigma & features(:,4) < 1.1*sigma;
%indX=features(:,5) < 1.0;
%indY=features(:,6) < 1.0;
%indF=indS & indX & indY;
%features=features(indF,:);
