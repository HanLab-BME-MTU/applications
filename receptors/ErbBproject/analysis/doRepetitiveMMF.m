function [ prmVec prmStd C res J status] = doRepetitiveMMF(img,prmVec0,varargin)
%DOREPETITIVEMMF fit Gaussian mixture model to 2D image data
%   This function tries to fit a Gaussian mixture model repetitively to
%   given image data. It first uses the obvious number of located maxima as
%   centroids for a Gaussian mixture model, and then tries to add one
%   Gaussian at a time to the current model in order to improve the fit.
%   This code is based on 'detectSubResFeatures2D'.
%  
%   Input:
%       img: image array
%       prmVec0: vector with initial parameter values
%       {'mode',mode}: any from 'xyAsc', default 'xyAc'
%       {'alphaR',alphaR}: alpha value for residual test, default 0.05
%       {'alphaA',alphaA}: alpha value for amplitude test, default 0.05
%
%   Output:
%       prmVec: vector of optimized parameter values
%       prmStd: error in optimized parameter values
%       C: variance-covariance matrix
%       res: structure with following fields
%           .data: residual
%           .pval: p value of Kolmogorov-Smirnov test (normal distributed)
%           .mean: mean of residuals
%           .std: standard deviation of residuals
%           .RSS: residual sum-of-squares
%       J: Jacobian
%
% 2011/07/27, US
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;
ip.addRequired('img',@isnumeric);
ip.addRequired('prmVec0',@isnumeric);

ip.addParamValue('mode','xyAc',@ischar);
ip.addParamValue('alphaA',0.05,@isscalar);
ip.addParamValue('alphaD',0.05,@isscalar);
ip.addParamValue('alphaF',0.05,@isscalar);
ip.addParamValue('alphaR',0.05,@isscalar);

ip.addOptional('doMMF',false,@islogical);

ip.parse(img,prmVec0,varargin{:});
mode=ip.Results.mode;
alphaA=ip.Results.alphaA;
alphaD=ip.Results.alphaD;
alphaF=ip.Results.alphaF;
alphaR=ip.Results.alphaR;

doMMF=ip.Results.doMMF;

% which parameters are present?
estIdx=regexpi('xyAsc', ['[' mode ']']);

% level for amplitude test (cf. fitGaussians2D.m)
kLevel=norminv(1.0-alphaA/2.0,0,1.0); 

% get dimensions of cropped region
[sizeX,~]=size(img);
% expected sigma of PSF
sig=prmVec0(4);

% by default, this routine returns failed signal
status=false;

% number of obvious potential maxima
nMax=numel(prmVec0(1:end-2))/3;
% number of valid pixels
nVal=sum(~isnan(img(:)));
% options for fitting routine
maxIter=500;
epsAbs=1e-6;
epsRel=1e-6;
fitOptions=[maxIter epsAbs epsRel];

% should a maxima been added to the current fit?
addP=doMMF;

prmVecOld=prmVec0;
% first fit with all maxima, isolated signals first
if nMax == 1
    [prmVecOld prmStdOld Cold resOld Jold]=...
        fitGaussian2D(img,prmVec0,mode,fitOptions);
else
% if nMax == 1
%     % do first anisotropic Gaussian fit
%     % see documentation for fitAnisoGaussian2D
%     prmVec0=[prmVec0(1:3) prmVec0(4) prmVec0(4) 0.0 prmVec0(5)];
%     [prmVecOld prmStdOld Cold resOld Jold]=...
%         fitAnisoGaussian2D(img,prmVec0,'xyarstc',fitOptions);
%     sx=prmVecOld(4);
%     sy=prmVecOld(5);
%     %phi=prmVecOld(6);
%     
%     % prepare prmVecOld for further processing or for return
%     % prmVecOld now contains: x, y, A, sigma and bg
%     prmVecOld=[prmVecOld(1:3) (sx+sy)/2.0 prmVecOld(end)];
%     
%     % check optimized parameters for size and symmetry
%     % add more potential maxima if the calculated PSF sigma (for x/y)
%     %  i) exceeds the theoretical PSF sigma considerably, or
%     % ii) deviates from symmetric signal
%     
%     % PSF sigma is too large
%     if sx > 1.5*sig || sy > 1.5*sig
%         prmVecOld(4)=sig;
%         addP=true;
%     else
%         % PSF sigma is symmetric?
%         if sx/sy > 0.8 && sx/sy < 1.2
%             prmVecOld(4)=sig;
%             [prmVecOld prmStdOld Cold resOld Jold]=...
%                 fitGaussian2D(img,prmVecOld,'xyAc',fitOptions);
%             addP=false;
%         else
%             addP=true;
%         end
%     end

    [prmVecOld prmStdOld Cold resOld Jold]=...
         fitGaussianMixture2D(img,prmVec0,mode,fitOptions);
end

% initial number of degrees of freedom
degFold=nVal-1-numel(prmVecOld)+1;
rssOld=resOld.RSS;

% get amplitudes for additional t-test
% n=(numel(prmVecOld)-2)/3;
% allAmp=NaN(n,1);
% allAmpStd=NaN(n,1);
% for k=1:n
%     allAmp(k)=prmVecOld(3*k);
%     allAmpStd(k)=prmStdOld(3*k);
% end
% tStat=allAmp./allAmpStd;
% pv=1.0-tcdf(tStat,degFold);
% indA=find(pv < alphaA);
% 
% if isempty(indA)
%     addP=false;
% else
%     % delete maxima that have failed test
%     indM=false(size(prmVecOld));
%     indM(end-1:end)=true;
%     for k=indA
%         indM(3*(k-1)+1)=true;
%         indM(3*(k-1)+2)=true;
%         indM(3*(k-1)+3)=true;
%     end
%     prmVecOld=prmVecOld(indM);
%     prmStdOld=prmStdOld(indM);
% end

while addP    
    % get location of maximal residual and corresponding raw data value
    [maxRes,indMaxRes]=max(resOld.data(:));
    [k1,k2]=ind2sub(size(img,1),indMaxRes);
    newMax=[k2-sizeX/2 k1-sizeX/2 maxRes];
    
    prmVecNew=[prmVecOld(1:end-2) newMax prmVecOld(end-1:end)];
    
    [prmVecNew prmStdNew Cnew resNew Jnew]=...
        fitGaussianMixture2D(img,prmVecNew,mode,fitOptions);
    
    degFnew=nVal-1-numel(prmVecNew)+1;
    rssNew=resNew.RSS;
    
    tStat=(rssNew/degFnew)/(rssOld/degFold);
    pVal=fcdf(tStat,degFnew,degFold);
    
    % repeat adding a maxima as long as F test fails
    if pVal < alphaR
        degFold=degFnew;
        rssOld=rssNew;
        prmVecOld=prmVecNew;
        prmStdOld=prmStdNew;
        Cold=Cnew;
        Jold=Jnew;
        resOld=resNew;
    else
        addP=false;
    end    
end

% discard all maxima that lie outside the fitting region and whose
% amplitude is not significantly above background

% number of fitted maxima
nMaxFit=(numel(prmVecOld)-2)/3;
% all fitted parameters go in one matrix
allPrm=reshape(prmVecOld(1:end-2),[3,nMaxFit])';
tmp=repmat(prmVecOld(end-1:end),nMaxFit,1);
allPrm(:,4:5)=tmp;

allStd=zeros(size(allPrm));
% reshape std vector to fit properly into allStd
% assume mode == 'xyac'
if nMaxFit > 1
  tmp=reshape(prmStdOld(1:end-2),[3,nMaxFit])';
  tmp(:,4)=prmStdOld(end);
  allStd(:,estIdx)=tmp;
else
    allStd(estIdx)=prmStdOld;
end

% does fitted maxima lie in image region?
imgSizeX=size(img,1)/2.0;
imgSizeY=size(img,2)/2.0;

idX=allPrm(:,1) > -imgSizeX & allPrm(:,1) < imgSizeX;
idY=allPrm(:,2) > -imgSizeY & allPrm(:,2) < imgSizeY;

id=idX & idY;
allPrm=allPrm(id,:);

% remove all maxima with negative amplitude
idAmp=allPrm(:,3) > 0.0;
allPrm=allPrm(idAmp,:);
allStd=allStd(idAmp,:);
% number of valid maxima
nMaxVal=size(allPrm,1);

idx=true(nMaxVal,1);

% test maxima for significance
A=allPrm(:,3);
dA=allStd(:,3);

sigma_r=resOld.std;
SE_sigma_r=resOld.std/sqrt(2*nVal-1);

df2=(nVal-1)*(dA.^2+(SE_sigma_r*kLevel)^2)./(dA.^4+(SE_sigma_r*kLevel)^4);
scomb=sqrt((dA.^2+(SE_sigma_r*kLevel)^2)/nVal);
T=(A-resOld.std*kLevel)./scomb;

pval_A=1.0-tcdf(T,df2);
hval_A=pval_A < alphaA;

%  results of Anderson-Darling test on resiudals
hval_AD=resOld.hAD;

rssVec=zeros(nMaxVal,7);
rssVec(:,1)=resOld.RSS;
rssVec(:,2)=sigma_r;
rssVec(:,3)=SE_sigma_r;
rssVec(:,4)=hval_AD;
rssVec(:,5)=pval_A;
rssVec(:,6)=hval_A;
rssVec(:,7)=nMaxVal;

%prmVec=prmVecOld;
%prmStd=prmStdOld;
prmVec=allPrm;
prmStd=[allStd,rssVec];
C=Cold;
res=resOld;
J=Jold;

status=true;

% if pVal < alphaR
%         degFold=degFnew;
%         rssOld=rssNew;
%         % maybe the following assignments are not necessary
%         %prmVecOld=prmVecNew;
%         %rmStdOld=prmStdNew;
%         %Cold=Cnew;
%         %resOld=resNew;
%         %Jold=Jnew;
%         
%         % get amplitudes for additional t-test
%         n=(numel(prmVecNew)-2)/3;
%         allAmp=NaN(n,1);
%         allAmpStd=NaN(n,1);
%         for k=1:n
%             allAmp(k)=prmVecNew(3*k);
%             allAmpStd(k)=prmStdNew(3*k);
%         end
%         tStat=allAmp./allAmpStd;
%         pv=1.0-tcdf(tStat,degFnew);
%         indA=find(pv < alphaA);
%         
%         if isempty(indA)
%             addP=false;
%         else
%             % delete maxima that have failed test
%             indM=false(size(prmVecNew));
%             indM(end-1:end)=true;
%             for k=indA
%                 indM(3*(k-1)+1)=true;
%                 indM(3*(k-1)+2)=true;
%                 indM(3*(k-1)+3)=true;
%             end
%             prmVecNew=prmVecNew(indM);
%             prmStdNew=prmStdNew(indM);
%         end
%         % also delete maxima that are outside of cropped region
% %         indM=false(size(prmVecNew));
% %         indM(end-1:end)=true;
% %         n=(numel(prmVecNew)-2)/3;
% %         for k=1:n
% %             xx=prmVecNew((k-1)*3+1);
% %             yy=prmVecNew((k-1)*3+1);
% %             
% %             if xx > -sizeX/2 && xx < sizeX/2 ...
% %                 && yy > -sizeY/2 && yy < sizeY/2
% %                 indM((k-1)*3+1)=true;
% %                 indM((k-1)*3+2)=true;
% %                 indM((k-1)*3+3)=true;
% %             end
% %         end
%         prmVecOld=prmVecNew;
%         prmStdOld=prmStdNew;
%         resOld=resNew;
%         Cold=Cnew;
%         Jold=Jnew;
%     else
%         addP=false;
%     end