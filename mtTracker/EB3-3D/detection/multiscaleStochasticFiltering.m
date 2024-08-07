function [LM,scaleVol]=multiscaleStochasticFiltering(vol,XZRatio)
% Created by Philippe Roudot, 2017.  
% Variables cleared out to streamline memory usage by Kevin Dean, 2017.
%
%% Set Scale Range for Detection
% Adhesions setup - [1.1:0.1:4]
scales=[2:0.5:4];

% Myosin setup
% scales=[1.2:0.1:8];

%% Preallocate Cell Arrays and Measure Response in Parallel
LoGs=cell(1,length(scales));
LMLoc=cell(1,length(scales));
LMScaleIdx=cell(1,length(scales));
LMLoG=cell(1,length(scales));

parfor sIdx=1:length(scales)
    [mask, imgLM, imgLoG]=pointSourceStochasticFiltering(vol,[scales(sIdx) XZRatio*scales(sIdx)]);
    LoGs{sIdx}=mask.*imgLoG;
    lmIdx = find(imgLM~=0);
    [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    LMLoc{sIdx}=[lmy,lmx,lmz];
    LMScaleIdx{sIdx}=sIdx*(ones(length(lmIdx),1));
    LMLoG{sIdx}=imgLM(imgLM~=0);
end
clear mask imgLM imgLoG lmy lmx lmz XZRatio

%%
LMLoc=vertcat(LMLoc{:});
LMLoG=vertcat(LMLoG{:});
LMScaleIdx=vertcat(LMScaleIdx{:});

%%
R = 10;
D = createSparseDistanceMatrix(LMLoc, LMLoc, R);
bestLMScaleIdx=zeros(size(LMScaleIdx));

for lmIdx=1:size(D,1)
    localLMIdx=find(D(lmIdx,:)~=0);
    [~,bestLocalLMIdx]=max(LMLoG(localLMIdx));
    bestLMScaleIdx(localLMIdx(bestLocalLMIdx))=LMScaleIdx(localLMIdx(bestLocalLMIdx));
end
clear LMLoG LMScaleIdx LMScaleIdx lmIdx D localLMIdx bestLocalLMIdx R

%%
bestLMLoc=zeros(sum(bestLMScaleIdx~=0),3);
bestLMLoc(:)=LMLoc(bestLMScaleIdx~=0,:);
bestLMScaleIdx=bestLMScaleIdx(bestLMScaleIdx~=0);
clear LMLoc

%%
LM=zeros(size(vol));
LM(sub2ind(size(vol),bestLMLoc(:,1),bestLMLoc(:,2),bestLMLoc(:,3)))=bestLMScaleIdx;
clear bestLMScaleIdx bestLMLoc vol

%% Dilate LM
[x,y,z] = ndgrid(-4:4,-4:4,-3:3);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);
LM=imdilate(LM,se);

clear se x y z
%%
LoGPrevMax=LoGs{1};
LoGLocaLMax=zeros(size(LoGPrevMax));

scaleVol=ones(size(LoGPrevMax));
scaleVol(LoGPrevMax==0)=0;

for sIdx=2:(length(scales)-1)
    % Test if semi-global maximum (until this scale)
    [maxInd]=LoGPrevMax<LoGs{sIdx};
    LoGPrevMax(maxInd)=LoGs{sIdx}(maxInd);
    
    % Test if local maximum
    lmaxInd=(LoGs{sIdx}>LoGs{sIdx+1})&maxInd;   % is this local maximum ?
    LoGLocaLMax(lmaxInd)=LoGs{sIdx}(lmaxInd);
    scaleVol(lmaxInd)=sIdx;
end
LoGmask=LoGLocaLMax;
LM(LoGmask==0)=0;
clear LoGs LoGPrevMax LoGLocaLMax sIdx maxInd LoGPrevMax lmaxInd LoGmask scales

