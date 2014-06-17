function [intensityRatioAve,intensityLocalAve,intensityBgAve,randRatioAve,...
    intensityLocalInd,intensityRatioInd] = colocalMeasurePt2Cnt(radius,percent,...
    randomRuns,detectionFile,firstImageFileCnt,firstImageFilePt,channelCnt,channelPt)
% COLOCALMEASUREPT2CNT measures colocalization for two channels where only one is punctate and the other is continuous
%
% Synopsis[intensityRatioAve,intensityLocalAve,intensityBgAve,randRatioAve,...
%    intensityLocalInd,intensityRatioInd] = colocalMeasurePt2Cnt(radius,percent,...
%    randomRuns,detectionFile,firstImageFileCnt,firstImageFilePt,channelCnt,channelPt)
%
%Input:
%
%   radius: radius of local area around detection points to be analysed.
%
%   percent: Top percentage of detections to be used.
%
%   randomRuns: the number of runs for simulating random detection points.
%
%   detectionFile: File name, including full path, where detection analysis is stored.
%    Ex: '/home2/avega/Desktop/CD36 Fyn/+TSP-1/detections.mat'
%
%   firstImageFileCnt: Name, including full path, of first continuum image tiff file.
%    Ex: '/home2/avega/Desktop/CD36 Fyn/+TSP-1/image_0001.tif'
%    Optional. User will be prompted to choose file if not supplied.
%
%   firstImageFilePt: Name, including full path, of first punctate image tiff file.
%    Optional. User will be prompted to choose file if not supplied.
%
%   channelCnt: Position of contiuum image channel in a multi-tiff file.
%               Use 1 if single-tiff file. Optional. Default" 1.
%
%   channelPt: Position of puctate image channel in a multi-tiff file.
%              Use 1 if single-tiff file. Optional. Default" 1.
%
% Output:
%   intensityRatioAve: vector of the ratio of intensityLocalAve to intensityBgAve values for each image
%
%   intensityLocalAve: vector of the average intensity of local area surrounding detections for each image
%
%   intensityBgAve: vector of the average intensity of everything within cell not including
%   detections for each image
%
%   randRatioAve: vector of ratio of average intensity of random detections
%   within cell to background intensity for each image
%
%   intensityLocalInd: 
%
%   intensityRatioInd:


%% Input

if nargin < 5 || isempty(firstImageFileCnt)
    [fName,dirName] = uigetfile('*.tif','PLEASE SPECIFY FIRST CONTINUUM IMAGE IN STACK');
else
    if iscell(firstImageFileCnt)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFileCnt{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFileCnt)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFileCnt);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end
outFileListCnt = getFileStackNames([dirName,fName]);
numFiles = length(outFileListCnt);

if nargin < 6 || isempty(firstImageFilePt)
    [fName,dirName] = uigetfile('*.tif','PLEASE SPECIFY FIRST PUNCTATE IMAGE IN STACK');
else
    if iscell(firstImageFilePt)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFilePt{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFilePt)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFilePt);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end
outFileListPt = getFileStackNames([dirName,fName]);

if nargin < 7 || isempty(channelCnt)
    channelCnt = 1;
end

if nargin < 8 || isempty(channelPt)
    channelPt = 1;
end

%% Analysis

load(detectionFile)

%initialize output variables
[intensityLocalAve,intensityBgAve,intensityRatioAve] = deal(zeros(numFiles,2));

%go over all images
for a = 1:numFiles
    
    % Read in continuum image
    imageCnt  = outFileListCnt{a};
    ICnt = double(imread(imageCnt,channelCnt));
    
    %read in punctate image
    imagePt = outFileListPt{a};
    IPt = double(imread(imagePt,channelPt));
    
    %Index points from punctate image
    xIndex = movieInfo(a).xCoord(:,1);
    yIndex = movieInfo(a).yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Calculate Cell Mask - MOVE TO OUTSIDE OF FUNCTION
    [maskList,mask] = calcCellBoundary(ICnt);
    
    %Find detections in QP that lie inside maskList
    lia1 = ismember(round(QP),maskList,'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    QP(:,1) = lia1.*QP(:,1);
    QP(:,2) = lia1.*QP(:,2);
    QP( ~any(QP,2), : ) = [];
    
    %Create mask of points using QP points
    roundedQP = [round(QP(:,1)) round(QP(:,2))];
    localMask = zeros(size(ICnt));
    indexQP = sub2ind(size(localMask),roundedQP(:,1),roundedQP(:,2));
    localMask(indexQP) = 1;
    
    % Filter Higher Intensity Detections
    %     filterImage = immultiply(localMask, refImage);
    %     [row,col,flocIntensities] = find(filterImage);
    %     top = prctile(flocIntensities,percent);
    %     topIntensities = flocIntensities >= top;
    %     col = col.*topIntensities;
    %     row = row.*topIntensities;
    %     newQP = [row,col];
    %     newQP( ~any(newQP,2), : ) = [];
    %
    %     localMask = zeros(size(I,1),size(I,2));
    %     indexQP = sub2ind(size(localMask),newQP(:,1),newQP(:,2));
    %     localMask(indexQP) = 1;
    
    % Dilate to create local environment
    localMaskTest = bwdist(localMask)<=radius;
    
    % Everything in cell not detected
    cellBgMask = immultiply(~localMaskTest,mask);
    cellBgCnt = immultiply(cellBgMask,ICnt); %Intensity of rest of cell in continuum channel
    cellBgPt = immultiply(cellBgMask,IPt); %Intensity of rest of cell in punctate channel
    
    %calculate average background intensity
    intensityBgAve(a,:) = [mean(cellBgCnt(cellBgCnt~=0)) mean(cellBgPt(cellBgPt~=0))];
    
    %get the intensity around each detection
    tmp = regionprops(localMaskTest,ICnt,'MeanIntensity','Area');
    intensityLocalInd{a,1} = vertcat(tmp.MeanIntensity);
    areaInd = vertcat(tmp.Area);
    tmp = regionprops(localMaskTest,IPt,'MeanIntensity');
    intensityLocalInd{a,2} = vertcat(tmp.MeanIntensity);
    
    %normalize local intensities with background intensity
    intensityRatioInd{a,1} = intensityLocalInd{a,1}/intensityBgAve(a,1);
    intensityRatioInd{a,2} = intensityLocalInd{a,2}/intensityBgAve(a,2);
    
    %calculate average intensity and ratio
    %do a weighted average to account for different local area sizes
    localWeights = areaInd/sum(areaInd);
    intensityLocalAve(a,:) = [intensityLocalInd{a,1}'*localWeights intensityLocalInd{a,2}'*localWeights];
    intensityRatioAve(a,:) = [intensityRatioInd{a,1}'*localWeights intensityRatioInd{a,2}'*localWeights];
    
    %     % Determine intensity variance of image
    %     testVariance = immultiply(mask,I); %Whole cell, no BG
    %     [~,~,totalIntensity] =find(testVariance);
    %     cellVariance(a,1) = var(totalIntensity);
    
    
    %JUST A PLACE-HOLDER FOR NOW
    randRatioAve = NaN;
    
    %     %Sample random positions to get intensities
    %     for j =1:randomRuns
    %
    %         cellBg = datasample(maskList, length(roundedQP),'Replace',false);
    %         bgMask = zeros(size(I,1),size(I,2));
    %         indexBg = sub2ind(size(bgMask),cellBg(:,1),cellBg(:,2));
    %         bgMask(indexBg) = 1;
    %
    %         %Determine avg intensity of random detection
    %         bgMask = bwdist(bgMask)<=radius;
    %         bgMask = immultiply(bgMask,I);
    %         [~,~,intensities] = find(bgMask);
    %         randBgArea = immultiply(~bgMask,mask);
    %         randBgArea = immultiply(randBgArea,I);
    %         [~,~,randBgIntensities]= find(randBgArea);
    %         randBg =mean(randBgIntensities);
    %         randRatioAve(j,1) = mean(intensities);
    %         randRatioAve(j,1) = randRatioAve(j,1)/randBg; %Comparison to random sample
    %
    %
    %     end
    
    %             minInt = min(refImage(:));
    %             maxInt = max(refImage(:));
    %             lim = [minInt maxInt];
    %             displayDet = immultiply(~localMaskTest,I);
    %             displayDetT = immultiply(~binImage,I);
    %             Fig2=figure;
    %             subplot(3,2,1); imagesc(refImage); title('Original CD36');
    %             subplot(3,2,2); imagesc(I); title('Original Fyn');
    %             subplot(3,2,3); imagesc(displayDet,lim); title('CD36 Detection');
    %             subplot(3,2,4); imagesc(mask); title('Cell Segmentation');
    %             subplot(3,2,5); imagesc(displayDetT,lim); title('Nicolas CD36 Detection');
    %             subplot(3,2,6); imagesc(mask2);  title('Nicolas Cell Segmenation');
    
    
end



