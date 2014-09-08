function [intensityRatioAve,intensityLocalAve,intensityBgAve,randRatioAve,...
    intensityLocalInd,intensityRatioInd,ratioFit,bayInfoCriterion, cellDensity ] = colocalMeasureMultiChannel(radius,...
    randomRuns,detectionFile,firstImageFileCnt,firstImageFilePt,firstImageFileSeg,channelCnt,channelPt,channelSeg, maskingFile,segmentFile)
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
%   firstImageFileSeg: Name, including full path, of first segmented image tiff file.
%    Optional. User will be prompted to choose file if not supplied. DEFINE
%    better later
%
%   channelCnt: Position of contiuum image channel in a multi-tiff file.
%               Use 1 if single-tiff file. Optional. Default" 1.
%
%   channelPt: Position of puctate image channel in a multi-tiff file.
%              Use 1 if single-tiff file. Optional. Default" 1.
%
%   channelSeg: Position of segmented image channel in a multi-tiff file.
%              Use 1 if single-tiff file. Optional. Default" 1.
%
%   maskingFile: Full path where cell mask and list of pixels in mask is
%   found.
%
%   segmentFile: Full path where segmented area of third channel is found.
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
%   intensityLocalInd: cell array containing intensity values of each
%   indiviual local environment (first column: Continuum channel, second column: Punctate Channel)
%
%   intensityRatioInd:cell array containing intensity ratios of each
%   indiviual local environment to corresponding background (first column: Continuum channel, second column: Punctate Channel)
%
%   ratioFit: slope and intercept of a linear fit where x values are
%   individual ratio values in the punctate channel and y values are individual
%   ratio values in the continuum channel
%
%   randRatioFit: slope and intercept of a linear fit where x values are
%   

%% Input

if nargin < 4 || isempty(firstImageFileCnt)
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

if nargin < 5 || isempty(firstImageFilePt)
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

if nargin < 6 || isempty(firstImageFileSeg)
    [fName,dirName] = uigetfile('*.tif','PLEASE SPECIFY FIRST SEGMENTED IMAGE IN STACK');
else
    if iscell(firstImageFileSeg)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFileSeg{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFileSeg)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFileSeg);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end
outFileListSeg = getFileStackNames([dirName,fName]);

if nargin < 7 || isempty(channelCnt)
    channelCnt = 1;
end

if nargin < 8 || isempty(channelPt)
    channelPt = 1;
end

if nargin < 9 || isempty(channelSeg)
    channelSeg = 1;
end

%% Analysis

load(detectionFile)
load(maskingFile)
load(segmentFile)
%initialize output variables
[intensityLocalAve,intensityBgAve,intensityRatioAve] = deal(zeros(numFiles,2));
randRatioAve = zeros(numFiles,randomRuns);
cellDensity = zeros(30,1);
%go over all images
for a = 1:numFiles
    
    % Read in continuum image
    imageCnt  = outFileListCnt{a};
    ICnt = double(imread(imageCnt,channelCnt));
    
    %read in punctate image
    imagePt = outFileListPt{a};
    IPt = double(imread(imagePt,channelPt));
    
    %Index points from punctate image
% %     %CHANGE THIS BACK ONCE DONE! a=>numFile or which ever has greatest
% exposure in cd36 channel
    xIndex = movieInfo(a).xCoord(:,1);
    yIndex = movieInfo(a).yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    lia1 = ismember(round(QP),maskList{a},'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    QP(:,1) = lia1.*QP(:,1);
    QP(:,2) = lia1.*QP(:,2);
    QP( ~any(QP,2), : ) = [];
    
    %Create mask of points using QP points
    roundedQP = [round(QP(:,1)) round(QP(:,2))];
    cellDensity(a) = length(roundedQP)/length(maskList{a}); 
    localMask = zeros(size(ICnt));
    indexQP = sub2ind(size(localMask),roundedQP(:,1),roundedQP(:,2));
    localMask(indexQP) = 1;
    
    % Filter Higher Intensity Detections
    %     filterImage = immultiply(localMask, imagePt);
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
%     localMaskTest = localMask;
    % Everything in cell not detected
    cellBgMask = immultiply(~localMaskTest,mask(:,:,a));
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
    

    
        %Sample random positions to get intensities 
        %NOTE: Currently we will just try one run per image
        for j =1:randomRuns
    
            cellBg = datasample(maskList{a}, length(roundedQP),'Replace',false);
            bgMask = zeros(size(ICnt,1),size(ICnt,2));
            indexBg = sub2ind(size(bgMask),cellBg(:,1),cellBg(:,2));
            bgMask(indexBg) = 1;
    
            %Determine avg intensity of random detection
            bgMask = bwdist(bgMask)<=radius;
            bgCellArea = immultiply(~bgMask,mask(:,:,a));
            bgCellAreaCnt = immultiply(bgCellArea,ICnt); %Use this to get bgCnt
            bgCellAreaPt = immultiply(bgCellArea,IPt); %Use this to get bgPt
            randIntensityBgAve(a,:) = [mean(bgCellAreaCnt(bgCellAreaCnt~=0)) mean(bgCellAreaPt(bgCellAreaPt~=0))];
            
            tmp = regionprops(bgMask,ICnt,'MeanIntensity','Area'); %CT
            randIntensityInd{a,1} = vertcat(tmp.MeanIntensity);

            tmp = regionprops(bgMask,IPt,'MeanIntensity','Area'); %PT
            randIntensityInd{a,2} = vertcat(tmp.MeanIntensity);
            
            randIntensityRatioInd{a,1} = randIntensityInd{a,1}/randIntensityBgAve(a,1);
            randIntensityRatioInd{a,2} = randIntensityInd{a,2}/randIntensityBgAve(a,2);
            
% %             randBgArea = immultiply(~bgMask,mask(:,:,a));
% %             randBgArea = immultiply(randBgArea,ICnt);
% %             [~,~,randBgIntensities]= find(randBgArea);%%
% %             randBg =mean(randBgIntensities);
% %             randRatioAve(a,j) = mean(intensities)/randBg; %Comparison to random sample
    
    
        end
    
    
end

%% Fitting Data
%test = cell2mat(intensityRatioInd);
%%TEMP: Remove image 13 from no TSP
% intensityRatioInd(13,:) = [];
% figure;
s=length(intensityRatioInd);
c= linspace(1,s,s);
model = @(xm,a) a(1)*xm+ a(2);
% Remove outliers------------------------------------------------------
    test = cell2mat(intensityRatioInd(:,:));
    intensityRatioIndNew = cell(length(intensityRatioInd),2);

    %Discard outliers
    [~, D]= knnsearch(test,test,'K',5);
    allDist = D(:,2:end);
    avgDist = mean(allDist,2);
    [row,~,~] = find(avgDist>prctile(avgDist,95));
    test(row,:) = NaN;
    
    for i = 1: length(intensityRatioInd)
        intensityRatioIndNew{i,1}=test(1:length(intensityRatioInd{i,1}),1); 
        intensityRatioIndNew{i,2}= test(1:length(intensityRatioInd{i,1}),2);
% %         bgIntensityCell{i,1} = intensityBgAve(i,1)*ones(length(intensityRatioIndNew{i,1}),1);
% %         bgIntensityCell{i,2}=test(1:length(intensityRatioInd{i,1}),1);
        test(1:length(intensityRatioInd{i,1}),:)=[];
    end
figure;
k=1;
% for k = 1:length(intensityRatioIndNew)
    test = cell2mat(intensityRatioIndNew(:,:));
    %Get rid of Nans
    keepInd = find(isnan(test(:,1)));
    test(keepInd,:)=[];
    sTest = sortrows(test,2);
    [Err, P] = fit_2D_data(sTest(:,2), sTest(:,1),'no');
%     scatter(sTest(:,2),sTest(:,1),'*');
    scatter(sTest(:,2),sTest(:,1),10,c(k)*ones(length(sTest),1));
    hold on
%     [ErrTLS,P] = numerFminS(model,2,[ 0 -10], [ 10 1], sTest(:,2), sTest(:,1));
%     [ErrTLS,P] = numerFminS(model,3,[-0.4 0.5 -1], [0.4 1.3 1], sTest(:,2), sTest(:,1))
    xInt = roots(P);
    f = polyval(P,sTest(:,2));
    f2 = polyval(P,0:0.1:3);
%     hold on
    plot(0:0.1:3,f2,'r');
% %     title(strcat('CD36-Actin +TSP  Fit:',num2str(P(1)),'x +',num2str(P(2)),'    X-Intercept: ',num2str(xInt(1))))
% %     ylabel('Actin Ratio')
% %     xlabel( 'CD36 Ratio')
% %     axis([0 2.5 0 2.5])
        y = sTest(:,1);
        res = y-f;
        sMin = sum(res.^2);
        D = length(y);
        bayInfoCriterion(k) = D*log(sMin/D)+ log(D)*2;
% %     %Alternative Fit
% %     xTest(:,1) =sTest(:,2);
% %     xTest(:,2) =sTest(:,1);
% %     [coeff,score,root] = pca(xTest);
% %     [n,p] = size(xTest);
% %     meanX = mean(xTest,1);
% %     Xfit1 = repmat(meanX,n,1) + score(:,1)*coeff(:,1)';
% %     dirVect = coeff(:,1);
% %     t = [min(score(:,1))-.2, max(score(:,1))+.2];
% %     endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];
% %     figure;
% %     plot(xTest(:,1),xTest(:,2),'bo');
% %     axis([0 2.5 0 2.5])
% %     hold on; plot(endpts(:,1),endpts(:,2),'k-');
    %-------------------------------------------------------------------
    ratioFit(k,:) = [P(1) P(2)];
% end
    ylabel('Actin Ratio')
    xlabel( 'CD36 Ratio')
    axis([0 2.5 0 2.5])
% end

end



