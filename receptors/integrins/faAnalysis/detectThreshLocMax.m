function [detectedFeatures,pixelPos,maskBlobs] = detectThreshLocMax(image,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,...
    alphaLocMax,plotRes,mask,maxSize,maskThresh,alphaSubRes,psfSigma)
%detectThreshLocMax combines blob segmentation with local maxima detection
%
%SYNOPSIS maskComb = detectThreshLocMax(image,thresholdMethod,methodValue,...
%    filterNoise,filterBackground,minSize,alphaLocMax,plotRes,mask)
%
%INPUT  image     : 2D image to be segmented.
%       thresholdMethod: 
%                   'otsu' for Otsu.
%                   'rosin' for Rosin.
%                   'minmax' for first minimum after first maximum.
%                   'prct' to use a certain percentile.
%                   'user' to use a threshold input by user.
%                   Optional. Default: 'otsu'.
%       methodValue: Needed only if thresholdMethod = 'prct' or 'user'.
%                    If 'prct', then this is the percentile to use.
%                    Optional. Default: 90.
%                    If 'user', then this is the threshold value to use.
%                    Optional. Default: 0.9.
%       filterNoise: Either 0 to not filter noise or filter sigma > 0 to
%                    filter noise.
%                    Optional. Default: 1.
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 10.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       alphaLocMax: Alpha-value for judging the significance of local
%                   maxima.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskBlobs : Mask of blobs. 1 inside blobs, 0 outside.
%
%Khuloud Jaqaman January 2013

%% Output
maskComb = [];


%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%get image size
[imageSizeX,imageSizeY] = size(image);

%thresholding method
if nargin < 2 || isempty(thresholdMethod)
    thresholdMethod = 'Otsu';
end

%method value
if nargin < 3 || isempty(methodValue)
    switch thresholdMethod
        case 'prct'
            methodValue = 90;
        case 'user'
            methodValue = 0.9;
        otherwise
            methodValue = [];
    end
end

%noise filtering
if nargin < 4 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 5 || isempty(filterBackground)
    filterBackground = 10;
end

%minimum blob size
if nargin < 6 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 7 || isempty(plotRes)
    plotRes = 0;
end

%mask
if nargin < 8 || isempty(mask)
    mask = ones(imageSizeX,imageSizeY);
end

if ~logical(mask)
    error('Mask must be a logical image.');
end

%% Image preparation

%make sure that image is in double format
image = double(image);

%crop image
image = image .* mask;

%find nonzero values (due to masking)
nzInd = find(image);

%get minumum and maximum pixel values in image
minSignal = min(image(nzInd));
maxSignal = max(image(nzInd));

%normalize nonzero value between 0 and 1
imageNorm = zeros(size(image));
imageNorm(nzInd) = (image(nzInd) - minSignal) / (maxSignal - minSignal);


%estimate background by filtering image with a wide Gaussian
if filterBackground > 0
    imageBackground = filterGauss2D(imageNorm,filterBackground);
else
    imageBackground = zeros(imageSizeX,imageSizeY);
end

%remove background from image
imageMinusBackground = imageNorm - imageBackground;
 

%remove noise by filtering image-background with a narrow Gaussian
if filterNoise > 0
    imageMinusBackgroundFiltered = filterGauss2D(imageMinusBackground,filterNoise);
else
    imageMinusBackgroundFiltered = imageMinusBackground;
end
imageMinusBackgroundFilteredNorm = imageMinusBackgroundFiltered;
%estimate noise per pixel
% imageMinusBackgroundFiltered = imageMinusBackground - imageMinusBackgroundFiltered;

%% Thresholding

%estimate the intensity level to use for thresholding the image
switch thresholdMethod
    case 'otsu'
        try
            level = graythresh(imageMinusBackgroundFilteredNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (otsu, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'rosin'
        try
            [dummy,level] = cutFirstHistMode(imageMinusBackgroundFilteredNorm(nzInd),0); %#ok<ASGLU>
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (rosin, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'minmax'
        try
            level = thresholdFluorescenceImage(imageMinusBackgroundFilteredNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (minmax, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'prct'
        try
            level = prctile(imageMinusBackgroundFilteredNorm(nzInd),methodValue);
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (prct, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'user'
        try
            level = methodValue;
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (user, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
end

%threshold the image
imageThresholded = im2bw(imageMinusBackgroundFilteredNorm,level);

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

%go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
stats = regionprops(labels, 'Area'); %#ok<MRPBW>
idx = find([stats.Area] > minSize);

%output final blob mask from thresholding
maskBlobs = ismember(labels, idx);
testMask = regionprops(maskBlobs,'Area');
if find(vertcat(testMask.Area)>100*maxSize) %Quick fix, find better way to handle bad threshold and just exclusions
    if isempty(maskThresh)
        disp('First image failed to mask properly, check filtering size or method')
        return
    else
    maskBlobs = maskThresh;
    end
end
%Select centroids instead of local maxima for position 
%stats = regionprops(maskBlobs, 'Centroid');

%% Local maxima

% %estimate local background and noise statistics
[bgMean,bgStd] = ...
    spatialMovAveBG(imageNorm,imageSizeX,imageSizeY);

% %estimate background/noise statistics %% This was added BACK, TONY
% % bgIntDistr = imageMinusBackgroundFiltered(~maskBlobs);
% % [bgMean,bgStd] = robustMean(bgIntDistr);
% % bgStd = max(bgStd,eps);

%estimate background/noise statistics
% bgMean = 0;
% [~,bgStd] = robustMean(imageMinusBackgroundFiltered(~maskBlobs));

%call locmax2d to get local maxima in filtered image
fImg = locmax2d(imageMinusBackgroundFiltered,[3 3],1);



%get positions and amplitudes of local maxima
localMax1DIndx = find(fImg);
[localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);
testMaxAmp = image(localMax1DIndx);

% %get background values corresponding to local maxima%% This was added BACK, TONY
bgMeanMax = bgMean(localMax1DIndx);
bgStdMax = bgStd(localMax1DIndx);

%calculate the p-value corresponding to the local maxima's amplitudes
%assume that background intensity is normally
%distributed with mean bgMeanMax and standard deviation bgStdMax
pValue = 1 - normcdf(localMaxAmp,bgMeanMax,bgStdMax);


%calculate the threshold to distinguish significant local maxima
[~,threshLocMax] = cutFirstHistMode(localMaxAmp,0);
%------------------------------------------------------------------------
% Keep indices either normally or using bgImageDir 
% % if bgImageDir
% %     bgImg = double(imread(bgImageDir));
% % % %     nzInd = find(bgImg);
% % % %     bgImgNorm = zeros(size(bgImg));
% % % %     bgImgNorm(nzInd) = (bgImg(nzInd) - minSignal) / (maxSignal - minSignal);
% % 
% % % %     [bgMeanAbs,bgStdAbs] = robustMean(bgImg(:));
% % % %     bgStdAbs = max(bgStdAbs,eps);
% %     pValueAbs = 1 - normcdf(testMaxAmp,bgMeanAbs,bgStdAbs);
% %     indxKeep = pValue < alphaLocMax & pValueAbs < alphaLocMaxAbs;
% % else
%retain only those maxima with significant amplitude
indxKeep = pValue < alphaLocMax;
% % end
%------------------------------------------------------------------------
% % indxKeep = localMaxAmp > threshLocMax;
localMax1DIndx = localMax1DIndx(indxKeep);
localMaxAmp = localMaxAmp(indxKeep);
localMaxPosX = localMaxPosX(indxKeep);
localMaxPosY = localMaxPosY(indxKeep);

%make a mask image from the local maxima
maskLocMax = zeros(imageSizeX,imageSizeY);
maskLocMax(localMax1DIndx) = 1;
% SE = strel('disk',1);
% maskLocMax = imdilate(maskLocMax,SE);

%% Final mask

maskComb = maskBlobs | maskLocMax;

%% Work Space---------------------------------------------------
% Sub-Res Section----------------------------------
numLocMax = length(localMaxPosX);
cands = repmat(struct('status',1,'IBkg',[],...
'Lmax',[],'amp',[],'pValue',[]),numLocMax,1);
bgMeanMax = bgMeanMax(indxKeep);
pValue = pValue(indxKeep);
for iMax = 1 : numLocMax
cands(iMax).IBkg = bgMeanMax(iMax);
cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
cands(iMax).amp = localMaxAmp(iMax);
cands(iMax).pValue = pValue(iMax);
end
testAlpha = struct('alphaR',alphaSubRes,'alphaA',alphaSubRes,'alphaD',alphaSubRes,'alphaF',0);
featuresInfo = detectSubResFeatures2D_V2(image, cands,psfSigma,testAlpha,[],0,[],[],mean(bgStd(:)));%0.59

testCoordOrg = [featuresInfo.yCoord(:,1) featuresInfo.xCoord(:,1)];
testX = round(featuresInfo.xCoord);
testY = round(featuresInfo.yCoord);
% testAmp = featuresInfo.amp;
testCoord = [testY(:,1),testX(:,1)];
% testSizeInfo = NaN(size(testCoord,1),size(testCoord,2));
c = unique(testCoord,'rows');
locmaxMapSubR = false(size(labels));
for i = 1:length(c)
locmaxMapSubR(c(i,1),c(i,2)) = true;
end
%---------------------------------------------
%First determines connected components found by segmentation
[labels,nLabels] = bwlabel(maskBlobs,4); %This will be the threshold image
s = regionprops(labels, imageMinusBackgroundFiltered, {'WeightedCentroid','PixelValues','BoundingBox',...
    'Eccentricity','Area','PixelList','MeanIntensity'});

%Create binary mask where true values represent coordinates of an ROI
N = length(localMaxPosY);
locmaxMap = false(size(labels));
for i =1:N
locmaxMap(localMaxPosX(i,1),localMaxPosY(i,1)) = true; 
end
% % figure;
% imshowpair(locmaxMap,maskBlobs);
% % hold on
%first test: Is there one or more LM in a single object?
countN =1;
countS =1;
maskMultLM = zeros(imageSizeX,imageSizeY); %mask with objects that contain many LM
maskSingLM = zeros(imageSizeX,imageSizeY); %mask with objects that contain single LM
maskMyst = zeros(imageSizeX,imageSizeY);
for k =1:nLabels
    
    cellBlob = (labels == k); %Loop through objects
    matches = locmaxMapSubR.*cellBlob;
    [row,col]  = find(matches);
    matchNum = length(row);
%     matches = locmaxMap(cellBlob); %Are there LM's in that object?
%     matchNum = sum(matches); %How many?
    
    if  matchNum > 1
    s(k).StandardDev = std(double(s(k).PixelValues));

    maskMultLM(cellBlob) = 1;
    testArrayN(countN:countN+length(row)-1,1) = row;
    testArrayN(countN:countN+length(row)-1,2) = col;
    elim = ismember(testCoord,testArrayN(countN:countN+length(row)-1,:),'rows');
    testCoord(elim,:) = [];
    testCoordOrg(elim,:) = [];
    countN = countN + length(row);
    elseif matchNum ==1
        testArrayS(countS,:) =[row,col];
        maskSingLM(cellBlob) = 1;
      s(k).StandardDev = 0;
        elim = ismember(testCoord,[row,col],'rows');
        testCoord(elim,:) = [];
        testCoordOrg(elim,:) = [];
     [~, edgePixel] = min(s(k).PixelValues);
      x = [round(s(k).WeightedCentroid(1)) s(k).PixelList(edgePixel,1)];
      y = [round(s(k).WeightedCentroid(2)) s(k).PixelList(edgePixel,2)];
      c = improfile(image, x,y);
      s(k).Spread = max(c)/min(c);
      countS = countS+1;
    else
      s(k).StandardDev = 0;
      [~, edgePixel] = min(s(k).PixelValues);
      x = [round(s(k).WeightedCentroid(1)) s(k).PixelList(edgePixel,1)];
      y = [round(s(k).WeightedCentroid(2)) s(k).PixelList(edgePixel,2)];
      c = improfile(image, x,y);
      s(k).Spread = max(c)/min(c);
       maskMyst(cellBlob) = 1;
    end                
    
end
maskMultLM = logical(maskMultLM);
maskSingLM = logical(maskSingLM);
maskMyst = logical(maskMyst);
spreadThreshold = 0;%mean(vertcat(s.Spread))-std(vertcat(s.Spread));

listVar =vertcat(s.StandardDev);
indxVar = find(listVar);        

%Remove maxima that are too close
for i = 1:length(indxVar)
    
    [r(:,1),r(:,2)] = find(labels == indxVar(i)); % all points in object, verify x and y
    lia1 = ismember(testArrayN,r,'rows');
    lmPos = find(lia1);
    pathMeasure = zeros(((length(lmPos)*(length(lmPos)-1))/2),1);
    point = zeros(length(lmPos),2);
    count = 1;
    %lmPoints = [testArray(lmPos,:);round(s(indxVar(i)).Centroid)];
    for j = 1:(length(lmPos)-1)
        for k = (j+1):length(lmPos)
            x = [testArrayN(lmPos(j),1) testArrayN(lmPos(k),1)];
            %x = [lmPoints(j,1) lmPoints(k,1)];
            y = [testArrayN(lmPos(j),2) testArrayN(lmPos(k),2)];
            %y = [lmPoints(j,2) lmPoints(k,2)];
            c = improfile(image, x,y);

             pathMeasure(count,1) = (max(c))/(min(c));   
            
            % hold on; improfile(image, x,y), grid on; 
            point(j,count+1)= pathMeasure(count,1);
            point(j,1)= image(testArrayN(lmPos(j),1),testArrayN(lmPos(j),2));
            point(k,count+1) = pathMeasure(count,1);
            point(k,1) = image(testArrayN(lmPos(k),1),testArrayN(lmPos(k),2));
            count = count+1;
        end
    end
    % Good var - bad var = #number of points we want
    
    checkList = pathMeasure <= spreadThreshold;
    checkList = checkList.*pathMeasure;
    [~,~,badPath] = find(checkList);
    if length(badPath) ==length(pathMeasure)
        testArrayN(lmPos,:) = testArrayN(lmPos,:).*zeros(length(lmPos),2);
        testArrayN(lmPos(1),:) = [round(s(indxVar(i)).WeightedCentroid(2)),round(s(indxVar(i)).WeightedCentroid(1))];
    else
        for m = 1:length(badPath)
          testA = ismember(point,badPath(m));
          [row, ~,~] = find(testA);
          if isempty(row) ==1
              continue
          end
          [dismiss,~,~] = find(point == min(point(row,1)));
          cleanUp = find(point(dismiss,:));
          point(:,cleanUp(2:end)) = 0;
          testArrayN(lmPos(dismiss),:) = testArrayN(dismiss,:).*zeros(length(dismiss),2);
          clear dismiss
        end
    end
    
    clear r pathVariance 
    %hold off;
end
% testArray(lmPos(discard),:) = [];

if ~isempty(indxVar)
    fixedTestArray = testArrayN;%round(newPos)];%round(newPosB)]; %Add centroid locations
    fixedTestArray( ~any(fixedTestArray,2), : ) = [];
    %Create new mask, dilate to PSF around LM
    linearInd = sub2ind(size(image),fixedTestArray(:,1),fixedTestArray(:,2));
    maskLM = zeros(imageSizeX,imageSizeY);
    maskLM(linearInd) = 1;
    SE = strel('disk',2);
    maskLMDil = imdilate(maskLM,SE);
    maskMultLM = logical(maskMultLM+maskLMDil);
    %Implement Watershed segmentation to correct/center positions
    imageQ = imcomplement(image);
    markers = zeros(size(image,1),size(image,2));
    for q = 1: length(fixedTestArray)
        markers(fixedTestArray(q,1),fixedTestArray(q,2)) = 1;
    end
    markers = logical(markers);
    % maskComb = bwdist(maskComb)<=1;
    wsInput = imimposemin(imageQ, markers);
    wsInput(~maskMultLM) = Inf;
    wsOutput = watershed(wsInput);
    % % postWS = regionprops(wsOutput, image, {'Centroid','Area','Eccentricity','PixelIdxList','MeanIntensity'});
    % % postWS(1) = [];

    wsOutputNoBg = immultiply(wsOutput,maskMultLM);
    BW = logical(wsOutputNoBg);
    % wsOutputNoBg(wsOutput == 1) = 0;
    % cc = connectedComponents.label2conncomp(wsOutputNoBg);
    % cc.PixelIdxList = cc.PixelIdxList(2:end);
    % cc.NumObjects = cc.NumObjects-1;
    % cc_dilated = connectedComponents.ccDilate(cc,strel('disk',1));
    rp = regionprops(BW,image,{'WeightedCentroid','Area','Eccentricity','PixelIdxList','MeanIntensity'});
    test = vertcat(rp.Area);
    testE = vertcat(rp.Eccentricity);
    rp(test > maxSize & testE <= 0.7) = [];
    test = vertcat(rp.Area);
    rp(test < minSize) = [];
end
SE = strel('disk',2);
%Getting information from sub res features not in maskBlobs
% subR = isnan(testSizeInfo(:,1));
subRInd = sub2ind(size(image),testCoord(:,1),testCoord(:,2));
imagePSF = filterGauss2D(image,psfSigma);%0.59
Vq = interp2(imagePSF,testCoord(:,2),testCoord(:,1));
rp2Intensity = Vq;
rp2Size = 13.*ones(length(Vq),1);
% [lmIndR,lmIndC] = find(maskLocMax);
% [IDX, ~]= knnsearch([lmIndR,lmIndC],testCoord,'K',1);
% maskTest = zeros(imageSizeX,imageSizeY);
% % lmInd = sub2ind([350,350],lmIndR(IDX),lmIndC(IDX));
% maskTest(subRInd) = 1;
%     cc = connectedComponents.label2conncomp(maskTest);
%     cc_dilated = connectedComponents.ccDilate(cc,strel('disk',2));
% % SE = strel('disk',2);
% % maskTestDil = imdilate(maskTest,SE);
% %Method to correct close LM's-----------
% % SE = strel('square',3);
% % maskTestDil2 = imdilate(maskTest,SE);
% % newMaskLM = (maskTestDil2.*maskLocMax);
% SE = strel('disk',2);
% maskTestDil = imdilate(maskTest,SE);
% %-------------------------------
% wsInput = imimposemin(imageQ, logical(maskTest));
% wsInput(~maskTestDil) = Inf;
% wsOutput = watershed(wsInput);
% wsOutputNoBg = immultiply(wsOutput,logical(maskTestDil));
% BW2 = logical(wsOutputNoBg);
% rp2 = regionprops(BW2,imageMinusBackgroundFiltered,{'WeightedCentroid','Area','Eccentricity','PixelIdxList','MeanIntensity'});
%Geting information from single Mask Blob features

%Create new mask, dilate to PSF around LM
linearInd = sub2ind(size(image),testArrayS(:,1),testArrayS(:,2));
maskLMS = zeros(imageSizeX,imageSizeY);
maskLMS(linearInd) = 1;
imageQ = imcomplement(image);
maskLMSDil = imdilate(maskLMS,SE);
maskSingLM = logical(maskSingLM+maskLMSDil);
wsInput = imimposemin(imageQ, logical(maskLMS));
wsInput(~maskSingLM) = Inf;
wsOutput = watershed(wsInput);
wsOutputNoBg = immultiply(wsOutput,logical(maskSingLM));
BW3 = logical(wsOutputNoBg);
% BW3 = BW3 + maskMyst;
rp3 = regionprops(BW3,image,{'WeightedCentroid','Area','Eccentricity','PixelIdxList','MeanIntensity'});
% testSizeInfo(subR,1) = vertcat(rp2.MeanIntensity);
% testSizeInfo(subR,2) = vertcat(rp2.Area);

clear fixedTestArray
% dispCentroid = round(vertcat(rp.Centroid));
if ~isempty(indxVar)
finalPositions = [vertcat(rp.WeightedCentroid);vertcat(rp3.WeightedCentroid);[testCoordOrg(:,2),testCoordOrg(:,1)]];
finalSize = [vertcat(rp.MeanIntensity),vertcat(rp.Area);vertcat(rp3.MeanIntensity),vertcat(rp3.Area);rp2Intensity,rp2Size];
finalEccentricity = [vertcat(rp.Eccentricity);vertcat(rp3.Eccentricity);ones(size(testCoordOrg,1),1)];
pixelPos(:,1) = {rp.PixelIdxList,rp3.PixelIdxList,subRInd};
else
 finalPositions = [vertcat(rp3.WeightedCentroid);[testCoordOrg(:,2),testCoordOrg(:,1)]];
 finalSize = [vertcat(rp3.MeanIntensity),vertcat(rp3.Area);rp2Intensity,rp2Size];
 finalEccentricity = [vertcat(rp3.Eccentricity);ones(size(testCoordOrg,1),1)]; 
 pixelPos(:,1) = {rp3.PixelIdxList,subRInd};
end
localMaxPosZ = zeros(length(finalPositions),1); %Find better way to ignore axis
varPosX = 0.5*ones(length(finalPositions),1);
varPosY = 0.5*ones(length(finalPositions),1);
varAmp = zeros(length(finalPositions),1);
%% Plotting

if plotRes

%     %figure 1: the different analysis steps
%     figure
% 
%     %subplot 1: original image
%     subplot(2,2,1)
%     imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
%     
%     %subplot 2: background image
%     subplot(2,2,2)
%     imshow(imageBackground,[prctile(imageBackground(:),1) prctile(imageBackground(:),99)]);
%     
%     %subplot 3: bandpass-filtered image
%     subplot(2,2,3)
%     imshow(imageMinusBackgroundFiltered,[prctile(imageMinusBackgroundFiltered(:),1) prctile(imageMinusBackgroundFiltered(:),99)])
%     
%     %subplot 4: blob edges + local maxima
%     subplot(2,2,4)
%     imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
%     hold on
%     maskBounds = bwboundaries(maskBlobs);
%     cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
%     plot(localMaxPosY,localMaxPosX,'go')
    
    %figure 2: final mask
% %     h = figure;
% %     imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
% %     hold on
% %     maskBounds = bwboundaries(labelmatrix(cc_dilated) > 0);
% %     cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
% %     scatter(dispCentroid(:,1), dispCentroid(:,2));
    %Record x, y positions
    %positions = regionprops(maskComb,'Centroid');
end
    detectedFeatures.xCoord = [finalPositions(:,1) varPosX];
    detectedFeatures.yCoord = [finalPositions(:,2) varPosY];%Understand this!
    detectedFeatures.zCoord = [localMaxPosZ localMaxPosZ];
    detectedFeatures.amp = finalSize;
    detectedFeatures.ecc = finalEccentricity;
    detectedFeatures.size = [finalSize(:,2) varAmp];
    
    %[positions(:).Centroid(2)];
    


%%-------------------------------------------------------------------------
end

%% ~~~ the end ~~~

