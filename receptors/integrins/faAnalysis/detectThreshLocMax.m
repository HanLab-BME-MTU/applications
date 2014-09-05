function [maskComb,imageMinusBackground,detectedFeatures,h,pixelPos] = detectThreshLocMax(image,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,...
    alphaLocMax,plotRes,mask, bgImageDir,alphaLocMaxAbs)
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
imageMinusBackground = [];
h= [];

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

%estimate background by filtering image with a wide Gaussian
if filterBackground > 0
    imageBackground = filterGauss2D(image,filterBackground);
else
    imageBackground = zeros(imageSizeX,imageSizeY);
end

%remove background from image
imageMinusBackground = image - imageBackground;

%remove noise by filtering image-background with a narrow Gaussian
if filterNoise > 0
    imageMinusBackgroundFiltered = filterGauss2D(imageMinusBackground,filterNoise);
else
    imageMinusBackgroundFiltered = imageMinusBackground;
end

%estimate noise per pixel
imageMinusBackgroundNoise = imageMinusBackground - imageMinusBackgroundFiltered;

%crop image
imageMinusBackgroundFiltered = imageMinusBackgroundFiltered .* mask;

%find nonzero values (due to masking)
nzInd = find(imageMinusBackgroundFiltered);

%get minumum and maximum pixel values in image
minSignal = min(imageMinusBackgroundFiltered(nzInd));
maxSignal = max(imageMinusBackgroundFiltered(nzInd));

%normalize nonzero value between 0 and 1
imageMinusBackgroundFilteredNorm = zeros(size(imageMinusBackgroundFiltered));
imageMinusBackgroundFilteredNorm(nzInd) = (imageMinusBackgroundFiltered(nzInd) - minSignal) / (maxSignal - minSignal);

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

%Select centroids instead of local maxima for position 
%stats = regionprops(maskBlobs, 'Centroid');

%% Local maxima

% %estimate local background and noise statistics
% [bgMean,bgStd] = ...
%     spatialMovAveBG(imageFilteredMinusBackground,imageSizeX,imageSizeY);

% %estimate background/noise statistics %% This was added BACK, TONY
bgIntDistr = imageMinusBackgroundFiltered(~maskBlobs);
[bgMean,bgStd] = robustMean(bgIntDistr);
bgStd = max(bgStd,eps);

%estimate background/noise statistics
% % bgMean = 0;
% % [~,bgStd] = robustMean(imageMinusBackgroundNoise(~maskBlobs));

%call locmax2d to get local maxima in filtered image
fImg = locmax2d(imageMinusBackgroundFiltered,[3 3],1);

%get positions and amplitudes of local maxima
localMax1DIndx = find(fImg);
[localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);

% %get background values corresponding to local maxima%% This was added BACK, TONY
% bgMeanMax = bgMean(localMax1DIndx);
% bgStdMax = bgStd(localMax1DIndx);

%calculate the p-value corresponding to the local maxima's amplitudes
%assume that background intensity is normally
%distributed with mean bgMeanMax and standard deviation bgStdMax
% pValue = 1 - normcdf(localMaxAmp,bgMeanMax,bgStdMax);
pValue = 1 - normcdf(localMaxAmp,bgMean,bgStd);

%calculate the threshold to distinguish significant local maxima
[~,threshLocMax] = cutFirstHistMode(localMaxAmp,0);
%------------------------------------------------------------------------
% Keep indices either normally or using bgImageDir 
if bgImageDir
    bgImg = double(imread(bgImageDir));
% %     bgMeanAbs = mean(bgImg(:));
% %     bgStdAbs = std(bgImg(:));
    [bgMeanAbs,bgStdAbs] = robustMean(bgImg(:));
    bgStdAbs = max(bgStdAbs,eps);
    pValueAbs = 1 - normcdf(localMaxAmp,bgMeanAbs,bgStdAbs);
    indxKeep = pValue < alphaLocMax & pValueAbs < alphaLocMaxAbs;
else
%retain only those maxima with significant amplitude
indxKeep = pValue < alphaLocMax;
end
%------------------------------------------------------------------------
% % indxKeep = localMaxAmp > threshLocMax;
localMax1DIndx = localMax1DIndx(indxKeep);
localMaxAmp = localMaxAmp(indxKeep);
localMaxPosX = localMaxPosX(indxKeep);
localMaxPosY = localMaxPosY(indxKeep);

%make a mask image from the local maxima
maskLocMax = zeros(imageSizeX,imageSizeY);
maskLocMax(localMax1DIndx) = 1;
SE = strel('square',3);
maskLocMax = imdilate(maskLocMax,SE);

%% Final mask

maskComb = maskBlobs | maskLocMax;

%% Work Space---------------------------------------------------


%First determines connected components found by segmentation
[labels,nLabels] = bwlabel(maskBlobs,4); %This will be the threshold image
s = regionprops(labels, image, {'Centroid','PixelValues','BoundingBox',...
    'Eccentricity','Area','PixelList'});

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

for k =1:nLabels
    
    cell = (labels == k); %Loop through objects
    matches = locmaxMap(cell); %Are there LM's in that object?
    matchNum = sum(matches); %How many?
    
    if  matchNum > 1
    s(k).StandardDev = std(double(s(k).PixelValues));
% %     text(s(k).Centroid(1),s(k).Centroid(2), ...
% %         sprintf('%2.1f', s(k).Eccentricity), ...
% %         'EdgeColor','b','Color','r');
% % if s(k).Eccentricity >= 0.9
% %     s(k).StandardDev = 0;

% % end

    else
      s(k).StandardDev = 0;
     [~, edgePixel] = min(s(k).PixelValues);
      x = [round(s(k).Centroid(1)) s(k).PixelList(edgePixel,1)];
      y = [round(s(k).Centroid(2)) s(k).PixelList(edgePixel,2)];
      c = improfile(image, x,y);
      s(k).Spread = max(c)/min(c);
    end                
    
end

spreadThreshold = mean(vertcat(s.Spread));%-std(vertcat(s.Spread));
listVar =vertcat(s.StandardDev);
indxVar = find(listVar);        
testArray(:,2) = localMaxPosX;
testArray(:,1) = localMaxPosY;
% discardPile = [];
%Remove maxima that are too close
for i = 1:length(indxVar)
    
    [r(:,2),r(:,1)] = find(labels == indxVar(i)); % all points in object, verify x and y
    lia1 = ismember(testArray,r,'rows');
    lmPos = find(lia1);
    pathMeasure = zeros(((length(lmPos)*(length(lmPos)-1))/2),1);
    point = zeros(length(lmPos),2);
    count = 1;
    %lmPoints = [testArray(lmPos,:);round(s(indxVar(i)).Centroid)];
    for j = 1:(length(lmPos)-1)
        for k = (j+1):length(lmPos)
            x = [testArray(lmPos(j),1) testArray(lmPos(k),1)];
            %x = [lmPoints(j,1) lmPoints(k,1)];
            y = [testArray(lmPos(j),2) testArray(lmPos(k),2)];
            %y = [lmPoints(j,2) lmPoints(k,2)];
            c = improfile(image, x,y);

             pathMeasure(count,1) = (max(c))/(min(c));   
            
            % hold on; improfile(image, x,y), grid on; 
            point(j,count+1)= pathMeasure(count,1);
            point(j,1)= image(testArray(lmPos(j),1),testArray(lmPos(j),2));
            point(k,count+1) = pathMeasure(count,1);
            point(k,1) = image(testArray(lmPos(k),1),testArray(lmPos(k),2));
            count = count+1;
        end
    end
    % Good var - bad var = #number of points we want
    
    checkList = pathMeasure <= spreadThreshold;
    checkList = checkList.*pathMeasure;
    [~,~,badPath] = find(checkList);
    if length(badPath) ==length(pathMeasure)
        testArray(lmPos,:) = testArray(lmPos,:).*zeros(length(lmPos),2);
        testArray(lmPos(1),:) = round(s(indxVar(i)).Centroid);
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
          testArray(lmPos(dismiss),:) = testArray(dismiss,:).*zeros(length(dismiss),2);
          clear dismiss
        end
    end
    
    
% %     numMolecules = length(find(pathMeasure > 1.25));
% %     if numMolecules == 0
% %         numMolecules = 1;
% %         sortedPoint = sort(point(:,2));
% %         keepMolecules = sortedPoint(end-(numMolecules-1):end);
% %         filterMolecules = ismember(point(:,2),keepMolecules);
% %         discard  = find(filterMolecules == 0);
% %         %testArray(lmPos(discard),:) = [];
% %         keep  = find(filterMolecules == 1);
% %         if length(keep) >1
% %         testArray(lmPos(keep(1)),:) = round(s(indxVar(i)).Centroid);
% %         testArray(lmPos([keep(2:end); discard]),:) = [];
% %         else
% %         testArray(lmPos(keep),:) = round(s(indxVar(i)).Centroid);
% %         testArray(lmPos(discard),:) = [];
% %         end
% %         
% %        
% %     else
% %     gbRatio = numMolecules-(length(pathMeasure)-numMolecules);
% %     if gbRatio <= 0
% %         numMolecules = 2;
% %     else
% %     numMolecules = round(sqrt(numMolecules*2));
% %     end
% %     sortedPoint = sort(point(:,1));
% %     keepMolecules = sortedPoint(end-(numMolecules-1):end);
% %     filterMolecules = ismember(point(:,1),keepMolecules);
% %     discard  = find(filterMolecules == 0);
% %     testArray(lmPos(discard),:) = []; 
% %     end
    clear r pathVariance 
    %hold off;
end
% testArray(lmPos(discard),:) = []; 
fixedTestArray = [testArray];%round(newPos)];%round(newPosB)]; %Add centroid locations
fixedTestArray( ~any(fixedTestArray,2), : ) = [];
%Implement Watershed segmentation to correct/center positions
imageQ = imcomplement(image);
markers = zeros(512,512);
for q = 1: length(fixedTestArray)
    markers(fixedTestArray(q,2),fixedTestArray(q,1)) = 1;
end
markers = logical(markers);
% maskComb = bwdist(maskComb)<=1;
wsInput = imimposemin(imageQ, ~maskComb | markers);
wsOutput = watershed(wsInput);
postWS = regionprops(wsOutput, image, {'Centroid','Area','Eccentricity','PixelIdxList','MeanIntensity'});
postWS(1) = [];
% % test = vertcat(postWS.Eccentricity);
% % discard = find(test>0.6);
% % postWS(discard) = [];

clear fixedTestArray
fixedTestArray = round(vertcat(postWS.Centroid));
objectAmp = round(vertcat(postWS.MeanIntensity));
localMaxPosZ = zeros(length(fixedTestArray),1); %Find better way to ignore axis
varPosX = 0.5*ones(length(fixedTestArray),1);
varPosY = 0.5*ones(length(fixedTestArray),1);
varAmp = zeros(length(fixedTestArray),1);
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
% % %     h = figure;
% % %     imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
% % %     hold on
%     maskBounds = bwboundaries(wsOutput);
%     cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
    %scatter(localMaxPosY, localMaxPosX);
% % %     scatter(fixedTestArray(:,1), fixedTestArray(:,2));
    %Record x, y positions
    %positions = regionprops(maskComb,'Centroid');
    detectedFeatures.xCoord = [fixedTestArray(:,1) varPosX];
    detectedFeatures.yCoord = [fixedTestArray(:,2) varPosY];%Understand this!
    detectedFeatures.zCoord = [localMaxPosZ localMaxPosZ];
    detectedFeatures.amp = [objectAmp varAmp];
    pixelPos(:,1) = {postWS.PixelIdxList};
    detectedFeatures.size = [vertcat(postWS.Area) varAmp];
    
    %[positions(:).Centroid(2)];
    


%%-------------------------------------------------------------------------
end

%% ~~~ the end ~~~

