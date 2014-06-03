function [intensityRatio,localArea, randRatio,bgArea] = colocalMeasurePt2Cnt(radius,percent,channel,randomRuns, detectionDir, imageDir)
% COLOCALMEASUREPT2CNT measures colocalization for two channels where only one is punctate and the other is continuous
% Synopsis
% [intensityRatio,localArea, randArea,bgArea] = colocalMeasurePt2Cnt(radius,percent,randomRuns. detectionDir, imageDir)
%Input:
%   radius: radius of local area around detection points to be analysed 
%   
%
%   percent: Top percentage of detections to be used
%
%   randomRuns: the number of runs for simulating random detection points
%
%   detectionDir: path to directory where detection analysis is stored
%    Ex: '/home2/avega/Desktop/CD36 Fyn/+TSP-1/detections.mat'
%
%   imageDir: path to directory where image tiff files are stored
%    Ex: '/home2/avega/Desktop/CD36 Fyn/+TSP-1/
%
% Output:
%   intensityRatio: vector of the ratio of localArea to bgArea values for each image
% 
%   localArea: vector of the average intensity of local area surrounding detections for each image 
%  
%
%   randRatio: vector of ratio of average intensity of random detections
%   within cell to background intensity for each image
%
%   bgArea: vector of the average intensity of everything within cell not including
%   detections for each image



load(detectionDir)


% Define size of  results. CONSIDER PUTTING IN STRUCT
[localArea,intensityRatio, randRatio,bgArea]=deal(zeros(length(movieInfo),1));

if nargin < 5 || isempty(imageDir)
    [fName,dirName] = uigetfile('*.tif','specify first image in the stack');
else
    if iscell(imageDir)
        [fpath,fname,fno,fext]=getFilenameBody(imageDir{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(imageDir)
        [fpath,fname,fno,fext]=getFilenameBody(imageDir);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end   
      outFileList = getFileStackNames([dirName,fName]);
      
      numFiles = length(outFileList);  



 for a = 1:numFiles
    
   % Read in second channel to generate cell mask and first channel to test
   % intensity correlation
    image  = outFileList{a};
%     refImage = imread(image,1);
%     refImage = double(refImage);
    I = imread(image,channel);
    I = double(I);

    %Index points from punctate image
    xIndex = movieInfo(a).xCoord(:,1);
    yIndex = movieInfo(a).yCoord(:,1);
    QP = [yIndex xIndex];
 
 
     %Calculate Cell Mask  
     [~, maskList,mask] = calcStateDensity(I); 

    %Find detections in QP that lie inside maskList    
     lia1 = ismember(round(QP),maskList,'rows');   

    %Multipling lia (binary vector) by QP will take replace coord outside
    %boundary with zero, last line removes all zeros from vector
     QP(:,1) = lia1.*QP(:,1);
     QP(:,2) = lia1.*QP(:,2);
     QP( ~any(QP,2), : ) = [];


    %Create Mask using QP points
    roundedQP = [round(QP(:,1)) round(QP(:,2))];
    localMask = zeros(size(I,1),size(I,2));
    indexQP = sub2ind(size(localMask),roundedQP(:,1),roundedQP(:,2));
    localMask(indexQP) = 1;

    
    %% Filter Higher Intensity Detections
    filterImage = immultiply(localMask, refImage);
    [row,col,flocIntensities] = find(filterImage);
    top = prctile(flocIntensities,percent);
    topIntensities = flocIntensities >= top;
    col = col.*topIntensities;
    row = row.*topIntensities;
    newQP = [row,col];
    newQP( ~any(newQP,2), : ) = [];
    
    localMask = zeros(size(I,1),size(I,2));
    indexQP = sub2ind(size(localMask),newQP(:,1),newQP(:,2));
    localMask(indexQP) = 1;
    
        % Dilate to create local environment and multiply mask by both channels
    localMaskTest = bwdist(localMask)<=radius;
    
    %% Generate Masks for Analysis
    
    cellBgMask = immultiply(~localMaskTest,mask);% Everything in cell not detected
    localMask = immultiply(localMaskTest,I); %Intensities of detections in fyn
%     compareIntensity = immultiply(localMaskTest,refImage); %Same but for cd36
    
%     % Determine intensity variance of image
%     testVariance = immultiply(mask,I); %Whole cell, no BG
%     [~,~,totalIntensity] =find(testVariance); 
%     cellVariance(a,1) = var(totalIntensity);
    
    % Determine significance of local environment intensity---------------
    [~,~,locIntensities] = find(localMask);
    localArea(a,1) = mean(locIntensities);
%     deviation(a,1) = std(locIntensities);

    %Now find the intensity of everything in cell but not near molecule
    cellBg = immultiply(cellBgMask,I); %Intensity of rest of cell, comparison
    [~,~,bgIntensities] = find(cellBg);
    bgArea(a,1) = mean(bgIntensities);
    intensityRatio(a,1) = localArea(a,1)/bgArea(a,1);
    
    

    

    %Sample random positions to get intensities
    for j =1:randomRuns

        cellBg = datasample(maskList, length(roundedQP),'Replace',false);  
        bgMask = zeros(size(I,1),size(I,2));
        indexBg = sub2ind(size(bgMask),cellBg(:,1),cellBg(:,2));
        bgMask(indexBg) = 1;

    %Determine avg intensity of random detection
        bgMask = bwdist(bgMask)<=radius;
        bgMask = immultiply(bgMask,I);
        [~,~,intensities] = find(bgMask);
        randBgArea = immultiply(~bgMask,mask);
        randBgArea = immultiply(randBgArea,I);
        [~,~,randBgIntensities]= find(randBgArea);
        randBg =mean(randBgIntensities);
        randRatio(j,1) = mean(intensities);
        randRatio(j,1) = randRatio(j,1)/randBg; %Comparison to random sample

       
    end
    
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



