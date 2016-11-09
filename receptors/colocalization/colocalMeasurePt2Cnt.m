function [ratioAve,localAve,bgAve,randRatioAve,ratioInd,localInd,randRatioInd,clusterDensity,cellIntensity] = colocalMeasurePt2Cnt(radius,...
    randomRuns,detectionData,imagePt,imageCnt,maskingFile)
% COLOCALMEASUREPT2CNT measures colocalization for two channels where only one is punctate and the other is continuous
%
% Synopsis: [intensityRatioAve,intensityLocalAve,intensityBgAve,randRatioAve,...
%    intensityLocalInd,intensityRatioInd] = colocalMeasurePt2Cnt(radius,percent,...
%    randomRuns,detectionFile,firstImageFileCnt,firstImageFilePt,channelCnt,channelPt)
%
%Input:
%
%   radius: radius of local area around detection points to be analysed.
%
%   randomRuns: the number of runs for simulating random detection points.
%
%   detectionData: movieInfo output from detection process containing x,y
%   positions of detected particles
%
%   imagePt: image of punctate channel
%
%   imageCnt: image of continuum channel
%
%   maskingFile: output of Masking process, binary matrix separating
%   background and foreground
%
%
% Output:
%   ratioAve: vector of the ratio of intensityLocalAve to intensityBgAve values for each image
%
%   localAve: vector of the average intensity of local area surrounding detections for each image
%
%   bgAve: vector of the average intensity of everything within cell not including
%   detections for each image
%
%   randRatioAve: vector of ratio of average intensity of random detections
%   within cell to background intensity for each image
%
%   ratioInd:cell array containing intensity ratios of each
%   individual local environment to corresponding background (first column: Continuum channel, second column: Punctate Channel)
%
%   localInd: cell array containing intensity values of each
%   individual local environment (first column: Continuum channel, second column: Punctate Channel)
%
%   randRatioInd:cell array containing randomized intensity ratios of each
%   individual local environment to corresponding background (first column: Continuum channel, second column: Punctate Channel)
%
%   clusterDensity: density of punctate objects in cell
%
%   clusterIntensity: intensity of punctate objects in cell
%
%
% Anthony Vega 09/2014   



%% Use mask to separate data
    maskingFile = logical(maskingFile);
    % Read in continuum and punctate image
    ICnt = double(imageCnt);
    orgCnt = double(imageCnt); 
    IPt = double(imagePt);

    %Index points from punctate image
    xIndex = detectionData.xCoord(:,1);
    yIndex = detectionData.yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    [row, col] = find(maskingFile); 
    maskList = [row, col];
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
    
    
    %% Image Processing
    % Background Subtraction and Uneven illumination correction
    corrValue = mean(ICnt(maskingFile==0));
    if isnan(corrValue)
        corrValue = 0;
    end
    corrMask = corrValue*ones(size(ICnt,1),size(ICnt,2)); 
    
    compValue = mean(ICnt(maskingFile~=0));
    compMask = compValue*ones(size(ICnt,1),size(ICnt,2));%
    nImage = filterGauss2D(ICnt,10);
    ICnt = ICnt-nImage;
    ICnt = ICnt+ compMask;
    ICnt = ICnt- corrMask;
    
    corrValue = mean(IPt(maskingFile==0));
    if isnan(corrValue)
        corrValue = 0;
    end
    corrMask = corrValue*ones(size(IPt,1),size(IPt,2)); 
    compValue = mean(IPt(maskingFile~=0));
    compMask = compValue*ones(size(IPt,1),size(IPt,2));
    nImage = filterGauss2D(IPt,10);
    IPt = IPt-nImage;
    IPt = IPt+ compMask;
    IPt = IPt- corrMask;
    
    %Slight erosion to avoid mask border
    se = strel('disk',2);
    correctOGMask = imerode(maskingFile,se);
    
    %% Create masks to further separate data
    % Use dilated objects for individual detections
    [localMaskInd,~] = bwlabel(localMask);
    cc = connectedComponents.label2conncomp(localMaskInd);
    cc_dilated = connectedComponents.ccDilate(cc,strel('disk',radius));

    % Dilate to create local environment %---------------------------------
    se = strel('disk',radius);
    localMaskTest = imdilate(localMask,se);
    localMaskTest = logical(localMaskTest.*correctOGMask);
    
    % Everything in cell not detected
    cellBgMask = immultiply(~localMaskTest,double(maskingFile));
    cellBgCnt = immultiply(cellBgMask,ICnt); %Intensity of rest of cell in continuum channel
    cellBgPt = immultiply(cellBgMask,IPt); %Intensity of rest of cell in punctate channel
    
    %% Data extraction from masks
    % Get cluster density and intensity in cell
    areaBG = find(maskingFile); %area of cell
    areaLM = find(localMask); %number of clusters
    clusterDensity = length(areaLM)/length(areaBG);
    cellNoBg = immultiply(correctOGMask,IPt);
    values = cellNoBg(cellNoBg~=0);
    cellIntensity = mean(values); %Mean cell intensity for punctate channel

    %calculate average background intensity
    intensityBgAve = [mean(cellBgCnt(cellBgCnt~=0)) mean(cellBgPt(cellBgPt~=0))];
 
    %get the intensity around each detection
    tmp = regionprops(cc_dilated,ICnt,'MeanIntensity','Area');
    intensityLocalInd(:,1) = vertcat(tmp.MeanIntensity);

    tmp = regionprops(cc_dilated,IPt,'MeanIntensity');
    intensityLocalInd(:,2) = vertcat(tmp.MeanIntensity);
    
    %normalize local intensities with background intensity
    intensityRatioInd(:,1) = intensityLocalInd(:,1)/intensityBgAve(1,1);
    intensityRatioInd(:,2) = intensityLocalInd(:,2)/intensityBgAve(1,2);
    
    %calculate average intensity and ratio
    intensityLocalAve = mean(intensityLocalInd);
    intensityRatioAve = mean(intensityRatioInd);
    
    
    %% Repeat procedures for randomized control 
        randIntensityRatioInd = [];
    for j =1:randomRuns

        cellBg = datasample(maskList, length(roundedQP),'Replace',false);
        bgMask = zeros(size(ICnt,1),size(ICnt,2));
        indexBg = sub2ind(size(bgMask),cellBg(:,1),cellBg(:,2));
        bgMask(indexBg) = 1;

        [bgMaskInd,~] = bwlabel(bgMask);
        cc = connectedComponents.label2conncomp(bgMaskInd);
        cc_dilated = connectedComponents.ccDilate(cc,strel('disk',2));%disk

        %Determine avg intensity of random detection
        se = strel('disk',radius);
        bgMaskTest = imdilate(bgMask,se);
        bgMaskTest = logical(bgMaskTest.*correctOGMask);
        bgCellArea = immultiply(~bgMaskTest,double(maskingFile));
        bgCellAreaCnt = immultiply(bgCellArea,ICnt); %Use this to get bgCnt
        bgCellAreaPt = immultiply(bgCellArea,IPt); %Use this to get bgPt
        randIntensityBgAve(1,:) = [mean(bgCellAreaCnt(bgCellAreaCnt~=0)) mean(bgCellAreaPt(bgCellAreaPt~=0))];

        tmp = regionprops(cc_dilated,ICnt,'MeanIntensity','Area'); %CT
        randIntensityInd(:,1) = vertcat(tmp.MeanIntensity)./randIntensityBgAve(1,1);

        tmp = regionprops(cc_dilated,IPt,'MeanIntensity','Area'); %PT
        randIntensityInd(:,2) = vertcat(tmp.MeanIntensity)./randIntensityBgAve(1,2);

        randIntensityRatioInd = [randIntensityRatioInd;randIntensityInd];
        clear randIntensityInd
%         randIntensityRatioInd(:,2) = [;randIntensityInd(:,2)/randIntensityBgAve(1,2)];

    end
 %% Store Results
randRatioAve = mean(randIntensityRatioInd);
ratioAve= intensityRatioAve;
localAve= intensityLocalAve; 
bgAve= intensityBgAve; 
ratioInd= intensityRatioInd;
localInd=intensityLocalInd ; 
randRatioInd= randIntensityRatioInd;        
        
end






