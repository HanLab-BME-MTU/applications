function [ratioAve,localAve,bgAve,randRatioAve,ratioInd,localInd,randRatioInd,clusterDensity,cellIntensity,backgroundStd] = colocalMeasurePt2Cnt(radius,...
    randomRuns,detectionData,imagePt,imageCnt,maskingFile)
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



%% Analysis
    maskingFile = logical(maskingFile);
    a =1; %Too lazy to change right now
    % Read in continuum image
    ICnt = double(imageCnt);
    orgCnt = double(imageCnt); 
    IPt = double(imagePt); 
    
    %Index points from punctate image
    xIndex = detectionData.xCoord(:,1);
    yIndex = detectionData.yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    [row, col] = find(maskingFile); %Verify
    maskList = [row, col];
    lia1 = ismember(round(QP),maskList,'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    QP(:,1) = lia1.*QP(:,1);
    QP(:,2) = lia1.*QP(:,2);
    QP( ~any(QP,2), : ) = [];
    
    %Create mask of points using QP points GOOD I THINK
    roundedQP = [round(QP(:,1)) round(QP(:,2))];
    localMask = zeros(size(ICnt));
    indexQP = sub2ind(size(localMask),roundedQP(:,1),roundedQP(:,2));
    localMask(indexQP) = 1;
    
        %% NEW Create masks to address non-uniform background and overgenerous
    %masking
    %new background subtraction
    corrValue = mean(ICnt(maskingFile==0));
%     corrValue = 1000;
    corrMask = corrValue*ones(size(ICnt,1),size(ICnt,2)); 
    
    compValue = mean(ICnt(maskingFile~=0));
    compMask = compValue*ones(size(ICnt,1),size(ICnt,2));%Should be ammenable to different sizes!
    nImage = filterGauss2D(ICnt,10);
    ICnt = ICnt-nImage;
    ICnt = ICnt+ compMask;
%     ICnt = ICnt- corrMask;
    
    corrValue = mean(IPt(maskingFile==0));
% %         corrValue = 1000;
    corrMask = corrValue*ones(size(IPt,1),size(IPt,2)); 
    compValue = mean(IPt(maskingFile~=0));
    compMask = compValue*ones(size(IPt,1),size(IPt,2));%Should be ammenable to different sizes!
    nImage = filterGauss2D(IPt,10);
    IPt = IPt-nImage;
    IPt = IPt+ compMask;
%     IPt = IPt- corrMask;

% %     meanMask = subImage>0;
% %     meanMaskComb = meanMask.*maskingFile;
% %     meanMaskErd = bwmorph(meanMaskComb,'clean');
    se = strel('disk',2);
    correctOGMask = imerode(maskingFile,se);
%     localMaskTest = localMaskTest.*correctOGMask;
% Use for individual detections
    [localMaskInd,~] = bwlabel(localMask);
    cc = connectedComponents.label2conncomp(localMaskInd);
    cc_dilated = connectedComponents.ccDilate(cc,strel('disk',2));%disk 2

    % Dilate to create local environment %---------------------------------
%     localMaskTest = bwdist(localMask)<=radius;
    se = strel('disk',radius);
    localMaskTest = imdilate(localMask,se);
    localMaskTest = logical(localMaskTest.*correctOGMask);
    % Everything in cell not detected
    cellBgMask = immultiply(~localMaskTest,double(maskingFile));
    cellBgCnt = immultiply(cellBgMask,ICnt); %Intensity of rest of cell in continuum channel
    cellBgPt = immultiply(cellBgMask,IPt); %Intensity of rest of cell in punctate channel
    backgroundStd = std(orgCnt(logical(cellBgMask)));
    %% Get cluster density in cell
    areaBG = find(maskingFile); %area of cell
    areaLM = find(localMask); %number of clusters
    density = length(areaLM)/length(areaBG);
    clusterDensity= (density/(0.0081));
    cellIntensityMask = orgCnt(maskingFile);
    cellIntensity = mean(cellIntensityMask(:));
    %% Back to orginal
% %     r=0.78;
% %     area = pi*r^2;
    %calculate average background intensity
    intensityBgAve = [mean(cellBgCnt(cellBgCnt~=0)) mean(cellBgPt(cellBgPt~=0))];
 
    %get the intensity around each detection
    tmp = regionprops(cc_dilated,ICnt,'MeanIntensity','Area');
    intensityLocalInd(:,1) = vertcat(tmp.MeanIntensity);
% %      tmp = regionprops(cc_dilated,ICnt,'PixelValues');
% %      mycell = struct2cell(tmp);
% %     mycell = mycell';
% %     myresult = cellfun(@(x) sum(x), mycell);
% %     intensityLocalInd(:,1) = myresult./area;   
% % %     areaInd = vertcat(tmp.Area);
    tmp = regionprops(cc_dilated,IPt,'MeanIntensity');
    intensityLocalInd(:,2) = vertcat(tmp.MeanIntensity);
% %      tmp = regionprops(cc_dilated,IPt,'PixelValues');
% %       mycell = struct2cell(tmp);
% %     mycell = mycell';
% %     myresult = cellfun(@(x) sum(x), mycell);
% %     intensityLocalInd(:,2) = myresult./area;   
    
    %normalize local intensities with background intensity
    intensityRatioInd(:,1) = intensityLocalInd(:,1)/intensityBgAve(a,1);
    intensityRatioInd(:,2) = intensityLocalInd(:,2)/intensityBgAve(a,2);
    
    %calculate average intensity and ratio
    intensityLocalAve = mean(intensityLocalInd);
    intensityRatioAve = mean(intensityRatioInd);
    %do a weighted average to account for different local area sizes
% % %     localWeights = areaInd/sum(areaInd);
% % %     intensityLocalAve = [intensityLocalInd(:,1)'*localWeights intensityLocalInd(:,2)'*localWeights];
% % %     intensityRatioAve = [intensityRatioInd(:,1)'*localWeights intensityRatioInd(:,2)'*localWeights];
    

    
        %Sample random positions to get intensities 
        %NOTE: Currently we will just try one run per image
        for j =1:randomRuns
    
            cellBg = datasample(maskList, length(roundedQP),'Replace',false);
            bgMask = zeros(size(ICnt,1),size(ICnt,2));
            indexBg = sub2ind(size(bgMask),cellBg(:,1),cellBg(:,2));
            bgMask(indexBg) = 1;
            
            [bgMaskInd,~] = bwlabel(bgMask);
            cc = connectedComponents.label2conncomp(bgMaskInd);
            cc_dilated = connectedComponents.ccDilate(cc,strel('disk',2));%disk
            
            %Determine avg intensity of random detection
%             bgMaskTest = bwdist(bgMask)<=radius;
            se = strel('disk',radius);
            bgMaskTest = imdilate(bgMask,se);
            bgMaskTest = logical(bgMaskTest.*correctOGMask);
            bgCellArea = immultiply(~bgMaskTest,double(maskingFile));
            bgCellAreaCnt = immultiply(bgCellArea,ICnt); %Use this to get bgCnt
            bgCellAreaPt = immultiply(bgCellArea,IPt); %Use this to get bgPt
            randIntensityBgAve(a,:) = [mean(bgCellAreaCnt(bgCellAreaCnt~=0)) mean(bgCellAreaPt(bgCellAreaPt~=0))];
            
            tmp = regionprops(cc_dilated,ICnt,'MeanIntensity','Area'); %CT
% %             tmp = regionprops(cc_dilated,ICnt,'PixelValues');
% %             mycell = struct2cell(tmp);
% %             mycell = mycell';
% %             myresult = cellfun(@(x) sum(x), mycell); 
% %             randIntensityInd(:,1) = myresult./area;
            randIntensityInd(:,1) = vertcat(tmp.MeanIntensity);
            
            tmp = regionprops(cc_dilated,IPt,'MeanIntensity','Area'); %PT
% %             tmp = regionprops(cc_dilated,IPt,'PixelValues');
% %             mycell = struct2cell(tmp);
% %             mycell = mycell';
% %             myresult = cellfun(@(x) sum(x), mycell); 
% %             randIntensityInd(:,2) = myresult./area;
            randIntensityInd(:,2) = vertcat(tmp.MeanIntensity);
            
            randIntensityRatioInd(:,1) = randIntensityInd(:,1)/randIntensityBgAve(a,1);
            randIntensityRatioInd(:,2) = randIntensityInd(:,2)/randIntensityBgAve(a,2);

            randRatioAve = mean(randIntensityRatioInd);
%             randRatioAve = mean(randIntensityRatioInd,1);
    
    
        end
ratioAve= intensityRatioAve;
localAve= intensityLocalAve; 
bgAve= intensityBgAve; 
ratioInd= intensityRatioInd;
localInd=intensityLocalInd ; 
randRatioInd= randIntensityRatioInd;        
        

% % %% Fitting Data
% % % Make optional
% % %test = cell2mat(intensityRatioInd);
% % %%TEMP: Remove image 13 from no TSP
% % % intensityRatioInd(13,:) = [];
% % % randIntensityRatioInd(13,:) = [];
% % figure;
% % s=length(intensityRatioInd);
% % c= linspace(1,s,s);
% % model = @(xm,a) a(1)*xm+ a(2);
% % % Remove outliers------------------------------------------------------
% % % %     test = cell2mat(intensityRatioInd(:,:));
% %     intensityRatioIndNew = intensityRatioInd;
% % % %     intensityRatioIndNew = cell(length(intensityRatioInd),2);
% % % % 
% % % %     %Discard outliers
% % % %     [~, D]= knnsearch(test,test,'K',5);
% % % %     allDist = D(:,2:end);
% % % %     avgDist = mean(allDist,2);
% % % %     [row,~,~] = find(avgDist>prctile(avgDist,95));
% % % %     test(row,:) = NaN;
% % % %     
% % % %     for i = 1: length(intensityRatioInd)
% % % %         intensityRatioIndNew{i,1}=test(1:length(intensityRatioInd{i,1}),1); 
% % % %         intensityRatioIndNew{i,2}= test(1:length(intensityRatioInd{i,1}),2);
% % % % % %         bgIntensityCell{i,1} = intensityBgAve(i,1)*ones(length(intensityRatioIndNew{i,1}),1);
% % % % % %         bgIntensityCell{i,2}=test(1:length(intensityRatioInd{i,1}),1);
% % % %         test(1:length(intensityRatioInd{i,1}),:)=[];
% % % %     end
% % % figure;
% % % k=1;
% % for k = 1:length(intensityRatioIndNew)
% %     test = cell2mat(intensityRatioIndNew(k,:));
% %     %Get rid of Nans
% %     keepInd = find(isnan(test(:,1)));
% %     test(keepInd,:)=[];
% %     sTest = sortrows(test,2);
% %     [Err, P] = fit_2D_data(sTest(:,2), sTest(:,1),'no');
% % %     scatter(sTest(:,2),sTest(:,1),'*');
% %     scatter(sTest(:,2),sTest(:,1),10,c(k)*ones(length(sTest),1));
% %     hold on
% % %     [ErrTLS,P] = numerFminS(model,2,[ 0 -10], [ 10 1], sTest(:,2), sTest(:,1));
% % %     [ErrTLS,P] = numerFminS(model,3,[-0.4 0.5 -1], [0.4 1.3 1], sTest(:,2), sTest(:,1))
% %     xInt = roots(P);
% %     f = polyval(P,sTest(:,2));
% %     f2 = polyval(P,0:0.1:3);
% % %     hold on
% % % %     plot(0:0.1:3,f2,'r');
% % % %     title(strcat('CD36-Actin +TSP  Fit:',num2str(P(1)),'x +',num2str(P(2)),'    X-Intercept: ',num2str(xInt(1))))
% % % %     ylabel('Actin Ratio')
% % % %     xlabel( 'CD36 Ratio')
% % % %     axis([0 2.5 0 2.5])
% %         y = sTest(:,1);
% %         res = y-f;
% %         sMin = sum(res.^2);
% %         D = length(y);
% %         bayInfoCriterion(k) = D*log(sMin/D)+ log(D)*2;
% % % %     %Alternative Fit
% % % %     xTest(:,1) =sTest(:,2);
% % % %     xTest(:,2) =sTest(:,1);
% % % %     [coeff,score,root] = pca(xTest);
% % % %     [n,p] = size(xTest);
% % % %     meanX = mean(xTest,1);
% % % %     Xfit1 = repmat(meanX,n,1) + score(:,1)*coeff(:,1)';
% % % %     dirVect = coeff(:,1);
% % % %     t = [min(score(:,1))-.2, max(score(:,1))+.2];
% % % %     endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];
% % % %     figure;
% % % %     plot(xTest(:,1),xTest(:,2),'bo');
% % % %     axis([0 2.5 0 2.5])
% % % %     hold on; plot(endpts(:,1),endpts(:,2),'k-');
% %     %-------------------------------------------------------------------
% %     ratioFit(k,:) = [P(1) P(2)];
% % end
% %     ylabel('Actin Ratio')
% %     xlabel( 'CD36 Ratio')
% %     axis([0 4 0 4])
% % % end

end






