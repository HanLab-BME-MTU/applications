function [ maskAdhesion2 ] = refineAdhesionSegmentation(maskAdhesion, I, xNA, yNA, varargin)
%[ maskAdhesion ] = refineAdhesionSegmentation(maskAdhesion,I,xNA,yNA)
% chops off over-segmented segmentation (e.g. series of closely-located
% nascent adhesions) to several non-overlapping individual masks based on
% NA locations provided by [xNA yNA]. 
%   The function also filter out segmentation that has very weak signal.
% input:    maskAdhesion    mask of segmentation
%           I               raw image
%           xNA, yNA        x y coordinates of point sources

% Go with each island of segmentation
Adhs = regionprops(maskAdhesion,I,'Area','PixelList','PixelIdxList','MaxIntensity');
numAdhs = numel(Adhs);
numNAs=numel(xNA);
% NAsInAdhs=zeros(1,numAdhs);
%% Debug purpose
% for ii=1:numAdhs
%     curAdh = Adhs(ii);
%     % Check if curAdh contains only one NA
%     idxIn = false(1,numNAs);
%     for jj=1:numNAs
%         idxIn(jj)=any(round(xNA(jj))==curAdh.PixelList(:,1) & round(yNA(jj))==curAdh.PixelList(:,2));
%     end
%     curNumIn = sum(idxIn);
%     if curNumIn>1
%         disp([num2str(ii) ' has ' num2str(curNumIn) '.'])
%     end
% end
% iFiltered = filterGauss2D(I,1);

%% Divide the big segmentation if it contains more than one NAs
% tic
bwInteriorBDAll=false(size(maskAdhesion));
allIndex = cell(numAdhs,1);
% progressText(0,'watershedding adhesion mask')
% for ii=1:numAdhs
parfor ii=1:numAdhs
    curAdh = Adhs(ii);
    % Check if curAdh contains only one NA
    idxIn = false(1,numNAs);
    for jj=1:numNAs
        idxIn(jj)=any(round(xNA(jj))==curAdh.PixelList(:,1) & round(yNA(jj))==curAdh.PixelList(:,2)); %#ok<PFBNS>
    end
    curNumIn = sum(idxIn);
%     NAsInAdhs(ii)=curNumIn;
    tempMask=false(size(maskAdhesion));
    seedMask=false(size(maskAdhesion));
    if curNumIn>1 % there lies more than one NA in the current segmentation
        %% Chop it with watershed
        % Make own mask with NAs
        tempMask(curAdh.PixelIdxList)=true;
        tempMask=bwmorph(tempMask,'close');
        for kk=find(idxIn)
            seedMask(round(yNA(kk)),round(xNA(kk)))=true;
        end
        curDist = bwdist(seedMask);
        tempMaskDilate=bwmorph(tempMask,'dilate',1);
        curPreWS = tempMaskDilate.*(curDist);
%         curPreWS = tempMask.*(-iFiltered);

        curPreWS(~tempMaskDilate)=-Inf;

%         curWS = watershed(curPreWS,4);
        curWS = watershed(curPreWS,4);
        
        %% Now they are divided. Seeing if the interior boundary should be merged or not
        % Find the interior boundaries
%         B=bwboundaries(tempMask); % This should generate only one array of boundary
%         curBdIdx=B{1};
        % fixing current curWS with real boundary of tempMask
        curWS(~tempMask)=1;
        bwPerim=bwperim(tempMask);
        curWS(bwPerim)=0;
        % Now where is the interior boundaries?
        bwInteriorBD=~bwPerim & (curWS==0);
        % Get the pixels in exterior boundary that is adjacent to interior
        % ones, in order to clearly separate among interior parts.
        % Get the most exterior interior boundary pixel, then propagate it
        % is connected to two interior boundary pixels.
        
        periIntPixels = bwdist(bwPerim)<1.5 & bwInteriorBD;
        % Per each periIntPixel, see if it is connected to two
        % periIntPixels
        [rowPP,colPP] = find(periIntPixels);
        numPP=length(rowPP);
        distFromInt=bwdist(bwInteriorBD);
        for pp=1:numPP
            curRowPP=rowPP(pp);
            curColPP=colPP(pp);
            curNeigh=periIntPixels(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1);
            curNumNeiPP=sum(curNeigh(:));
            if curNumNeiPP<3 % This means it's not fully connected
                % cost matrix - 1. Closest neighbors get priority
                cropDistFromInt=distFromInt(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1);
                cropDistFromInt(cropDistFromInt==0)=3; % self and existing interior boundary has highest cost
                % 2. Pixel that is opposite side of existing interior
                % boundary gets priority
                cropIntPixels = bwInteriorBD(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1);
                cropIntPixelCost = ~cropIntPixels+1;
                cropIntPixelCost = rot90(cropIntPixelCost,2);
                % 3. If there is outside area in this neighboring box, it
                % has the highest priority
                outsidePixels = tempMask(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1);
                outsidePixelCost = outsidePixels+0.1;
                % 4. Exterior boundary definitely definitely has priority
                cropBwPerim=bwPerim(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1);
                cropBwPerimCost = ~cropBwPerim+0.5;
                % Cost
                costMatrix = cropDistFromInt.*cropIntPixelCost.*outsidePixelCost.*cropBwPerimCost;
                nPixelToConvert = 3-curNumNeiPP;
                [~,indSortedCost]=sort(costMatrix(:));
                cropIntPixels(indSortedCost(1:nPixelToConvert))=true;
                % Cancel interior boundary conversion if it is from outside
                % area criterion (#3)
                cropIntPixels(~outsidePixels)=false;
                bwInteriorBD(curRowPP-1:curRowPP+1,curColPP-1:curColPP+1)=cropIntPixels;
%                 bwInteriorBDAll(bwInteriorBD)=true;
            end
        end
        %% Local thresholding: I decided to do local thresholding for
        % individual interior boundaries 05/30/17
        inAdhPixels = I(tempMask); %#ok<PFBNS>
        intIndBDSeg=regionprops(bwconncomp(bwInteriorBD,8),I,'MeanIntensity','PixelIdxList');
        meanIntenIdvBDs=arrayfun(@(x) x.MeanIntensity,intIndBDSeg);
        thresLocal = mean(double(inAdhPixels))-0.3*std(double(inAdhPixels));
        
        idxSurvivedIdvBDs = meanIntenIdvBDs<thresLocal;
        % now have your true bwInteriorBD
        indInterPixIDs=arrayfun(@(x) x.PixelIdxList,intIndBDSeg(idxSurvivedIdvBDs),'uniformoutput',false);
        indInterPixIDsMat=cell2mat(indInterPixIDs);
        allIndex{ii}=(indInterPixIDsMat);
    else
        allIndex{ii}=[];
    end
%     progressText(ii/numAdhs)
end
% toc
%% combining all after parfor
allIndexAll=cell2mat(allIndex);
bwInteriorBDAll(allIndexAll)=true;
%% Go over each interior boundary and delete some of them if their intensity is significantly higher than background level
% % Look at the background intensity distribution
% inCelloutAdhMask = mask & ~maskAdhesion;
% bgInten = I(inCelloutAdhMask);
% % figure, histogram(bgInten)
% % There is clearly two distinct distribution. First one should be
% % the dark background near the actin-identified cell edges. The
% % second mode should be really distribution of cell inside. We need
% % to isolate this distribution.
% % opts = statset('Display','final','MaxIter',200);
% % try
% %     objBgInten = gmdistribution.fit(double(bgInten), 2, 'Options', opts);
% %     % Get the second distribution
% %     [~,smallerMu]=min(objBgInten.mu);
% %     muCytoInten = objBgInten.mu(smallerMu);
% %     stdCytoInten = sqrt(objBgInten.Sigma(:,:,smallerMu));
% %     thresInten = muCytoInten+4*stdCytoInten;
% % catch
%     try
%         thresInten = graythresh(double(bgInten)/double(max(bgInten)))*max(bgInten);
%         muCytoInten = mean(bgInten);
%         stdCytoInten = std(bgInten);
%     catch
%         try
%             [~, thresInten] = cutFirstHistMode(bgInten,0); %Rosin
%         catch
%             thresInten = quantile(bgInten,0.95); %Rosin
%             disp('95% Quantile of the background is used for interior boundary threshold.');
%         end
%         muCytoInten = mean(bgInten);
%         stdCytoInten = std(bgInten); %sqrt(objBgInten.Sigma(:,:,largerMu));
%     end
% % end
% %         posteriorBG=objBgInten.posterior(bgInten);
% %         maskDouble=double(mask);
% %         maskDouble(find(inCelloutAdhMask(:)))=posteriorBG(:,1);
%% Go really with interior bd segmentation
% intBDSeg=regionprops(bwconncomp(bwInteriorBDAll,4),I,'MeanIntensity','PixelIdxList');
% % Drop down bds whose intensity is bigger than thresInten
% meanIntenIntBDAll = arrayfun(@(x) x.MeanIntensity,intBDSeg);
% idxSurvivedBDs = meanIntenIntBDAll<thresInten;
% % now have your true bwInteriorBD
% allInterPixIDs=arrayfun(@(x) x.PixelIdxList,intBDSeg(idxSurvivedBDs),'uniformoutput',false);
% allInterPixIDs=cell2mat(allInterPixIDs);
% bwInteriorBD2 = false(size(mask));
% bwInteriorBD2(allInterPixIDs) = true;
% % bwInteriorBD2 = bwInteriorBDAll;
%% Make bwInteriorBD2 effective in maskAdhesion
maskAdhesion2 = maskAdhesion & ~bwInteriorBDAll;
%% Filter the weak segmenation
% % Look at the distribution of area and intensity
% survivedSegs=regionprops(bwconncomp(maskAdhesion2,4),I,'Area','MeanIntensity','PixelIdxList');
% 
% areaAll = arrayfun(@(x) x.Area, survivedSegs);
% intenAll = arrayfun(@(x) x.MeanIntensity, survivedSegs);
% meanArea =mean(areaAll);
% stdArea = std(areaAll);
% thresInten2 = muCytoInten+1*stdCytoInten;
% largeWeakAdhs = (areaAll>(meanArea+2*stdArea)) & (intenAll<thresInten2);
% allInterPixIDs=arrayfun(@(x) x.PixelIdxList,survivedSegs(~largeWeakAdhs),'uniformoutput',false);
% allInterPixIDs=cell2mat(allInterPixIDs);
% maskAdhesion3 = false(size(mask));
% maskAdhesion3(allInterPixIDs) = true;
% % figure, histogram(areaAll)
% % figure, plot(areaAll,intenAll,'ro')
% % curtosis or Thoman's spottyness
% % First drop down large area that has low intensity

% maskAdhesion3=maskAdhesion2;
end

