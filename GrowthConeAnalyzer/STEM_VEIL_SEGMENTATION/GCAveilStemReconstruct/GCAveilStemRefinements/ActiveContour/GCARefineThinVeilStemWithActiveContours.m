function [ analInfo ] = GCARefineThinVeilStemWithActiveContours(img,veilStemMaskC,idxEnterNeurite,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% was testLowerScalesActiveContour until 20141207



    
    
    idxEnterNeurite = analInfo(frameList(iFrame)).idxEnterNeurite;
    
    
    img = double(imread( [listOfImages{frameList(iFrame),2} filesep listOfImages{frameList(iFrame),1}] ));
    [ny,nx] = size(img);
    [maskBack,backMu,backSig] = estimateBackgroundArea(img);
    scales = 1:10;
    [maxResLarge, maxThLarge ,maxNMSLarge ,scaleMapLarge]= multiscaleSteerableDetector(img,4,scales);
    
    
    maxNMSLarge = maxNMSLarge.*~maskBack ; % again only take responses with high fluorescence/image intensity
    values = maxNMSLarge(maxNMSLarge~=0);
    %         analInfo(iFrame).bodyEst.maxNMSLarge = maxNMSLarge;
    cutoff = prctile(values,50); % note arbitrarily set this cut-off to the 25th percentile
    ridgeCand = maxNMSLarge>cutoff;
    ridgeCand(ridgeCand>0)=1;
    
    ridgeCand = bwmorph(ridgeCand,'thin','inf');
    
    % combine all strong ridge candidates and the original mask
    %combine = ridgeCand | cMask;
    % take out those not connected
    %toErode = getLargestCC(combine);
    %toErode = logical(toErode);
    % imshow(toErode,[])
    % plot the scales measured
    
    setFigure(nx,ny,'off');
    %
    imagesc(ridgeCand.*scaleMapLarge);
    colorbar
    %
    if ~isdir([saveDir filesep 'Frames_Scales']);
        mkdir([saveDir filesep 'Frames_Scales']);
    end
    %
    saveas(gcf,[saveDir filesep 'Frames_Scales' filesep  num2str(frameList(iFrame),'%03d') '.tif']);
    saveas(gcf,[saveDir filesep 'Frames_Scales' filesep num2str(frameList(iFrame),'%03d') '.fig']); 
    % close gcf
    
    setFigure(nx,ny,'off');
    
    imshow(-img,[]);
    hold on
    ridgeCandLow = ridgeCand;
    ridgeCandHiBD = ridgeCand ;
    ridgeCandHiBD(scaleMapLarge<=4)=0; % filter out low scales
    ridgeCandLow(scaleMapLarge>4) = 0 ; % filter out high scales
    % dilate ridgeCandLow
    % ridgeCandHiDil = imdilate(ridgeCandHi,strel('disk',2));
    % Remove the current mask - assume it was accurate - combine in next
    % step. 
    ridgeCandHiBD(cMask==1)=0 ; 
    
    % Dilate the veil stem ridge detected
    ridgeCandHi = imdilate(ridgeCandHiBD,strel('disk',4)); % 
    
    
        
    % use an active contour method with a large smoothing factor- initiate
    % with the estimated dilated larger scale ridge region. 
    expand = activecontour(img,ridgeCandHi,100,'Chan-Vese',1);
    backbone=  analInfo(frameList(iFrame)).bodyEst.backbone;
    veilStem = (expand|cMask);
    makePlot =1 ; 
    if makePlot == 1
        imshow(-img,[]) 
        hold on 
        spy(ridgeCandHiBD,'r')
        roiYXC = bwboundaries(cMask); 
        cellfun(@(x) plot(x(:,2),x(:,1),'r'),roiYXC);
        roiYXInitiate = bwboundaries(ridgeCandHi); 
        cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYXInitiate); 
        roiYXContour = bwboundaries(expand); 
        cellfun(@(x) plot(x(:,2),x(:,1),'k'),roiYXContour); 
         if ~isdir([saveDir filesep 'Overlay']);
        mkdir([saveDir filesep 'Overlay']);
    end
        
        saveas(gcf,[saveDir filesep 'Overlay' filesep num2str(iFrame,'%03d') '.fig']); 
    end
    
    
    
    
    notBody = backbone.*~veilStem;
    
    %toErode = (ridgeCandHiDil | cMask | ridgeCandLow) ;
    nn = padarrayXT(double(notBody~=0), [1 1]);
    sumKernel = [1 1 1];
    nn = conv2(sumKernel, sumKernel', nn, 'valid');
    forJunctBreak = notBody~=0;
    % for now just make sure not f-ing up the original backbone
    %
    %forJunctBreak(backboneInfo(frameList(iFrame)).backboneSeedMask==1) = 0;
    %forJunctBreak(backboneInfo.bodyReconstruct.AfterConnect) = 0; %
    % remove anything that was originall in teh cleaned version (eventually
    % need to make it so that you are absolutely certain that there is no
    % residual junctions left in the cleaned ridge connect.
    % 2014 03 09-
    % 20141030 this might cause some problems in this capacity but I think
    % this was designed primarily to be for an iterative workup
    
    
    nn1 = (nn-1) .* forJunctBreak;
    junctionMask = nn1>2;
    % fix so that not taking out those previous connections you made
    
    
    
    
    
    
    backbone(cMask==1) = 0 ;
    backbone(expand==1) =0 ;
    backbone(junctionMask) = 0;
    
    
    toErode = (expand | cMask | backbone);
    
    
    
    
    
    [EPs,~,coords] = skel2graph2D(toErode);
    if ~isempty(vertcat(coords{:}))
        indEP = sub2ind([ny,nx],EPs(:,1),EPs(:,2));
        
        
        
        
        %idxSave = find(indEP == idxEnterNeurite); %% note sometimes bug here... should make so reiterate if this fails...
        %instead of doing find use intersect
        overlap = intersect(idxEnterNeurite,indEP);
        % save the entering neurite pieces from erosion
        if ~isempty(overlap)
            % maybe add to take it up or down a scale...but for now
            % just leave it
            idxSave = arrayfun(@(x) find(indEP == overlap(x)),1:length(overlap));
            idxSave = idxSave';
            
            EPs(idxSave,:) = [];
            coords(idxSave) = [];
        end
        
        
        
        if ~isempty(vertcat(coords{:}));
            coordBBOver = vertcat(coords{:});
            idxBBOver = sub2ind([ny,nx],coordBBOver(:,1),coordBBOver(:,2));
            idxEP = sub2ind([ny,nx],EPs(:,1),EPs(:,2));
            backbone([idxBBOver ; idxEP]) = 0;
        end
    end
    
    %% quick fix made 12-08-13
    %backbone2Dil = bwmorph(backbone2Dil,'bridge');
    backbone = bwmorph(backbone,'thin',inf);
    
    
    CCNotBody = bwconncomp(backbone);
    
    CCBody = bwconncomp(cMask);
    bodyLabels =  labelmatrix(CCBody);
    %
    for iCC = 1:numel(CCNotBody.PixelIdxList)
        % find the endpoints and test their labels
        % check if singleton
        if length(CCNotBody.PixelIdxList{iCC}) ==1
            idx = CCNotBody.PixelIdxList{iCC};
        else
            
            xy = getEndpoints(CCNotBody.PixelIdxList{iCC},size(notBody),0,0,4);
            
            idx = sub2ind(size(notBody),xy(:,2),xy(:,1));
        end
        test = zeros(size(notBody));
        test(idx) = 1;
        test = imdilate(test,strel('disk',2));
        labelsOverlap = bodyLabels(test==1);
        % test also if the piece entering the frame (ie overlapping with the
        % neurite body)
        overlapEnter = intersect(idxEnterNeurite,idx);
        
        % get rid of non overlapping labels
        labelsOverlap = labelsOverlap(labelsOverlap~=0) ;
        if (length(unique(labelsOverlap))==1 && isempty(overlapEnter)) ;
            % get rid of the body part
            backbone(CCNotBody.PixelIdxList{iCC})= 0;
            
        end
        
    end
    
    
    
    
    
    
    ridgeCandDilFinal = imdilate(backbone,strel('disk',2));
    %  spy(ridgeCandHi,'r');
    %  spy(ridgeCandLow,'b');
    cMaskFinal = (cMask | ridgeCandDilFinal| expand);
    cMaskFinal = logical(getLargestCC(cMaskFinal));
    % need to make sure to fill holes particularly important when finding
    % the neurite length- holes will tend to make it get stuck in the
    % calculation. 
    cMaskFinal = imfill(cMaskFinal,'holes'); 
    roiYXB = bwboundaries(cMask);
    
    roiYX = bwboundaries(cMaskFinal);
    cellfun(@(x) plot(x(:,2),x(:,1),'r','Linewidth',2),roiYXB);
    cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYX);
    
    if ~isdir([saveDir filesep  'Frames_Overlays'])
        mkdir([saveDir filesep 'Frames_Overlays']);
    end
    
    
    
    saveas(gcf,[saveDir filesep 'Frames_Overlays' filesep num2str(frameList(iFrame),'%03d') '.tif']);
    close gcf
    
    analInfo(frameList(iFrame)).masks.neuriteEdge = cMaskFinal; % replace
    
    %% also relabel the pixels associated with 'thick body'
    % this is unfortunately important just because of the way I set up
    % the filo reconstruct in the second step. I first label thick and
    % thin parts of the body so that we can collect filopodia from
    % these different regions. These definitions get a bit more hairy
    % when introduce the active contour.. ie regions at the tip could
    % be considered 'thin' regions by the scale measurement here -- as
    % they were eroded in the first step. For now just include them in
    % the 'thick regions' of body... as we are not really using these
    % definitions currently to delinate filopodia populations.
    % in the future you might want to work to reupdate more cleanly
    % this and thin pieces.. I think this might be done better though
    % in a step after the reconstruction is complete.. from the dist
    % trans information..
    % get new veil stem mask
    
    roiYXAll = bwboundaries(cMaskFinal);
    pixIndicesAll = sub2ind(size(cMaskFinal),roiYXAll{1}(:,1),roiYXAll{1}(:,2));
    
    % load the old thin body pixels
    pixIndicesThinBody =  analInfo(iFrame).bodyEst.pixIndThinBody;
    % anything for now that is not thin body is new thick body
    pixIndicesThickBodyNew = setdiff(pixIndicesAll,pixIndicesThinBody);
    
    thickBodyMask = zeros(size(cMaskFinal));
    thickBodyMask(pixIndicesThickBodyNew)=1;
    thickBodyMask = logical(thickBodyMask);
    % record
    analInfo(iFrame).masks.thickBodyMask = thickBodyMask;
    analInfo(iFrame).bodyEst.pixIndThickBody = pixIndicesThickBodyNew;
    
    
    save([saveDir filesep  'analInfoWithScaleSave.mat'],'analInfo');
end
end

