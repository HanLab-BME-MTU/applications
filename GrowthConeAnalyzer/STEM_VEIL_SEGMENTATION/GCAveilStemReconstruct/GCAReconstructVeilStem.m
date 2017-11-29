function [ veilStem,TSFigs ] = GCAReconstructVeilStem(listOfImages,backboneInfo,BBScales,varargin)
% GCAReconstructVeilStem: (STEP III of GCA PACKAGE)
% This function reconstructs the veil/stem of the neurite.
% Veil here is broad protrusions typically associated with branched actin network.
% while Stem are thinner consolidated regions of neurite bridging larger amorphous veil pieces)
%
%% INPUT:
%
%   listOfImages: (REQUIRED) :
%
%   backboneInfo: (REQUIRED) :
%

%   veilStem: (OPTIONAL) : input as a cell array of a struct
%
%  (only relavent if need to restart)
%
% %% PARAMS: Local Thresholding
%    'OutputDirectory' (PARAM) : character array
%     Default: pwd
%
%    'LocalThresholdPatchSize' (PARAM) : Positive Scalar or path
%       Default:
%
%
%
% %% PARAMS: Morphological Opening %%
%
%    'DiskSizeLarge' (PARAM) : Positive scalar
%       that specifies the radius of the disk (in pixels) to be used for
%       the removal of thin objects from the input mask (ie the filopodia)
%       Larger values will remove thicker structures. Note for the lifeAct
%       channels that have very strong filopodia signal very often these
%       gradients for crossing filopodia tend not to be well segmented so
%       practically larger disk sizes are used in these cases.  If using a membrane
%       marker the filpodia often exhibit weak signal relative to
%       entire image and are often not segmented at all.
%       gcaMorphologicalOpeningWithGeometryConstraints.m
%       Default: 6
%
%    'DiskSizeSmall' (PARAM) : Positive scalar
%
%     Default = 3 :
%       See gcaMorphologicalOpeningWithGeometryConstraints.m
%
%    'TSMovie' (PARAM) : logical
%       Default : False
%       Used for paper
%
%     'MaxRadiusBridgeRidges' (PARAM) :  Positive scalar
%      Default : 5 (pixels)
%      The endpoints of the current veil/stem and the closest point along other
%      ridge candidates will be bridged. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%    veilStem : structure with fields
%            .finalMask : an rxc logical array (ie binary mask)
%
%
%            .rmHiIntPieces : an rxc logical array (ie binary mask)
%            .cycleFlag : character array documenting one of the following
%               'noCycles'
%               'parallelCyclesResolved'
%               'nonParallelCyclesResolved'
%               'unResolvedCycle'
%
%            .backbone
%            .idxEnterNeurite
%            .hashTag
%            .timeStamp
%
%% INPUTPARSER
%%Input check
ip= inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;
%REQUIRED
ip.addRequired('listOfImages');
ip.addRequired('backboneInfo');
ip.addRequired('BBScales');
%
ip.addOptional('veilStem',[]);
ip.addOptional('paramsArchived',[]);

% PARAMETERS

% Input/Output
ip.addParameter('OutputDirectory',[],@(x) ischar(x));

% Restart Options
ip.addParameter('StartFrame',1,@(x) isscalar(x));
nFrames = size(listOfImages,1);
ip.addParameter('EndFrame',nFrames,@(x) isscalar(x));

% Archiving the software version 
ip.addParameter('getGITHashTag',false);

% Troubleshoot Overlay flags and options 
ip.addParameter('TSOverlays',true,@(x) islogical(x));
ip.addParameter('TSMovie',false,@(x) islogical(x));
ip.addParameter('writeTitles',true);
ip.addParameter('screen2png',false);

% Detection of amorphous veil pieces : Initial 
ip.addParameter('maskDirectory',[]); % if empty perform local thresholding
    % else load the  masks from the input directory
ip.addParameter('LocalThresholdPatchSize',75,@(x) isscalar(x));

% Detection of amorphous veil pieces : Morphological Opening 
ip.addParameter('DiskSizeLarge',6,@(x) isscalar(x)); % in pixels: 
ip.addParameter('DiskSizeSmall',3,@(x) isscalar(x)); % in pixels:
    % see gcaMorphologicalOpeningWithGeometryConstraints.m

% Stem Radius Estimate  
ip.addParameter('reDilationRadius',4) ; % in pixels : see gcaResolveVeilStemCycles 

% Extend Path Search Via Bridging 
ip.addParameter('MaxRadiusBridgeRidges',5,@(x) isscalar(x));

ip.parse(listOfImages,backboneInfo,BBScales,varargin{:});
p = ip.Results;
if ~isempty(p.veilStem)
    veilStem = p.veilStem{1};
    paramsArchived = p.paramsArchived{1};
end

pToSave = rmfield(p,{'backboneInfo','listOfImages','veilStem','BBScales','paramsArchived'});
%%
if ~isempty(ip.Results.maskDirectory)
    listOfMasks = searchFiles('.tif',[],ip.Results.maskDirectory,0);
end
%% Start Loop
for iFrame = ip.Results.StartFrame:ip.Results.EndFrame
    figCount = 1;
    if p.TSMovie == true
        countMovie = 1;
    end
    %% Load information
    fileName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
    img = double(imread(fileName));
    
    backbone = backboneInfo(iFrame).backboneSeedMask;
    [ny,nx] = size(backboneInfo(1).backboneSeedMask);
    
    xyEnterNeurite = backboneInfo(iFrame).coordsEnterNeurite;
    idxEnterNeurite = sub2ind(size(img),xyEnterNeurite(:,2), xyEnterNeurite(:,1));
    cleanedRidge = backboneInfo(iFrame).linkedRidgesFinal;
    cleanedRidge(backbone==1) = 0;
    
    cleanedRidge =  bwmorph(cleanedRidge,'spur');
    veilStem(iFrame).bridged = false;
    
    %% Get Amorphous Large Scale Veil/Pieces
    
    if isempty(ip.Results.maskDirectory)
        % Perform local thresholding
        [~,maskForErod] = gcaThresholdOtsu_local(img, ip.Results.LocalThresholdPatchSize,3);
    else
        maskName = [char(listOfMasks(iFrame,2)) filesep char(listOfMasks(iFrame,1))];
        maskForErod = logical(imread(maskName));
    end
    
    % Perform a morphological opening with some geometric constraints.
    [erodForBody,saveMask]  =   gcaMorphologicalOpeningWithGeometryConstraints(maskForErod,p);
    %% TS Overlays: Show Mask Before Morphological Opening, After Morphological Opening, and the mask that 'saved'
    % using some geometry constraints.
%     if ip.Results.TSOverlays == true
%         
%         [ny,nx] = size(img);
%         TSFigs1(figCount).h = setFigure(nx,ny,'off');
%         TSFigs1(figCount).name = 'Morphological Opening';
%         imshow(-img,[]);
%         hold on
%         roiYX  = bwboundaries(maskForErod);
%         roiYX2 = bwboundaries(erodForBody);
%         
%         cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYX);
%         spy(saveMask,'r');
%         cellfun(@(x) plot(x(:,2),x(:,1),'color','b','Linewidth',2),roiYX2);
%         
%         text(5,5,['Removed Small Scale Signal- Disk Size ' num2str(ip.Results.DiskSizeLarge)],'color','g','FontSize',10);
%         text(5,20,'Final Veil/Stem Pieces Before Reconstruct With Ridge Info','color','b','FontSize',10);
%         text(5,40,{['Shrink Structuring Element Size From ' num2str(ip.Results.DiskSizeLarge) ' to '] ; ...
%             [ num2str(ip.Results.DiskSizeSmall) ' in Region Due to Geometry Constraints']},'color','r','FontSize',10)
%         figCount= figCount+1;
%     end % TS Overlay Morphological Opening
    %%  Delete veil pieces that are located along the edge of the frame that are not
    % overlapping with the main backbone entrance. (avoids cycles...)
    
    CCBodyPieces = bwconncomp(erodForBody);
    imSize = size(img);
    boundaryMask = zeros(imSize);
    boundaryMask(1:ny,1) =1;
    boundaryMask(1:ny,nx)=1;
    boundaryMask(1,1:nx)= 1;
    boundaryMask(ny,1:nx) =1;
    pixBound = find(boundaryMask==1);
    
    % find pieces that overlap with boundary but do not overlap with
    % the backbone seed
    pixBB = find(backbone==1);
    floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCBodyPieces.PixelIdxList);
    boundIdxCC = cellfun(@(x) ~isempty(intersect(pixBound,x)),CCBodyPieces.PixelIdxList);
    toKeep = ~(floatingBodyIdxCC & boundIdxCC);
    % save the pieces of boundIdxCC for later just in case there is
    % an error where it is needed later.
    boundaryBodyPixelIdx = CCBodyPieces.PixelIdxList(boundIdxCC & floatingBodyIdxCC);
    erodForBody = zeros([ny,nx]);
    erodForBody(vertcat(CCBodyPieces.PixelIdxList{toKeep}))=1;
    
    %% TS Figure: Removing Border Veil Node to Avoid Cycles
    toRemove = ~toKeep;
    if sum(toRemove)~=0
        removed = zeros([ny,nx]);
        removed(vertcat(CCBodyPieces.PixelIdxList{toRemove})) = 1;
        TSFigs1(figCount).h = setFigure(nx,ny,'off');
        TSFigs1(figCount).name = 'BorderPixelsRemoved';
        imshow(-img,[])
        hold on
        spy(removed,'r');
        text(5,5,'Border Node Removed To Avoid Cycles', 'color','r','FontSize',10);
        figCount = figCount +1;
    end
    %% Small Method Movie Flag: Show Backbone
    if p.TSMovie == true
        imgLarge = img;
        
        saveDirMov =  [ip.Results.OutputDirectory filesep 'SmallMovie' filesep 'Frame' num2str(iFrame,'%03d') ];
        if ~isdir(saveDirMov)
            mkdir(saveDirMov)
        end
        
        [nyLarge,nxLarge] = size(imgLarge);
        setFigure(nxLarge,nyLarge,'off');
        idx = find(backboneInfo(iFrame).backboneSeedMask==1);
        [yback,xback] = ind2sub([ny,nx],idx);
        
        % start by plotting the backbone coord and then showing
        % interations
        imshow(-imgLarge,[]);
        
        for i = 1:2
            if ip.Results.writeTitles
                hText = text(5,5,'Raw Image','Fontsize',10);
            end
            
            if ip.Results.screen2png
                helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            else
                saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            end
            
            countMovie = countMovie +1;
        end
        
        % save one with a scale bar
        pixSizeMic = 0.216;
        pixels = round(10/pixSizeMic);
        plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
        saveas(gcf,[saveDirMov filesep 'RawWithScaleBar' ,'.eps'],'psc2');
        
        close gcf
        setFigure(nxLarge,nyLarge,'off');
        imshow(-imgLarge,[]);
        hold on
        c = brewermap(2,'dark2');
        scatter(xback,yback,10,c(2,:),'filled');
        
        hold on
        
        scatter(xyEnterNeurite(:,1),xyEnterNeurite(:,2),100,'k','Marker','*');
        
        hold on
        
        %         pixSizeMic = 0.216;
        %         pixels = round(10/pixSizeMic);
        %         plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
        if ip.Results.writeTitles
            text(5,5,'Backbone', 'Color',c(2,:),'FontSize',10);
        end
        for i = 1:2
            if ip.Results.screen2png
                helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            else
                saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            end
            countMovie = countMovie +1;
        end
    end % end TSMovie
    %%  Start Main
    % Initiate the reconstruction of the body
    erodfilo = erodForBody;
    
    % get the connected components corresponding to the amorphous veil
    % portions
    CCAllBody = bwconncomp(erodfilo); %
    
    % Get Backbone
    pixBB = find(backbone==1);
    
    % Find the veil/stem nodes that do not intersect with the putative backbone
    floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
    
    % Take out the floating veil/stem nodes
    CC = CCAllBody; % copy cc structure and truncate to make the new body mask
    CC.NumObjects = CCAllBody.NumObjects -length(CCAllBody.PixelIdxList(floatingBodyIdxCC));
    CC.PixelIdxList(floatingBodyIdxCC) = [];
    
    % Veil/Stem Nodes For Mask
    newBodyMask= labelmatrix(CC);
    %% TSMovie Option: Visualize Veil/Stem Nodes Added to Mask
    if ip.Results.TSMovie
        
        roiYXBN =  bwboundaries(newBodyMask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','b','linewidth',2),roiYXBN);
        if ip.Results.writeTitles
            text(5,15,'Veil Added','Color','b','fontsize',10);
        end
        for i = 1:2
            if ip.Results.screen2png
                helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            else
                saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            end
            countMovie = countMovie+1;
        end
    end
    %% If it the first iteration and you did not add any Veil/Stem nodes to your mask
    % it can be due to one of two reasons
    % 1. you neruite is extremely thin and is best described as a ridge
    % 2. you had an extremely bad initial backbone estimation for that frame
    % Here we (non-elegantly) use some previous temporal information to find the correct path of the neurite
    % if there is evidence that we have a severely truncated
    % segmentation for a given frame
    
    if sum(newBodyMask(:))==0 % you didn't add anything backbone might just be truncated
        % check to see if frame 1: if frame 1 can't do temporal back
        % information on the fly
        if iFrame ~=1
            if sum(erodfilo(:))== 0% this is likely the problem you have nothing there!
                % try putting back the boundary body pixels : if not
                % just have a thin neurite (need to make small function
                % to address)
                erodfilo(vertcat(boundaryBodyPixelIdx{:}))=1;  % either that or you can not take out if problem (this might be why it would be good to set up as a cycle problem in end)
            end
            %try getting the veil/stem mask from the last frame
            % it is likely very similar
            oldBody = veilStem(iFrame-1).finalMask;
            pixOldBody= find(oldBody==1);
            
            % try again: find veil/stem nodes pieces that overlap
            CC = bwconncomp(erodfilo);
            floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixOldBody,x)),CC.PixelIdxList);
            CC.NumObjects = CC.NumObjects -length(CC.PixelIdxList(floatingBodyIdxCC));
            CC.PixelIdxList(floatingBodyIdxCC) = [];
            newBodyMask= labelmatrix(CC);
            
            % plot old body : old backbone/new backbone
            if ip.Results.TSOverlays == true
                TSFigs1(figCount).h = setFigure(nx,ny,'off');
                TSFigs1(figCount).name = 'Initial Scan For Truncations';
                imshow(-img,[]);
                hold on
                roiYXOld=  bwboundaries(oldBody);
                cellfun(@(x) plot(x(:,2),x(:,1),'b'),roiYXOld);
                spy(backbone);
                roiYXNew = bwboundaries(newBodyMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYXNew);
                figCount  = figCount +1;
            end
            
            % update the current backbone information as well since it
            % was likely poor
            backbone = backboneInfo(iFrame-1).backboneSeedMask;
            
            backboneInfo(iFrame).backboneSeedMask = backbone;  % just in case you need to use it again.
            
            xyEnterNeurite = backboneInfo(iFrame-1).coordsEnterNeurite;
            
            idxEnterNeurite = sub2ind(size(img),xyEnterNeurite(:,2), xyEnterNeurite(:,1));
            backboneInfo(iFrame).coordsEnterNeurite = xyEnterNeurite;
            cleanedRidge = backboneInfo(iFrame-1).linkedRidgesFinal;
            backboneInfo(iFrame).linkedRidgesFinal = cleanedRidge;
            cleanedRidge(backbone==1) = 0;
            
            cleanedRidge =  bwmorph(cleanedRidge,'spur');
        end
        veilStem(iFrame).trunc = 1;
    else
        veilStem(iFrame).trunc = 0;
        
    end %%   sum(NewBodyMask(:)) == 0 (test for complete truncation)
    
    %% Test to see if all floatingVeilStem pieces have been accounted for
    newBodyMask = newBodyMask>0;
    moreVeilStemNodes = sum(floatingBodyIdxCC);
    if  moreVeilStemNodes == 0
        veilStem(iFrame).rmHiIntPieces = [];
        %         moreVeilStemNodes = 1;
    end
    %% If there are more have it enter the reconstruct...
    if moreVeilStemNodes > 0
        floatingBodyIdxCCOld = floatingBodyIdxCC;
        iter = 1;
        reconstruct = 1;
        
        while reconstruct == 1
            %% Get All Ridges Overlapping to New/Veil Stem to Explore Potential Paths
            pixBodySave = find(newBodyMask==1);
            CCRidgeBone = bwconncomp(cleanedRidge);
            
            % Delete all ridges that fail to connect to the current iteration of the veil/stem
            idxNoIntersectRidgeCC = cellfun(@(x) isempty(intersect(pixBodySave,x)),CCRidgeBone.PixelIdxList);
            CCRidgeBoneConnect = CCRidgeBone;
            CCRidgeBoneConnect.PixelIdxList(idxNoIntersectRidgeCC) = [];
            CCRidgeBoneConnect.NumObjects = CCRidgeBoneConnect.NumObjects - sum(idxNoIntersectRidgeCC);
            
            CCCleanedRidge = CCRidgeBone;
            CCCleanedRidge.PixelIdxList(~idxNoIntersectRidgeCC) =[];
            CCCleanedRidge.NumObjects = CCRidgeBone.NumObjects -sum(~idxNoIntersectRidgeCC);
            
            pixRidgeConn = vertcat(CCRidgeBoneConnect.PixelIdxList{:});
            % update the backbone information...
            backbone(pixRidgeConn) = 1;
            labelMatCandRidge = labelmatrix(CCCleanedRidge);
            cleanedRidge = labelMatCandRidge>0;
            
            %% TSMovie : Exploring paths
            if ip.Results.TSMovie
                [yBack,xBack]=  ind2sub([ny,nx],pixRidgeConn);
                scatter(xBack,yBack,10,c(2,:),'filled') % 'color','r','linewidth',2);
                if ip.Results.writeTitles
                    text(5,35,'Explore New Ridge Paths','color',c(2,:),'Fontsize',10);
                end
                for i = 1:2
                    if ip.Results.screen2png
                        helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
                    else
                        saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
                    end
                    countMovie = countMovie+1;
                end
            end % if TSMovie
            %% Try to bridge small gaps in seed with with other nearby ridges/BodyCCs
            % Find end points of pixBB. use these as seed points
            % Find all points either ridge or floating body with in 2 pixels
            % of the endpoint. do this from the new body mask
            bridge = 1;
            if bridge == 1
                if sum(cleanedRidge(:))~=0 % if there are more cleaned ridge paths left... try to bridge them
                    % with endpoints of those overlapping with the current
                    % iteration of the veil/stem mask.
                    EPCoordsBodyMask = getEndpoints([find(newBodyMask==1) ; pixBB ; pixRidgeConn],[ny,nx]);
                    newBodyMaskPreBridge = newBodyMask;
                    backboneBeforeBridge = pixBB;
                    
                    if isempty(EPCoordsBodyMask)
                        
                        EPCoordsBodyMask = getEndpoints(pixBB,[ny,nx]);
                        idxEPs = sub2ind([ny,nx],EPCoordsBodyMask(:,2),EPCoordsBodyMask(:,1));
                        
                        coordsIn = backboneInfo(iFrame).coordsEnterNeurite;
                        idxIn = sub2ind([ny,nx],coordsIn(:,2),coordsIn(:,1));
                        EPCoordsBodyMask(idxEPs==idxIn,:)= [];
                    end
                    
                    idxPixRidgeCand = find(cleanedRidge==1);
                    [xyCandRidge(:,2), xyCandRidge(:,1)] = ind2sub([ny,nx],idxPixRidgeCand); %
                    
                    % first test for ridge connections
                    [idxRidge,dRidge]  = KDTreeBallQuery(xyCandRidge,EPCoordsBodyMask,ip.Results.MaxRadiusBridgeRidges); % for now just do ridges and link by
                    
                    if ~isempty(vertcat(idxRidge{:})) % if it found any point within range
                        % make a repmat for each query point the length of number of
                        % possible connections
                        E = arrayfun(@(i) [repmat(i, [numel(idxRidge{i}) 1]) idxRidge{i}], 1:length(EPCoordsBodyMask(:,1)), 'UniformOutput', false);
                        E = vertcat(E{:});
                        dRidge = vertcat(dRidge{:});
                        % give the idxRidges it's own nodes
                        % CCCleanedRidge = CCRidgeBone;
                        
                        % fix the node labels for matching
                        % want all col one to be 1:nendpoints testing while
                        % want to associate each idx in the second col with an independent node
                        % label larger than the endpoint.
                        NNodeQuery = size(EPCoordsBodyMask,1);  % get the
                        nodeLabels = E(:,1);
                        [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two separate EPs are competing over the same body point
                        
                        % Do this just in case two filo are competing over the same seed point
                        nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
                        EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
                        numberNodes = length(inputLinks) + NNodeQuery;
                        
                        D = max(dRidge)-dRidge;
                        % normalize
                        D= D./max(D);
                        
                        costTotal = (D);
                        
                        % set up for graph match edge number, idx of path, and associated cost
                        
                        E = [E costTotal D ];
                        
                        if ~isempty(EFinal)
                            
                            M = maxWeightedMatching(numberNodes, EFinal, costTotal);
                            E = E(M,:);% get those edges that matched (from original indexing)
                        end   % if ~isempty(EFinal)
                        
                        if ~isempty(E) % might be empty now if all were repeats
                            
                            % add linear segments corresponding to linked endpoints
                           
                            for iMatch = 1:size(E,1)
                                % get coords
                                xCoordQuery = EPCoordsBodyMask(E(iMatch,1),1);
                                yCoordQuery = EPCoordsBodyMask(E(iMatch,1),2) ;
                                xCoordInput = xyCandRidge(E(iMatch,2),1);
                                yCoordInput = xyCandRidge(E(iMatch,2),2);
                                
                                iseg{iMatch} = gcaBresenham([xCoordQuery yCoordQuery],...
                                    [xCoordInput yCoordInput]);
                                % save the coords coorsponding to the pixels you are adding
                                % need this to potentially tell not to wreck the junction
                                % mask
                                segSave{iMatch}= sub2ind([ny,nx], iseg{iMatch}(:,2), iseg{iMatch}(:,1));
                                backbone(sub2ind([ny,nx], iseg{iMatch}(:,2), iseg{iMatch}(:,1))) = 1;% add to mask
                                labels2Keep{iMatch} = labelMatCandRidge(yCoordInput,xCoordInput);
                                backbone(vertcat(CCCleanedRidge.PixelIdxList{labels2Keep{iMatch}})) = 1;
                                cleanedRidge(vertcat(CCCleanedRidge.PixelIdxList{labels2Keep{iMatch}})) = 0; % get ridge of in the cleaned ridge cand list
                            end % for iMatch
                            
                            % Mark frame as bridged : This can flag potentially
                            % lower confidence segmentation frames to check.
                            veilStem(iFrame).bridged = true;
                            %% Fix this TS Overlay for bridging: potentially use the figure in the
                            % TSMovie Overlay
                            if  ip.Results.TSOverlays == true
                                TSFigs1(figCount).h = setFigure(nx,ny,'off');
                                TSFigs1(figCount).name = 'Bridging';
                                imshow(backbone,[]);
                                hold on
                                scatter(EPCoordsBodyMask(:,1),EPCoordsBodyMask(:,2),'r','filled');
                                scatter(xyCandRidge(:,1),xyCandRidge(:,2),'y','filled');
                                figCount = figCount +1;
                            end
                            %%
                            if ip.Results.TSMovie
                                % plot the end-points
                                scatter(EPCoordsBodyMask(:,1),EPCoordsBodyMask(:,2),15,'k','filled');
                                
                                % plot the ridge candidates
                                hCand = scatter(xyCandRidge(:,1),xyCandRidge(:,2),5,c(2,:),'filled');
                                
                                saveAll = vertcat(iseg{:});
                                hLinks = scatter(saveAll(:,1),saveAll(:,2),10,'k','filled');
                                if ip.Results.screen2png
                                    helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
                                end
                                
                                countMovie = countMovie +1;
                                
                                delete(hCand)
                                
                                labels2KeepAll = horzcat(labels2Keep{:});
                                idxCandKeep = vertcat(CCCleanedRidge.PixelIdxList{labels2KeepAll});
                                [yCandKeep,xCandKeep] = ind2sub(size(img),idxCandKeep);
                                scatter(xCandKeep,yCandKeep,10,c(2,:),'filled');
                                if ip.Results.screen2png
                                    helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
                                end
                                countMovie = countMovie+1;
                            end  % TSMovie
                            %%
                        end % isempty(E)
                    end % isempty(vertcat(idxRidge{:}))
                    clear xyCandRidge  segSave
                    % New backbone complete
                end % if sumCleanedRidge
            end % if bridge
            %% Test to see if any of the new ridge paths overlap with other veil stem nodes
            backbone = bwmorph(backbone,'thin','inf');
            % get the pixels associated with the new backbone bone
            pixBB = find(backbone==1) ;
            
            % Get the logical indices of those body pieces that intersect
            % with the backbone
            floatingBodyIdxCCNew = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
            
            newBodyMask(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCCNew}))= 1;
            %% TSMovie
            if ip.Results.TSMovie
                
                maskAdd = zeros(ny,nx);
                maskAdd(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCCNew}))=1;
                roiYXVeilAdd = bwboundaries(logical(maskAdd));
                cellfun(@(x) plot(x(:,2),x(:,1),'color','b','linewidth',2),roiYXVeilAdd);
                
                if ip.Results.screen2png
                    helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
                end
                countMovie = countMovie+1;
                
            end % if TSMovie
            %% Test to see if either 1) all veil stem nodes have been accounted
            % for or 2) the extended backbone failed to find any path
            noAdd = sum(~(floatingBodyIdxCCNew == floatingBodyIdxCCOld));
            if  (noAdd==0  || sum(floatingBodyIdxCCNew)== 0)
                if sum(floatingBodyIdxCCNew)>0 % extra pieces not added save for potential refinement steps
                    % save as a mask- though might want to reconsider what is the most
                    % effective method of saving
                    rmHiIntPieces = zeros(ny,nx);
                    rmHiIntPieces(vertcat(CCAllBody.PixelIdxList{floatingBodyIdxCCNew}))=1;
                    veilStem(iFrame).rmHiIntPieces = rmHiIntPieces;     
                else      
                    veilStem(iFrame).rmHiIntPieces = [];
                    break
                end
                reconstruct = 0;
            end % if NoAdd
            floatingBodyIdxCCOld = floatingBodyIdxCCNew;
            iter = iter +1;
            %display(num2str(iter))
            if iter>100
                error('FiloTracker:NeuriteBodyEst:FailureToTerminate','Too Many Iterations Check Code');
            end
        end
    end % moreVeilStemNodes
    %% Remove RidgePieces With Free Ends
    newBodyXY = bwboundaries(newBodyMask);
    idxBody = cellfun(@(x) sub2ind(size(newBodyMask),x(:,1),x(:,2)),newBodyXY,'uniformoutput',0);
    idxBody = vertcat(idxBody{:});
    bodyNoFill = zeros(size(img));
    bodyNoFill(idxBody) = 1;
    
    notBody = double(backbone).*~double(newBodyMask);
    notBody = bwmorph(notBody,'thin','inf');
    
    %% Should no longer be relavent
    % This step was needed when I was bridging ridge  paths and forming new junctions
    %  in the backbone should
    % fix so that not taking out those previous connections you made
    if exist('saveJunct','var')==1
        juncCC = bwconncomp(junctionMask);
        
        pixAdded = vertcat(saveJunct{:});
        % do any of the junctions overlap with any of the added segments if not
        % save Though this might bite you in the ass as well if you added shit...
        % if have overlap with a segment you added previously don't remove that
        % junction.
        saveJunctIdx =cellfun(@(x) ~isempty(intersect(pixAdded,x)),juncCC.PixelIdxList);
        if sum(saveJunctIdx)~=0
            veilStem(iFrame).saveJunct = 1;
            toSaveIdx = vertcat(juncCC.PixelIdxList{saveJunctIdx}) ;
            junctionMask(toSaveIdx) = 0;
        end
    end
    %%
    veilStem(iFrame).backbone = backbone;
    [EPs,~,coords,branchPtMask] = skel2graph2D(notBody|newBodyMask);
    
    % make sure the endpoint is not the neurite entrance...
    indEP = sub2ind(size(img),EPs(:,1),EPs(:,2));
    
    %idxSave = find(indEP == idxEnterNeurite); %% note sometimes bug here... should make so reiterate if this fails...
    %instead of doing find use intersect
    overlap = intersect(idxEnterNeurite,indEP);
    % save the entering neurite pieces from erosion
    if ~isempty(overlap)
        idxSave = arrayfun(@(x) find(indEP == overlap(x)),1:length(overlap));
        idxSave = idxSave';
        
        EPs(idxSave,:) = [];
        coords(idxSave) = [];
    end
    
    if ~isempty(vertcat(coords{:}))
        while ~isempty(EPs)
            
            if ~isempty(vertcat(coords{:}))
                coordBBOver = vertcat(coords{:});
                idxBBOver = sub2ind(size(img),coordBBOver(:,1),coordBBOver(:,2));
                idxEP = sub2ind(size(img),EPs(:,1),EPs(:,2));
                notBody([idxBBOver ; idxEP]) = 0;
                backbone([idxBBOver ;idxEP]) = 0;
                backbone(branchPtMask==1) = 1;
                backbone = bwmorph(backbone,'thin',inf);
                
                notBody(branchPtMask==1) = 1;
                notBody = bwmorph(notBody,'thin',inf);
            end
            new = (notBody|newBodyMask);
            new = logical(getLargestCC(new));
            
            [EPs,~,coords,branchPtMask] = skel2graph2D(new);
            
            indEP = sub2ind(size(img),EPs(:,1),EPs(:,2));
            
            overlap = intersect(idxEnterNeurite,indEP);
            % save the entering neurite pieces from erosion
            if ~isempty(overlap)
                idxSave = arrayfun(@(x) find(indEP == overlap(x)),1:length(overlap));
                idxSave = idxSave';
                
                EPs(idxSave,:) = [];
                coords(idxSave) = [];
            end
            
        end % while
    end % isempty
    
    newBodyMask= imfill(newBodyMask,'holes');
    
    %% This was another means for removing end piece by checking labels:
    % Check
    backbone2Dil = backbone;
    backbone2Dil(backbone==1 & newBodyMask==1) = 0;
    
    CCNotBody = bwconncomp(backbone2Dil);
    
    CCBody = bwconncomp(newBodyMask);
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
            backbone2Dil(CCNotBody.PixelIdxList{iCC})= 0;
        end
    end % for iCC
    %% TSMovie
    if ip.Results.TSMovie == true
        close gcf
        setFigure(nxLarge,nyLarge,'off')
        imshow(-imgLarge,[]);
        hold on
        roiYXNB = bwboundaries(newBodyMask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','b','linewidth',2),roiYXNB);
        
        idxDil =  find(backbone2Dil==1);
        [yBackD,xBackD] =  ind2sub([ny,nx],idxDil);
        scatter(xBackD,yBackD,10,c(2,:),'filled');
        %         pixels = round(10/pixSizeMic);
        %         plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
        for i = 1:2
            if ip.Results.screen2png
                helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            else
                saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            end
            countMovie = countMovie+1;
        end
        close gcf
    end % TSMovie
    %%  RE-DILATE the ridge regions
    %dilBB = imdilate(backbone2Dil,strel('disk',4)); % arbitrary... need to find small function redilate based on ridge estimation
    veilStemNodeMask = newBodyMask;
    backboneInfoC = backboneInfo(iFrame);
    [fullMask,cycleFlag,TSFigs2] = gcaResolveVeilStemCycles(backbone2Dil,veilStemNodeMask,backboneInfoC,img,p);
    %% Troubleshoot Overlays: Veil/Stem Cycles
    if ip.Results.TSOverlays
        TSFigs = [TSFigs1 TSFigs2];
        % Save any trouble shoot figures associated with the cycle resolutions
        % make the directories for the figures if required.
        for iFig = 1:length(TSFigs)
            if ~isdir([ip.Results.OutputDirectory filesep TSFigs(iFig).name])
                mkdir([ip.Results.OutputDirectory filesep TSFigs(iFig).name])
            end
        end
        type{1} = '.fig';
        type{2} = '.tif';
        
        if ~isempty(TSFigs)
            if ip.Results.screen2png
                arrayfun(@(x) helperScreen2png([ip.Results.OutputDirectory filesep x.name filesep ...
                    num2str(iFrame,'%02d') '.png'],'figureHandle',x.h),TSFigs);
            else
                
                for iType = 1:numel(type)
                    arrayfun(@(x) saveas(x.h,...
                        [ip.Results.OutputDirectory filesep x.name filesep num2str(iFrame,'%02d') type{iType}]),TSFigs);
                end
            end
        end
        
        close all
        clear TSFigs TSFigs1 TSFigs2
    end  % TS Overlays
    %% TS Movie
    if ip.Results.TSMovie == true
        close gcf
        setFigure(nxLarge,nyLarge,'off');
        imshow(-imgLarge,[]);
        hold on
        roiYXNB = bwboundaries(fullMask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color',[ 0.0039  ,  0.4264 ,   0.3848],'LineWidth',2),roiYXNB);
        % pixels = round(10/pixSizeMic);
        % plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
        if ip.Results.writeTitles
            text(5,5,'Veil/Stem Mask Final','color',[ 0.0039  ,  0.4264 ,   0.3848],'Fontsize',10);
        end
        for i = 1:2
            if ip.Results.screen2png
                helperScreen2png([saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            else
                saveas(gcf,[saveDirMov filesep num2str(countMovie,'%03d') '.png']);
            end
            
            countMovie = countMovie+1;
        end
    end % TSMovie
    %% Save the veilStem Information
    veilStem(iFrame).finalMask = fullMask;
    veilStem(iFrame).idxEnterNeurite = idxEnterNeurite;
    veilStem(iFrame).cycleFlag = cycleFlag;
    veilStem(iFrame).timeStamp = clock;
    
    if ip.Results.getGITHashTag
        hashTag =  gcaArchiveGetGitHashTag;
    else
        hashTag = NaN;
    end
    veilStem(iFrame).hashTag = hashTag;
    
    paramsArchived(iFrame) = pToSave;
    save([ip.Results.OutputDirectory filesep 'params.mat'],'paramsArchived');
    save([ip.Results.OutputDirectory filesep 'veilStem.mat'],'veilStem','-v7.3');
    
end % for iFrame
end

