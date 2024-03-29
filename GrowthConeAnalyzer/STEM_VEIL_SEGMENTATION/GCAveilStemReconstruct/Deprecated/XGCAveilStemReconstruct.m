function [ analInfo,hTroubleshoot ] = GCAveilStemReconstruct(img,maskForErod,backboneInfo,paramsIn,iFrame)
% GCAveilStemReconstruct

% INPUT:
%   img: a RxC array of the image where R and C correspond to the image size
%   maskForErod: a RxC array of the binary mask of the corresponding image where R and C correspond to the
%   image/mask size
%   backboneInfo:
%   paramsIn: a structure of input parameters including fields
%
% OUTPUT
% analInfo: A structure with perhaps a bit too much info- within it is the
% final masks


%% Check Parameters
cycleFlag = 0; % initiate cycle flag
% erod out filopodia like structures and ridges

% imopen(maskForErod,strel('disk',p.DiskSizeForErod,0));%

% performs 'smarter' opening that preserves spanning pieces of thin
% veil/stem
[erodForBody,saveMask]  =   findBreakageCCsFromErosion(maskForErod,6);

%%%%% Maybe delete
yxvalues = bwboundaries(erodForBody);
% might not keep this in the output
analInfo(iFrame).bodyEst.erodForBody = yxvalues;
%thinnedBodyAll = bwmorph(erodForBody,'thin','inf');
thinnedBodyAll = bwmorph(erodForBody,'skel','inf');
analInfo(iFrame).bodyEst.thinnedBodyAll = thinnedBodyAll;
clear yxvalues
%%%%%

%% Load Backbone
backbone = backboneInfo.backboneSeedMask;
%% 20140822 : might in the future need to change this option but for
% now. delete all the body pieces that are located along the edge of
% the frame that are NOT overlapping with the main entrance BB.

CCBodyPieces = bwconncomp(erodForBody);
[imSize] = size(img);
ny = imSize(1);
nx = imSize(2);
boundaryMask = zeros(imSize);
boundaryMask(1:ny,1) =1;
boundaryMask(1:ny,nx)=1;
boundaryMask(1,1:nx)= 1;
boundaryMask(ny,1:nx) =1;
pixBound = find(boundaryMask==1); % get pixel indices

% find pieces that overlap with boundary but do not overlap with
% the backbone seed
pixBB = find(backbone==1);

floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCBodyPieces.PixelIdxList);
boundIdxCC = cellfun(@(x) ~isempty(intersect(pixBound,x)),CCBodyPieces.PixelIdxList);
toKeep = ~(floatingBodyIdxCC & boundIdxCC);
% save the pieces of boundIdxCC for later just in case there is
% an error where the bacbone is completely fucked and you take out good portions
boundaryBodyPixelIdx = CCBodyPieces.PixelIdxList(boundIdxCC & floatingBodyIdxCC);
erodForBody = zeros([ny,nx]);
erodForBody(vertcat(CCBodyPieces.PixelIdxList{toKeep}))=1;   % a little redundant here with below clean up later

%%
xyEnterNeurite = backboneInfo.coordsEnterNeurite;
%     maxNMSLarge = analInfo(iFrame).bodyEst.maxNMSLarge; % note should maybe put these all into
%     cutoff = analInfo(iFrame).bodyEst.cutoff; %

idxEnterNeurite = sub2ind(size(img),xyEnterNeurite(:,2), xyEnterNeurite(:,1));
cleanedRidge = backboneInfo.bodyReconstruct.AfterConnect;
cleanedRidge(backbone==1) = 0;
% added 20140307
cleanedRidge =  bwmorph(cleanedRidge,'spur');
%% SIDE NOTE STEP: Make some movies of the iterative reconstruct if user desires
if paramsIn.makeMovie == 1;
    imgLarge = [ img ;img];
    
    [nyLarge,nxLarge] = size(imgLarge);
    setFigure(nxLarge,nyLarge,'off');
    idx = find(backboneInfo().backboneSeedMask==1);
    [yback,xback] = ind2sub([ny,nx],idx);
    
    % start by plotting the backbone coord and then showing
    % iterations
    imshow(-imgLarge,[]);
    hold on
    
    scatter(xback,yback,10,'r','filled');
    hold on
    scatter(xyEnterNeurite(:,1),xyEnterNeurite(:,2),50,'y','filled');
    hold on
    saveDirMov =  [paramsIn.OutputDirectory filesep 'ExampleMovie'];
    if ~isdir(saveDirMov)
        mkdir(saveDirMov)
    end
    pixSizeMic = 0.216;
    pixels = round(10/pixSizeMic);
    plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
    
    saveas(gcf,[saveDirMov filesep '01.png']);    % first part of iteration
end

%% MAIN STEP: USE THE BACKBONE INFO TO CONNECT EROSIONS
% INITIATE
stopFlagBodyReconstruct =0;
% Initiate the reconstruction of the body
erodfilo = erodForBody;
iter = 1;
% get the connected components corresponding to "fat body"
CCAllBody = bwconncomp(erodfilo); %

% START WHILE LOOP: Keep adding thick veil pieces until no longer find viable candidates
while stopFlagBodyReconstruct ==0
    pixBB = find(backbone==1);
    % find the pieces that do not intersect with the putative backbone
    floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
    %%
    % Make the new mask of only ridge connected body parts (ie take out
    % floating parts)
    CC = CCAllBody; % copy cc structure and truncate to make the new body mask
    CC.NumObjects = CCAllBody.NumObjects -length(CCAllBody.PixelIdxList(floatingBodyIdxCC));
    CC.PixelIdxList(floatingBodyIdxCC) = [];
    newBodyMask= labelmatrix(CC);
    if paramsIn.makeMovie == 1;
        roiYXBN =  bwboundaries(newBodyMask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','y','linewidth',2),roiYXBN);
        saveas(gcf,[saveDirMov filesep '02.png']);
    end
    %% Here we (non-elegantly) use some previous temporal information to find the correct path of the neurite
    % if there is evidence that we have a severely truncated
    % segmentation for a given frame - the algorithm currently
    % progresses forward - so the first frame does NOT benefit from the
    % temporal information
    if iter == 1 && sum(newBodyMask(:))==0 % you didn't add anything backbone might just be truncated
        % check to see if frame 1
        if iFrame ~=1
            %try last frame
            oldBody = analInfo(iFrame-1).masks.thickBodyMask;
            pixOldBody= find(oldBody==1);
            % try again
            % first test to make sure erodfilo is not equal to zero -
            % part of the problem is you might have had a very good
            % path estimate from the erosion and a very bad backbone
            % estimate this good path with connect at the boundary and
            % would be removed above this way check to bring it back
            % here
            if sum(erodfilo(:))== 0;% this is likely the problem you have nothing there!
                erodfilo(vertcat(boundaryBodyPixelIdx{:}))=1;  % either that or you can not take out if problem (this might be why it would be good to set up as a cycle problem in end)
                
            end
            CC = bwconncomp(erodfilo);
            floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixOldBody,x)),CC.PixelIdxList);
            
            CC.NumObjects = CC.NumObjects -length(CC.PixelIdxList(floatingBodyIdxCC));
            CC.PixelIdxList(floatingBodyIdxCC) = [];
            newBodyMask= labelmatrix(CC);
            
            
            
            backbone = backboneInfo(iFrame-1).backboneSeedMask;
            % make sure the idxEnterNeurite is updated such that you know you
            % are using the previous for saving the junction s
            
            backboneInfo(iFrame).backboneSeedMask = backbone;  % just in case you need to use it again.
            % maybe flag % added 20140321
            
            xyEnterNeurite = backboneInfo(iFrame-1).coordsEnterNeurite;
            %     maxNMSLarge = analInfo(iFrame).bodyEst.maxNMSLarge; % note should maybe put these all into
            %     cutoff = analInfo(iFrame).bodyEst.cutoff; %
            %
            idxEnterNeurite = sub2ind(size(img),xyEnterNeurite(:,2), xyEnterNeurite(:,1));
            backboneInfo(iFrame).coordsEnterNeurite = xyEnterNeurite;
            cleanedRidge = backboneInfo(iFrame-1).bodyReconstruct.AfterConnect;
            backboneInfo(iFrame).bodyReconstruct.AfterConnect = cleanedRidge;
            cleanedRidge(backbone==1) = 0;
            % added 20140307
            cleanedRidge =  bwmorph(cleanedRidge,'spur');
        end
        
        
    end
    newBodyMask(newBodyMask>0) =1;
    
    % first test if need to break
    
    if sum(floatingBodyIdxCC)==0; % all potential body pieces were assigned to neurite
        % ADDED 20141009 - previously we didn't want to introduce
        % competing paths because we had no mechanism to fix them.
        % now we would like to keep these competing paths just in
        % EVEN in the case when the estimated backbone picks up
        % ALL body pieces- as there are several cases where this
        % is not necessarily the optimal path
        pixBodySave = find(newBodyMask==1);
        CCRidgeBone = bwconncomp(cleanedRidge);
        % stopFlagBodyReconstruct = 1 ; % stop Fat Body Reconstruct proceed to backbone erosion
        idxNoIntersectRidgeCC = cellfun(@(x) isempty(intersect(pixBodySave,x)),CCRidgeBone.PixelIdxList);
        
        CCRidgeBone.PixelIdxList(idxNoIntersectRidgeCC) = [];
        CCRidgeBone.NumObjects = CCRidgeBone.NumObjects - sum(idxNoIntersectRidgeCC);
        backbone(vertcat(CCRidgeBone.PixelIdxList{:}))=1;
        analInfo(iFrame).bodyEst.rmHiIntPieces = [];
        break
        
        
        
        
    end % if sum % note eventually might want to look at scale information
    
    
    %% New Step: try to bridge small gaps in seed with with other nearby ridges/BodyCCs
    % Find end points of pixBB. use these as seed points
    % Find all points either ridge or floating body with in 2 pixels
    % of the endpoint. do this from the new body mask
    
    
    
    pixBodySave = find(newBodyMask==1);
    CCRidgeBone = bwconncomp(cleanedRidge);
    
    % Get rid of all ridges that fail to connect to the new body
    % parts
    idxNoIntersectRidgeCC = cellfun(@(x) isempty(intersect(pixBodySave,x)),CCRidgeBone.PixelIdxList);
    CCRidgeBoneConnect = CCRidgeBone;
    CCRidgeBoneConnect.PixelIdxList(idxNoIntersectRidgeCC) = [];
    CCRidgeBoneConnect.NumObjects = CCRidgeBoneConnect.NumObjects - sum(idxNoIntersectRidgeCC);
    % get endpoints of these ridges that connect to the body.
    CCCleanedRidge = CCRidgeBone;
    CCCleanedRidge.PixelIdxList(~idxNoIntersectRidgeCC) =[];
    CCCleanedRidge.NumObjects = CCRidgeBone.NumObjects -sum(~idxNoIntersectRidgeCC);
    
    pixRidgeConn = vertcat(CCRidgeBoneConnect.PixelIdxList{:});
    % backbone should now be all pixels that connect to the body...
    backbone(pixRidgeConn) = 1;
    labelMatCandRidge = labelmatrix(CCCleanedRidge);
    cleanedRidge = labelMatCandRidge>0;
    
    if paramsIn.makeMovie == 1
        [yBack,xBack]=  ind2sub([ny,nx],pixRidgeConn);
        scatter(xBack,yBack,10,'r','filled') % 'color','r','linewidth',2);
        saveas(gcf,[saveDirMov filesep '03.png']);
    end
    
    % pixBB = pixRidgeConn; % added 2014-01-29
    if sum(cleanedRidge(:))~=0
        EPCoordsBodyMask = getEndpoints([find(newBodyMask==1) ; pixBB ; pixRidgeConn],[ny,nx]);
        % get coords of all ridges and body pieces
        % ridges
        
        
        
        % for now only try is isempty
        if isempty(EPCoordsBodyMask)
            % try the endpoint of the backbone.
            
            
            
            EPCoordsBodyMask = getEndpoints(pixBB,[ny,nx]);
            
            idxEPs = sub2ind([ny,nx],EPCoordsBodyMask(:,2),EPCoordsBodyMask(:,1));
            
            
            coordsIn = backboneInfo(iFrame).coordsEnterNeurite;
            idxIn = sub2ind([ny,nx],coordsIn(:,2),coordsIn(:,1));
            EPCoordsBodyMask(idxEPs==idxIn,:)= [];
            
        end
        
        idxPixRidgeCand = find(cleanedRidge==1);
        [xyCandRidge(:,2), xyCandRidge(:,1)] = ind2sub([ny,nx],idxPixRidgeCand); % for some reason when put into a mat first output is x? opposite if just a single var
        
        %[yFloatBody,xFloatBody] = ind2sub([ny,nx],floatingBodyPix);
        %?allCand = [xyCandRidge;roiXYFloaters]; % eventually I would like to make this all candidates
        % first test for ridge connections
        [idxRidge,dRidge]  = KDTreeBallQuery(xyCandRidge,EPCoordsBodyMask,5); % for now just do ridges and link by
        
        if ~isempty(vertcat(idxRidge{:})); % if it found any point within range
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
            % label larger than the nendpoint.
            NNodeQuery = size(EPCoordsBodyMask,1);  % get the
            nodeLabels = E(:,1);
            [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two separate EPs are competing over the same body point
            
            % reason note: just in case two filo are competing over the same seed point
            nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
            EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
            numberNodes = length(inputLinks) + NNodeQuery;
            
            %Collect all distance values
            %D = arrayfun(@(i) vertcat(distNode{i}{:}),1:numel(distNode),'uniformoutput',0);
            
            D = max(dRidge)-dRidge;
            % normalize
            D= D./max(D);
            
            costTotal = (D); % eventially make this path dependent
            %
            % set up for graph match edge number, idx of path, and associated cost
            
            E = [E costTotal D ];
            
            if ~isempty(EFinal)
                %     idx = E(:,1) < E(:,2);
                
                %     E = E(idx,:); % remove redundancy
                
                M = maxWeightedMatching(numberNodes, EFinal, costTotal);
                % check for double labels
                % convertBack
                % E = [candFiloNodes(nodeLabels) inputLinks(nodeLabelsInput)]; % convert back to original indices of input and query points
                E = E(M,:);% get those edges that matched (from original indexing)
            end   % ise
            if ~isempty(E) % might be empty now if all were repeats
                
                % add linear segments corresponding to linked endpoints
                % actually this gets me into trouble if really want to link
                % these effectively need to consider the fluorescence intensity
                % into the cost of linking each point this can throw off my fits
                % it's a little stupid because we have some of this junction info
                % before the NMS but I throw it away
                %    goodConnect = iSeg(M);
                
                % convert to pixIdx
                %            pixGoodConnect = cellfun(@(i) sub2ind(dims,i(:,2),i(:,1)), goodConnect,'uniformoutput',0);
                %            out(vertcat(pixGoodConnect{:}))= 1;
                
                % TO DELETE already calculated for all b/c need intensities so
                % redundant
                for iMatch = 1:size(E,1)
                    % get coords
                    xCoordQuery = EPCoordsBodyMask(E(iMatch,1),1);
                    yCoordQuery =EPCoordsBodyMask (E(iMatch,1),2) ;
                    xCoordInput = xyCandRidge(E(iMatch,2),1);
                    yCoordInput = xyCandRidge(E(iMatch,2),2);
                    %iseg = gcaBresenham([xe(E(i,1)) ye(E(i,1))], [xe(E(i,2)) ye(E(i,2))]);
                    iseg = gcaBresenham([xCoordQuery yCoordQuery],...
                        [xCoordInput yCoordInput]);
                    % save the coords coorsponding to the pixels you are adding
                    % need this to potentially tell not to wreck the junction
                    % mask
                    segSave{iMatch}= sub2ind([ny,nx], iseg(:,2), iseg(:,1));
                    backbone(sub2ind([ny,nx], iseg(:,2), iseg(:,1))) = 1;% add to mask
                    labels2Keep = labelMatCandRidge(yCoordInput,xCoordInput);
                    backbone(vertcat(CCCleanedRidge.PixelIdxList{labels2Keep})) = 1;
                    cleanedRidge(vertcat(CCCleanedRidge.PixelIdxList{labels2Keep})) = 0; % get ridge of in the cleaned ridge cand list
                end
                
                % in some cases just make sure to thin and it will fix junction
                % problems
                backbone = bwmorph(backbone,'thin','inf');
                % update pixBB
                pixBB = find(backbone==1) ;
                % does the new piece of backbone help pick up any more body
                % pieces if it does save junction
                
                % get the logical indices of those body pieces that intersect
                % with the backbone
                floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
                % if
                if ~isempty(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCC}))
                    
                    saveJunct{iter}= vertcat(segSave{:});
                end
                newBodyMask(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCC}))= 1;
                if paramsIn.makeMovie == 1
                    roiYXNB = bwboundaries(newBodyMask);
                    cellfun(@(x) plot(x(:,2),x(:,1),'color','y','linewidth',2),roiYXNB);
                    saveas(gcf,[saveDirMov filesep '04.png']);
                end
                
                
                
            end % isempty(E)
            
            
            
        end % isempty(vertcat(idxRidge{:}))
        
        clear xyCandRidge  segSave
        % New backbone complete
        %%
    end % if sumCleanedRidge
    
    % in some cases just make sure to thin and it will fix junction
    % problems
    backbone = bwmorph(backbone,'thin','inf');
    % update pixBB
    pixBB = find(backbone==1) ;
    % does the new piece of backbone help pick up any more body
    % pieces if it does save junction
    
    % get the logical indices of those body pieces that intersect
    % with the backbone
    floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
    % if
    %            if ~isempty(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCC}))
    %
    %             saveJunct{iter}= vertcat(segSave{:});
    %            end
    newBodyMask(vertcat(CCAllBody.PixelIdxList{~floatingBodyIdxCC}))= 1;
    
    
    
    
    % get the logical indices of those body pieces that intersect
    % with the backbone
    floatingBodyIdxCC = cellfun(@(x) isempty(intersect(pixBB,x)),CCAllBody.PixelIdxList);
    
    
    if sum(floatingBodyIdxCC)==0; % all potential body pieces were assigned to neurite
        % stopFlagBodyReconstruct = 1 ; % stop Fat Body Reconstruct proceed to backbone erosion
        analInfo(iFrame).bodyEst.rmHiIntPieces = [];
        break
        
        
    end % if sum
    
    
    % if all the body pieces were not used iteratively test if the
    % ridges eminating from the body pieces
    
    % take new body mask and search for ridges
    % get pixels of the new body mask
    pixBodySave = find(newBodyMask==1);
    CCRidgeBone = bwconncomp(cleanedRidge);
    
    % Get rid of all ridges that fail to connect to the new body
    % parts
    idxNoIntersectRidgeCC = cellfun(@(x) isempty(intersect(pixBodySave,x)),CCRidgeBone.PixelIdxList);
    CCRidgeBoneConnect = CCRidgeBone;
    CCRidgeBoneConnect.PixelIdxList(idxNoIntersectRidgeCC) = [];
    CCRidgeBoneConnect.NumObjects = CCRidgeBoneConnect.NumObjects - sum(idxNoIntersectRidgeCC);
    
    
    
    
    
    
    % test again for overlap of these high confidence ridges with any
    % misc body parts
    CCBody = bwconncomp(erodForBody);% note might have to change this here.... i think this implementation is redundant
    
    testBodyCC = CCBody.PixelIdxList(floatingBodyIdxCC);
    testBodyPixels = vertcat(CCBody.PixelIdxList{floatingBodyIdxCC});
    connectedRidgePix = vertcat(CCRidgeBoneConnect.PixelIdxList{:});
    idxIntersectBody2 = cellfun(@(x) ~isempty(intersect(connectedRidgePix,x)),testBodyCC); % find which body pieces (if any) overlap with ridge
    idxRidgesFatBodyConnect = cellfun(@(x) ~isempty(intersect(testBodyPixels,x)), CCRidgeBoneConnect.PixelIdxList);
    
    newToKeep = testBodyCC(idxIntersectBody2);
    newBodyMask(vertcat(newToKeep{:})) = 1; % add pixels to mask
    backbone(vertcat(CCRidgeBoneConnect.PixelIdxList{idxRidgesFatBodyConnect})) =1; % add pixles to backbone
    % If all body parts have been confirmed OR could not assign any of the remaining based
    % on geometry/ridge path arguments then stop the iterations and record any floating body pieces
    % remaining- if the frame is marked a neurite length outlier- these pieces will be reconsidered for reattachment.
    % if piece attached and one still remaining
    % re-iterate.
    if (sum(idxIntersectBody2)/length(idxIntersectBody2) ==1 || sum(idxIntersectBody2) == 0 )% all body parts have been assigned or couldn't assign any. 08-03-2013 this is a weird stop flag? check again
        %
        if sum(idxIntersectBody2)==0
            
            % save as a mask- though might want to reconsider what is the most
            % effective method of saving
            rmHiIntPieces = zeros(ny,nx);
            rmHiIntPieces(vertcat(testBodyCC{:}))=1;
            
            analInfo(iFrame).bodyEst.rmHiIntPieces = rmHiIntPieces;
            
        else
            
            analInfo(iFrame).bodyEst.rmHiIntPieces = [];
        end
        
        
        
        
        stopFlagBodyReconstruct =1;
        
    end % if sum
    iter = iter +1;
    if iter>100
        error('GCAAnalyzer:NeuriteBodyEst:FailureToTerminate','Too Many Iterations Check Code');
    end
    
    
    
end % while stopFlagBodyReconstruct (don't add any more body parts)


% get the fat body masks with inner pixels removed for backbone
% body erosion
newBodyXY = bwboundaries(newBodyMask);
idxBody = cellfun(@(x) sub2ind(size(newBodyMask),x(:,1),x(:,2)),newBodyXY,'uniformoutput',0);
idxBody = vertcat(idxBody{:});
bodyNoFill = zeros(size(img));
bodyNoFill(idxBody) = 1;

% prune any junctions in the backbone
%          backbonePreErosion = backbone;
%           backbonePreErosion= bwmorph(backbonePreErosion,'thin','inf');
% nn = padarrayXT(double(backbonePreErosion~=0), [1 1]);
% sumKernel = [1 1 1];
% nn = conv2(sumKernel, sumKernel', nn, 'valid');
% nn1 = (nn-1) .* (backbonePreErosion~=0);
% junctionMask = nn1>2;
% backbonePreErosion(junctionMask) =0;
%  backbonePreErosion = backbonePreErosion | newBodyMask;
% backbonePreErosion = getLargestCC(backbonePreErosion);
% backbonePreErosion(newBodyMask==1)= 0;
% backbone = backbonePreErosion; % let's give it a try



% prune the backbone at the ends (only keep the intersections)
% break the junctions in notBody
%       nn = padarrayXT(double(notBody~=0), [1 1]);
%   sumKernel = [1 1 1];
%   nn = conv2(sumKernel, sumKernel', nn, 'valid');
%   nn1 = (nn-1) .* (notBody~=0);
%   junctionMask = nn1>2;
%   notBody(junctionMask) =0;
%   backbone(junctionMask) = 0 ;
notBody = double(backbone).*~double(newBodyMask);

% 12-11 note just make sure to thin it here. check why this thin step was
% miseed before ...though I previously implemented it...

notBody = bwmorph(notBody,'thin','inf');

%% 2013_12_08 check this - in theory should not have junctions if did this cleanly (rethink)
% Break junctions
% also don't want to break junctions in the original backbone...
nn = padarrayXT(double(notBody~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
forJunctBreak = notBody~=0;
% for now just make sure not f-ing up the original backbone
%
forJunctBreak(backboneInfo(iFrame).backboneSeedMask==1) = 0;
%quick fix test for junctions in original cleaned ridge
% add to save

cRidgeTest = backboneInfo(iFrame).bodyReconstruct.AfterConnect;
nnR = padarrayXT(double(cRidgeTest~=0), [1 1]);
nnR= conv2(sumKernel, sumKernel', nnR, 'valid');
nnR1 = (nnR-1).*cRidgeTest;
junctSaveMask = (nnR1>2);
CCTest = bwconncomp(junctSaveMask);
if CCTest.NumObjects~=0
    for i = 1:numel(CCTest.PixelIdxList)
        pixIdx = CCTest.PixelIdxList{i};
        if exist('saveJunct','var') ==1
            
            saveJunct{end+1} = pixIdx;
        else
            
            saveJunct{1} = pixIdx;
        end
    end
    
    
    
    
end


forJunctBreak(backboneInfo(iFrame).bodyReconstruct.AfterConnect) = 0; %
% remove anything that was originall in teh cleaned version (eventually
% need to make it so that you are absolutely certain that there is no
% residual junctions left in the cleaned ridge connect.
% 2014 03 09


nn1 = (nn-1) .* forJunctBreak;
junctionMask = nn1>2;
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
        toSaveIdx = vertcat(juncCC.PixelIdxList{saveJunctIdx}) ;
        
        
        junctionMask(toSaveIdx) = 0;
    end
end
%%

% junctionMask(vertcat(segSave{:})) = 0; % don't break junctions you made previous step
notBody(junctionMask) =0;
backboneTest = notBody;
backboneTest(junctionMask) = 0 ;



% %%i don't know what I was doing here...
test = backboneTest | newBodyMask;
CCTest = bwconncomp(test); % try to test this for the piece that is overlapping with the  body
idxBodyOverlap = cellfun(@(x) isempty(intersect(x,find(newBodyMask==1))),CCTest.PixelIdxList);
backbone(vertcat(CCTest.PixelIdxList{idxBodyOverlap==1})) = 0; % get rid of floaters
%
% one way to deal with it is to simpy test for overlap with junctions you
% already created

test = notBody|bodyNoFill;% get the backbone pieces and the no fill body for erosion
test2 = bwmorph(test,'thin',inf);
backbone(test & ~test2) = 0 ; % not  sure what I was doing here
%get rid of singletons : this might not be necessary if don't break
%junctions in the first step
CCTest2 = bwconncomp(test2);
csize = cellfun(@(x) length(x),CCTest2.PixelIdxList);
backbone(vertcat(CCTest2.PixelIdxList{csize==1})) = 0; % set these = to zero
CCTest2.PixelIdxList(csize==1) = [];
CCTest2.NumObjects = CCTest2.NumObjects - sum(csize==1);
test2= labelmatrix(CCTest2);
test2 = test2>0;
[EPs,~,coords] = skel2graph2D(test2);


%%
analInfo.bodyEst.backbone = backbone;



if ~isempty(vertcat(coords{:}))
    indEP = sub2ind(size(img),EPs(:,1),EPs(:,2));
    
    
    
    
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
        idxBBOver = sub2ind(size(img),coordBBOver(:,1),coordBBOver(:,2));
        idxEP = sub2ind(size(img),EPs(:,1),EPs(:,2));
        backbone([idxBBOver ; idxEP]) = 0;
    end
end

%%% take the new body mask and make some calculation corresponding
%%% to geometry initially we discussed we wanted this to be skel
%%% then the avg dist transformation of that.
%thinnedBody = bwmorph(newBodyMask,'thin','inf');
%analInfo(iFrame).bodyEst.skelFatBody = thinnedBody;
newBodyMask= imfill(newBodyMask,'holes');
%         distTrans = bwdist(~newBodyMask);
%         distTrans = distTrans.*thinnedBody;
%         valDist = distTrans(distTrans~=0);
%geoParamFatBody = mean(valDist).*0.216; % make sure to convert to microns
%analInfo(iFrame).bodyEst.geoParamFatBody = geoParamFatBody;
%area = bwarea(newBodyMask).*0.216*0.216;
%analInfo(iFrame).bodyEst.area = area;
thickBodyMask = newBodyMask; % save here in case you want to put into the protrusion software.


%%
backbone2Dil = backbone;
backbone2Dil(backbone==1 & newBodyMask==1) = 0;
% for now the dilation might be the best


%% before redilate the backbone find notBody pieces that overlap the body twice...
%(likely a very stupid way of doing this but can make it work for now)
% dilate the notBody part and see if two CCs

% give them a label
% get rid of all but endpoints
% dilate endpoints
% test for overlap with labelednewbodyMask
% if overlap only == 1 type same body piece discard, but save idxEnter

% goal is to see if the label of the CC body piece and the body blob are
% the same or different
% if they are the same it means the CC body piece is spanning one body
%% quick fix made 12-08-13
% backbone2Dil = bwmorph(backbone2Dil,'bridge');
% backbone2Dil = bwmorph(backbone2Dil,'thin',inf);


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
    
end
%%
if paramsIn.makeMovie ==1
    close gcf
    setFigure(nxLarge,nyLarge,'off')
    imshow(-imgLarge,[]);
    hold on
    roiYXNB = bwboundaries(newBodyMask);
    cellfun(@(x) plot(x(:,2),x(:,1),'color','y','linewidth',2),roiYXNB);
    
    idxDil =  find(backbone2Dil==1);
    [yBackD,xBackD] =  ind2sub([ny,nx],idxDil);
    scatter(xBackD,yBackD,10,'r','filled');
    pixels = round(10/pixSizeMic);
    plotScaleBar(pixels,pixels/20,'Label','10um','Color',[0 0 0]);
    saveas(gcf,[saveDirMov filesep '05.png']);
end


%%% FOR NOW JUST ARBITRARILY RE-DILATE
dilBB = imdilate(backbone2Dil,strel('disk',4));

roiYXthinBody= bwboundaries(dilBB);
roiYX1 = vertcat(roiYXthinBody{:});
roiYXthickBody = bwboundaries(newBodyMask);
roiYX2 = vertcat(roiYXthickBody{:});

if ~isempty(roiYX1)
    
    pixIndThinBody = sub2ind(size(img),roiYX1(:,1),roiYX1(:,2));
else
    pixIndThinBody = [];
end

if ~isempty(roiYX2)
    pixIndThickBody= sub2ind(size(img),roiYX2(:,1),roiYX2(:,2));
else
    pixIndThickBody = [];
end


fullMask = dilBB | newBodyMask;

% take largest cc and fill holes
fullMask = logical(getLargestCC(fullMask));
%% TEST FOR CYCLES AND CORRECT
prefill = fullMask; % added 20140819
fullMask = imfill(fullMask,'holes');
if ~isequal(prefill,fullMask);% you have cycle.
    cycleFlag = 1;
    display('you have a cycle');
    % deconstruct the body as a graph
    % label the new body mask
    labelsBody = bwlabel(newBodyMask);
    % get body nodes
    nodeNum = unique(labelsBody(labelsBody~=0));
    % Problem 20141009 becomes that dilation of 4 is can cross the
    % small body pieces resulting in a merging of the two paths
    % quick solution is use a smaller dilation for the label making
    % In the end we want to get a better estimate for the redilation
    % of the paths anywa
    %dilBBForLabels =  imdilate(backbone2Dil,strel('disk',2));
    
    
    
    %% Small test to mask sure dilation does not merge edge paths
    % the idea here is want to dilate for the intensity
    % integration reponse metrics but don't want this to be -
    % think about reworks in this coding for the final release as
    % it is a bit rough.
    
    CCPreDil = bwconncomp(backbone2Dil);
    CCEdges = bwconncomp(dilBB);
    stopFlagLowerDil = CCPreDil.NumObjects  > CCEdges.NumObjects;
    countDilDec = 1;
    while stopFlagLowerDil >0
        dilBB = imdilate(backbone2Dil,strel('disk',4-countDilDec));
        CCEdges = bwconncomp(dilBB);
        stopFlagLowerDil = CCPreDil.NumObjects  > CCEdges.NumObjects;
        countDilDec = 1 + countDilDec;
    end % while
    % END TEST 1
    %% 2nd Test for problems in the case dilation was too large.
    labelsC = bwlabel(dilBB);
    % for each label get the pixels that overlap body labels%
    CCEdges = bwconncomp(labelsC);
    edges = cellfun(@(x) unique(labelsBody(x)),CCEdges.PixelIdxList,'uniformoutput',0);
    % dilate each piece and give
    %conForLabel = imdilate(dilBB,strel('disk',3));
    edges =  cellfun(@(x) x(x~=0)',edges,'uniformoutput',0);
    
    % test for problems in the dilation
    numVertices = cellfun(@(x) length(x) ,edges);
    stopFlagLowerDil = sum(numVertices>2);% initiate stop flag
    % Also check to make sure that CC before dilation remains
    % consistent
    
    countDilDec = 1; % initiate count
    %if sum(numVertices > 2) ~=0 % test for problem cases where the dilation of 4 was too large
    % and spanned the body the small node...
    while stopFlagLowerDil >0
        dilBB = imdilate(backbone2Dil,strel('disk',4-countDilDec));
        labelsC = bwlabel(dilBB);
        CCEdges = bwconncomp(labelsC);
        edges = cellfun(@(x) unique(labelsBody(x)),CCEdges.PixelIdxList,'uniformoutput',0);
        % dilate each piece and give
        %conForLabel = imdilate(dilBB,strel('disk',3));
        edges =  cellfun(@(x) x(x~=0)',edges,'uniformoutput',0);
        
        
        numVertices = cellfun(@(x) length(x),edges);
        stopFlagLowerDil = sum(numVertices>2);
        countDilDec = 1 + countDilDec;
    end % while
    
    
    
    
    %% start calculating scores for the edges
    % likely candidates have high fluorescence intensity
    % intScore = cellfun(@(x) mean(img(x)),CCEdges.PixelIdxList);
    
    % most probable paths have high response steerable filter
    % response values.
    responseMap = backboneInfo(iFrame).maxNMSLarge; % currently do NOT save the full
    % response in the backboneInfo need to check if this is more
    % helpful.
    
    resScore = cellfun(@(x) responseMap(x),CCEdges.PixelIdxList,'uniformoutput',0);
    resScore = cellfun(@(x) mean(x(x~=0)),resScore); % take out zer values from NMS.
    
    resScore = resScore./max(resScore);
    
    
    % Make Trouble Shooting Plots for the Steerable Filter Response
    % Score (Mean of Path)
    if p.plots == 1
        TBResPath =  [p.OutputDirectory filesep 'TroubleShootResScores'];
        if ~isdir(TBResPath)
            mkdir(TBResPath)
        end
        setFigure(nx,ny,'on');
        
        imagesc(responseMap)
        hold on
        text(20,20,'Trouble Shoot Response Score');
        colorbar
        % plot the outline of the paths considered
        roiYX = bwboundaries(dilBB);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
        roiYXPieces = bwboundaries(newBodyMask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYXPieces);
        % document response score for each path
        centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
        arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(resScore(i),3),'color','w'),1:length(centers));
        saveas(gcf,[TBResPath filesep num2str(iFrame,'%03d') '.tif']);
        saveas(gcf,[TBResPath filesep num2str(iFrame,'%03d') '.fig']);
        close gcf
    end % p.plots == 1
    
    
    
    
    
    % most probable paths have larger scale ridges (at least when the
    % scale estimate is working correctly - currenlty there seems to
    % think the max response is defaulting to max scale tested due
    % to either potential bug or something have to actually work
    % out.
    scaleMap = backboneInfo(iFrame).scaleMapLarge ;
    scaleScore =  cellfun(@(x) mean(scaleMap(x)),CCEdges.PixelIdxList);
    scaleScore = scaleScore./max(scaleScore); % make between 0 and 1
    
    
    
    
    % Make Trouble Shooting Plots for the Steerable Filter Response
    % Score (Mean of Path)
    if p.plots == 1
        TBScalePath =  [p.OutputDirectory filesep 'TroubleShootScaleScores'];
        if ~isdir(TBScalePath)
            mkdir(TBScalePath)
        end
        setFigure(nx,ny,'on');
        
        imagesc(scaleMap)
        hold on
        text(20,20,'Trouble Shoot Scale Score');
        colorbar
        % plot the outline of the paths considered
        roiYX = bwboundaries(dilBB);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYXPieces);
        % document response score for each path
        centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
        arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(scaleScore(i),3),'color','w'),1:length(centers));
        saveas(gcf,[TBScalePath filesep num2str(iFrame,'%03d') '.tif']);
        saveas(gcf,[TBScalePath filesep num2str(iFrame,'%03d') '.fig']);
        close gcf
    end % p.plots == 1
    
    %% add intensity information
    intScore = cellfun(@(x) img(x),CCEdges.PixelIdxList,'uniformoutput',0);
    intScore = cellfun(@(x) mean(x(x~=0)),intScore); % take out zer values from NMS.
    
    intScore = intScore./max(intScore);
    
    % Make Trouble Shooting Plots for the Steerable Filter Response
    % Score (Mean of Path)
    if p.plots == 1
        TBIntPath =  [p.OutputDirectory filesep 'TroubleShootIntensity'];
        if ~isdir(TBIntPath)
            mkdir(TBIntPath)
        end
        setFigure(nx,ny,'on');
        
        imshow(-img,[])
        hold on
        text(20,20,'Trouble Shoot Scale Score');
        colorbar
        % plot the outline of the paths considered
        roiYX = bwboundaries(dilBB);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','y'),roiYX);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','r'),roiYXPieces);
        % document response score for each path
        centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
        arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(intScore(i),3),'color','w'),1:length(centers));
        saveas(gcf,[TBIntPath filesep num2str(iFrame,'%03d') '.tif']);
        saveas(gcf,[TBIntPath filesep num2str(iFrame,'%03d') '.fig']);
        close gcf
    end % p.plots == 1
    
    
    
    
    
    finalScore = scaleScore + resScore +intScore; % for now just use a scale Score
    %%
    % TEST FOR PARALLEL EDGES
    edgesStr  = cellfun(@(x) num2str(x),edges,'uniformoutput',0); % put into string to use unique
    % find all parallel edges - stupid way for sure but lets do it
    % for now.
    [test,ic,iR]  = unique(edgesStr);
    if length(edgesStr) ~= length(test); % test for repeats.
        display('Parallel Connections To Same Node Found: Fixing');
        
        uiR = unique(iR); % get the indexes of the potential repeats
        idxDiscard = cell(length(uiR),1);
        for iPotRepeat = 1:length(uiR)
            
            repeatTest = sum(iR == uiR(iPotRepeat)); % problem if there are two repeats Fixed 20141129.
            if repeatTest > 1; % parallel edge
                % get the score for each edge and choose the max % or
                % could possibly just
                idxTest  = find(iR == uiR(iPotRepeat)); % indices of th
                % getScores or repeats
                scoresRepeats = finalScore(iR==uiR(iPotRepeat));
                maxScoreInGroup = max(scoresRepeats);
                idxDiscard{iPotRepeat} = idxTest(scoresRepeats~=maxScoreInGroup);
            end  % repeat test
        end %
        idxDiscard = vertcat(idxDiscard{:}); % changed from horzcat... 20141207
        % if check
        check = 1;
        if check == 1
            imshow(img,[]);
            hold on
            edgeMask = zeros(size(img));
            edgeMask(vertcat(CCEdges.PixelIdxList{:})) = 1;
            spy(edgeMask,'b');
            edgeMask(vertcat(CCEdges.PixelIdxList{idxDiscard}))= 0 ;
            spy(edgeMask,'r');
        end
        
        % discard that edge
        edges(idxDiscard) = [];
        dilBB(vertcat(CCEdges.PixelIdxList{idxDiscard}))=0;
        CCEdges.PixelIdxList(idxDiscard) = [];
        CCEdges.NumObjects = CCEdges.NumObjects - length(idxDiscard);
        finalScore(idxDiscard) = [] ;
        % CCEdges.NumObjects = CCEdges.
        %end % repeat test
        % end % iPotRepeat
        % for iRepeat = 1:length(uiR)
        % end
        %if length(iRepeat
        
        % seem to have a problem with putting parallel edges
        % set up weights
        
        
    end
    prefill = (newBodyMask | dilBB);
    % try again to fix the cycle
    
    
    % retest for cycles.
    
    % set up as a min span tree (need to make a sparse array)
    
    % check for parallel edges by finding the edges repeats in
    % the edge map
    % get all repeats. choose the max score  path here.
    %Try to test again
    
    fullMask = imfill(prefill,'holes');
    diffMask= fullMask-prefill;
    numBodyNodes = cellfun(@(x) length(x),edges); % as I don't treat
    % the surrounding frame pixels as a body "node" sometimes these
    % an path can connect to only 1 body piece this will cause the
    % following steps to crash - therefore select for only those
    % paths that connect two well-defined body pieces)
    edges = edges(numBodyNodes ==2) ;
    finalScore= finalScore(numBodyNodes==2);
    CCEdges.PixelIdxList = CCEdges.PixelIdxList(numBodyNodes==2);
    CCEdges.NumObjects = CCEdges.NumObjects - sum(numBodyNodes==2) ;
    
    if (~isequal(fullMask,prefill) && sum(diffMask(:))>2 && ~isempty(edges)); % make the size
        % slightly larger as cycle test currently based on simple
        % fill criterion therefore not that stable. see if can make
        % more stable.
        % if ~isempty(edges) % again simply check if it is a viable cycle-
        % again this is a weakness in the way I implemented
        % this..
        cycleFlag = 2; % non-parallel cycles
        % perform minspanning tree..
        display('cylces NOT due to parallel paths: performing minspantree');
        % put the information into a sparse mtrix.
        vect1 = cellfun(@(x) x(1),edges);
        vect2 = cellfun(@(x) x(2),edges);
        % make so lower score more favorable but do not have a zero weight
        % as a sparse array as will remove that edge completely
        finalScore = max(finalScore) - finalScore + 0.01;
        UG = sparse([vect1 vect2], [vect2 vect1], [finalScore finalScore] ); % need to make undirected.
        gFinal =  graphminspantree(UG);
        weightsFinal = full(gFinal(gFinal~=0));
        % delete edges given min span tree output.
        weightDelete = setdiff(finalScore,weightsFinal);
        idxDiscard   = arrayfun(@(i) find(finalScore == weightDelete(i)),1:length(weightDelete));
        if p.plots == 1
            % create the directory if doesn't exist
            outPathTBTree = [p.OutputDirectory filesep 'MinSpan Tree Weights'];
            if ~isdir(outPathTBTree)
                mkdir(outPathTBTree)
            end
            
            setFigure(nx,ny,'on')
            imshow(-img,[]);
            hold on
            
            spy(dilBB,'b');
            hold on
        end
        % delete the edge from the dilated backbone mask
        dilBB(vertcat(CCEdges.PixelIdxList{idxDiscard}))= 0 ; %
        if p.plots == 1
            
            spy(dilBB,'r');
            cellfun(@(x) plot(x(:,2),x(:,1),'color','y'),roiYXPieces);
            saveas(gcf, [outPathTBTree filesep num2str(iFrame,'%03d') '.tif']);
            saveas(gcf,[outPathTBTree filesep num2str(iFrame,'%03d') '.fig']);
            
        end
        % discard the appropriate edges from the list (ccs of the dilBB mask )
        CCEdges.PixelIdxList(idxDiscard) = [];
        CCEdges.NumObjects = CCEdges.NumObjects - length(idxDiscard);
        
        fullMaskPreFill = newBodyMask | dilBB;
        % test for cycles one last time. - if still have cycles can flag
        fullMask = imfill(fullMaskPreFill,'holes');
        if ~isequal(fullMaskPreFill,fullMask)
            display('You still have a hole that the algorithm cannot currently resolve: this frame will be flagged as a low confidence neurite body reconstruct');
            cycleFlag= 3;
        end
    end % second check for cycles (before minSpanTree)
end % all tests for cycles
%%
roiYXAll = bwboundaries(fullMask);
pixIndAll =  sub2ind(size(img),roiYXAll{1}(:,1),roiYXAll{1}(:,2));
pixIndThickBody = intersect(pixIndAll,pixIndThickBody);
pixIndThinBody = intersect(pixIndAll,pixIndThinBody);

thickBodyMask = zeros(size(img));
thickBodyMask(pixIndThickBody) = 1;
thickBodyMask = imfill(thickBodyMask,'holes');
analInfo.masks.thickBodyMask = thickBodyMask;
%%
if paramsIn.makeMovie ==1
    close gcf
    setFigure(nxLarge,nyLarge,'off');
    imshow(-imgLarge,[]);
    hold on
    roiYXNB = bwboundaries(fullMask);
    cellfun(@(x) plot(x(:,2),x(:,1),'color','y'),roiYXNB);
    pixels = round(10/pixSizeMic);
    plotScaleBar(pixels,pixels/20,'Color',[0 0 0]);
    saveas(gcf,[saveDirMov filesep '06.png'])
end
%%
analInfo.masks.neuriteEdge = fullMask;
% update so can restart if crash (as for now not saving as .tif (could
% change that- wanted to keep the filenames adaptable..)
analInfo.idxEnterNeurite = idxEnterNeurite; % save this information to test for later;

analInfo.bodyEst.pixIndThickBody = pixIndThickBody;
analInfo.bodyEst.pixIndThinBody = pixIndThinBody; % think about how you want to store this.
analInfo.cycleFlag = cycleFlag;

clear intAvg
if paramsIn.plots ==1
    
    hTroubleshoot =  GCAVisualsMakeTroubleshootVeilStemOverlay(img,analInfo,backboneInfo);
else
    hTroubleshoot = [];
end
end

