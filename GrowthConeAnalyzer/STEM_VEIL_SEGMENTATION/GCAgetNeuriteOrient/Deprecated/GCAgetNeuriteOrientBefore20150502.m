function [ backboneInfo] = GCAgetNeuriteOrient(img,p,iFrame)
% GCAgetNeuriteOrient: (STEP I of GCA package) 
% This function finds the large scale ridge that is most likely to
% correspond to the point of entry of the neurite in the image.
%
%% INPUT:
%   img: (REQUIRED) a RxC array of the image where N and M correspond to the image size
%  
% 






% EXTENDED SUMMMARY:
% Neurite backbone seed candidates are obtained via a 4th order steerable
% filtering process followed by non-maximal suppression.  
% As these signals can still be somewhat noisy, some ridge cleaning/linking steps are
% employed. These steps slightly less elegant. 
% 
%%% Steps in ridge Cleaning
% 1. Ridges that are associated with a very conservative background
% estimate based on intensity are removed. 
% 2. Thresholding of Very Low Response is performed 
% 3. Junction break
% 4. A small size filter (potentialy to be made optional or to include all) and
%    a linking step to link co-linear ridges within a certain distance.
% 5. Small ridge linking step.
% 6. Another small size filter for ridges post ridge linking 

% All ridge candidates with an endpoint within a user defined radius 
% (default 10) pixels of the image boundary are considered backbone seed 
% candidates, if no candidates are
% found this search radius is widened to another user selected value
% (default 20 pixels)
%
% A simple length and intensity linear cost for each candidate
% ridge is calculated to determine the probability of being a true backbone
% seed:
%
% Higher average intensities and longer ridges are considered more probable
% true backbone seeds.
%
% Currently the algorithm for veil/stem estimation is iterative and builds
% upon this initial backbone seed to link pieces of
% thicker more amorphous veil structures.  backboneInfo.after



%
% OUTPUT
% backboneInfo: A structure with perhaps a bit too much information
% includes the mask of the chosen seed candidate.
% as well as other potential candidates so that temproal information via a
% majority voting can be used to facilitate the candidate seed selection
% the next step.






%% CHECK INPUT
if nargin<2
    % Set Defaults.
    
    p.plots = 1;
    % Specific
    p.BBScale = 5:10; % PERSONAL NOTE: not sure if this is optimal
    p.FilterOrderBB = 4;
    p.OutputDirectory = pwd;
    p.MaxRadiusLargeScaleLink = 10; 
    p.ThreshNMSResponse = 25; % 
    p.MinCCRidgeBeforeConnect = 3 ; % in pixels 
    p.MinCCRidgeAfterConnect = 5; % in pixels (not sure why I need this step try to remove?) 
   
    
    p.MinSizeForEntranceRidge = 10 ; % in pixels - is adaptive however and will shrink if no candidates within that size range are found 
    p.maxDistanceFromEdge1 = 10;  % in pixels 
    p.maxDistanceFromEdge2  = 20; % in pixels 
end
saveDir = p.OutputDirectory;
[ny,nx] = size(img);
imSize = size(img);
%%
linked =0 ; % initiate small flag for plotting links in the colinear link step

% Calculate LargeScale Ridges
[maxResLarge, maxThLarge ,maxNMSLarge ,scaleMapLarge]= gcaMultiscaleSteerableDetector(img,p.FilterOrderBB,p.BBScale);

% Estimate the background of the image based on intensity
[maskBack,~,~] = gcaEstimateBackgroundArea(img);

% only use responses that are associated with the fluorescence background
maxNMSLarge = maxNMSLarge.*~maskBack ; %  only take responses with high fluorescence/image intensity

backboneInfo.scaleMapLarge = scaleMapLarge;
if p.plots == 1
    figDirScale = [saveDir filesep 'ScalesOfRidges'];
    if ~isdir(figDirScale)
        mkdir(figDirScale)
    end
    figScale = setFigure(nx,ny,'off');
    
    imagesc(scaleMapLarge.*(maxNMSLarge>0));
    colorbar
    saveas(figScale,[figDirScale filesep num2str(iFrame,'%03d') '.tif']);
    close gcf
end
% Threshold the NonMaximalSuppression Ridge Response Values
values = maxNMSLarge(maxNMSLarge~=0);
backboneInfo.maxNMSLarge = maxNMSLarge;

cutoff = prctile(values,p.ThreshNMSResponse); % note arbitrarily set this cut-off to the 25th percentile- tried to fit the first mode of
backboneInfo.bodyEst.cutoff = cutoff;
% this distribution to a Gaussian and take X STD above the mean,
% however this was typically way to conservative- not much signal
% beyond that of the gaussian fit - therefore just use this for
% now to clean out small amounts of noise.
%%% NOTE TO SELF : % REVISIT ABOVE


% get ridge candidates
ridgeCand = maxNMSLarge>cutoff;
ridgeCand(ridgeCand>0)=1;
% break all junctions
ridgeCand = bwmorph(ridgeCand,'thin','inf');


%% Start Optional TroubleShoot Plot I: RidgeCand Before After Clean
if p.plots == 1
    figDir1 = [saveDir filesep 'RidgeCandBeforeAfterClean'];
    if ~isdir(figDir1);
        mkdir(figDir1);
    end
    
    [ny,nx]= size(img);
    figRidgeCand = setFigure(nx,ny,'off');
    imshow(img,[]);
    hold on
    spy(ridgeCand,'b');
    hold on
end
%% Make Boundary Mask
% there are typically edge effects around the border don't consider
% response at the border

boundaryMask = zeros(imSize);
boundaryMask(1:imSize(1),1) =1;
boundaryMask(1:imSize(1),imSize(2))=1;
boundaryMask(1,1:imSize(2))= 1;
boundaryMask(imSize(1),1:imSize(2)) =1;

%% CLEAN RIDGES (the maxNMS)
% STEP I: Cut Junctions in the maxNMS
nn = padarrayXT(double(ridgeCand~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn1 = (nn-1) .* (ridgeCand~=0);
junctionMask = nn1>2;
ridgeCand(junctionMask) =0;



%Consider and ridges along the border likely edge effects - Make
%sure to take care of before
ridgeCand(:,1:3) =  0;
ridgeCand(:,end-2:end) = 0;
ridgeCand(1:3,:) = 0;
ridgeCand(end-2:end,:) = 0;



% STEP II: Clean Pieces Less than X pixels
% Clean small pieces less than X pixels: again this might lead you to lose info- might want to reduce to singletons
 % with singletons it is hard to tell a direction for linking could include but currently better
% I think to exclude.
CCRidgeBone = bwconncomp(ridgeCand);
% filter out noise (defined here simply as CCs less than 5 pixels
% (arbitrary assignment)
csizeRidgeFirst = cellfun(@(x) length(x),CCRidgeBone.PixelIdxList);
CCRidgeBone.PixelIdxList(csizeRidgeFirst < p.MinCCRidgeBeforeConnect) = [];
CCRidgeBone.NumObjects = CCRidgeBone.NumObjects - sum(csizeRidgeFirst<X);
cleanedRidgeLabels = labelmatrix(CCRidgeBone);
cleanedRidge = cleanedRidgeLabels>0;




%%  Finish Optional TroubleShoot Plot I: RidgeCand Before After Clean
if p.plots == 1
    spy(cleanedRidge,'g');
    saveas(figRidgeCand, [figDir1 filesep num2str(iFrame,'%03d') '.tif']);
    close gcf
end


%% Set-up to Perform "Linear" Links-  MaxWeighted GraphMatch

CCRidgeBonePreLink = bwconncomp(cleanedRidge);
cleanedRidgeLabelsPreLink = labelmatrix(CCRidgeBonePreLink);

%%% Get endpoints  %%%- note singletons filtered here above so no
%%% need to check for them.  20140922 - fix this for final
EPCandidateSort = cellfun(@(x) getEndpoints(x,size(img),0,1),CCRidgeBonePreLink.PixelIdxList,'uniformoutput',0);
backboneInfo.beforeConnect= cleanedRidge;


%% Start Plot 2: Before and after connect
if p.plots == 1
    figDir2 = [saveDir filesep 'BeforeAndAfterConnect'];
    if ~isdir(figDir2)
        mkdir(figDir2)
    end
    figConnect = setFigure(nx,ny,'off');
    hold on
    imshow(img,[]) ;
    hold on
    % plot ridge before connect
    spy(cleanedRidge)
    hold on
    
    test = vertcat(EPCandidateSort{:});
    
end
%% Connect Linear Structures % NOTE to SELF as of 2014-01-30 : this part needs to be optimized
% main problems fixed 20141021

[cleanedRidge,linkMask,~,~,madeLinks] = gcaConnectLinearRidges(EPCandidateSort,maxThLarge,cleanedRidge,cleanedRidgeLabelsPreLink, p.MaxRadiusLargeScaleLink);%NEED to make a variable!
cleanedRidge = bwmorph(cleanedRidge,'thin');  % after do this type of connect always need to thin!
%% STEP TO POTENTIALLY REMOVE: 
% I believe I had this step here originally as it was after a junction
% break- all data run before 04-27-2015 will have included this step 
% STEP II: Clean Pieces Less than 5 pixels
% Clean small pieces less than 5 pixels: again this might lead you to lose info- might want to reduce to singletons
CCRidgeBone = bwconncomp(cleanedRidge);
% filter out noise (defined here simply as CCs less than 5 pixels
% (arbitrary assignment)
csizeRidgeFirst = cellfun(@(x) length(x),CCRidgeBone.PixelIdxList);
CCRidgeBone.PixelIdxList(csizeRidgeFirst < ccMinRemoveAfterConnect) = [];
CCRidgeBone.NumObjects = CCRidgeBone.NumObjects - sum(csizeRidgeFirst< p.MinCCRidgeAfterConnect);
cleanedRidgeLabels = labelmatrix(CCRidgeBone);
cleanedRidge = cleanedRidgeLabels>0;

%% Finish Plot 2: Before and after connect
if madeLinks == 1
    spy(linkMask,'r');
end
scatter(test(:,1),test(:,2),5,'y','filled'); % scatter end points
saveas(figConnect,[figDir2 filesep num2str(iFrame,'%03d') '.tif']);

%% Find the Ridge components that are near the boundary of the cell: choose the most likely candidate
% STEP I. Get connected components of all ridges and make label matrix
CCRidgeBoneLinked = bwconncomp(cleanedRidge);
cleanedRidgeLinkedLabelMat = labelmatrix(CCRidgeBoneLinked);
cleanedRidgeLinked = cleanedRidgeLinkedLabelMat>0;



backboneInfo.bodyReconstruct.AfterConnect = cleanedRidgeLinked;
%STEP III: Get all maxNMSRidge pixels that are along or very near border

% perform distance transfer from boundary
distTransBound =  bwdist(boundaryMask);
add = 0 ;
CCBorderRidgeCandPixels = [];

while isempty(CCBorderRidgeCandPixels)
    
    
    % get the indices of ridge pixels within 4 pixels of boundary
    idxRidgeBorder_Logic = (distTransBound<(10+add) & cleanedRidgeLinked ==1); % remember filtered 3 pixels of the maxNMS around border
    
    % find the connected component ridge labels associated with these pixels
    ridgeLabelsBorder =  unique(cleanedRidgeLinkedLabelMat(idxRidgeBorder_Logic));
    
    % test these tiny candidates
    
    
    
    % test these border ridges for liklihood of being the primary enterning
    % neurite.
    
    CCBorderRidgeCandPixels = CCRidgeBoneLinked.PixelIdxList(ridgeLabelsBorder);
    
    add = add + 1;
end

%
% get rid of tiny candidates

x =0;
while x == 0
    
    csizeCand = cellfun(@(x) length(x), CCBorderRidgeCandPixels);
    x = sum(csizeCand>minSizeForEntranceRidge);
    cutoff = cutoff-1;
end

CCBorderRidgeCandPixels = CCBorderRidgeCandPixels(csizeCand>cutoff);
% save candidateSeeds
candidateSeeds = zeros(size(img));
candidateSeeds(vertcat(CCBorderRidgeCandPixels{:})) = 1;
backboneInfo.candSeeds = candidateSeeds; % need to save these as they are possible new seeds
% quickFix
%         if isempty(CCBorderRidgeCandPixels)
%             CCBorderRidgeCandPixels = CCRidgeBoneLinked.PixelIdxList(ridgeLabelsBorder);
%            csizeCand = cellfun(@(x) length(x),CCBorderRidgeCandPixels);
%            CCBorderRidgeCandPixels = CCBorderRidgeCandPixels(csizeCand>15);
%         end
%% Start Plot III:
% plot candidates
if p.plots == 1
    figSeed = setFigure(nx,ny,'off');
    figDir3 = [saveDir filesep 'CandSeeds'];
    if ~isdir(figDir3)
        mkdir(figDir3)
    end
    hold on
    imshow(img,[]);
    hold on
    candidateSeeds = zeros(size(img));
    candidateSeeds(vertcat(CCBorderRidgeCandPixels{:})) = 1;
    spy(candidateSeeds,'g');
    
    
end

%% Find the Best Entering Ridge Candidate
% filter out those pieces less than 10 pixels

%
% STEP II: Check Average Intensity: Largest intensity likely candidate
avgIntCandBorderRidge = cellfun(@(x) mean(img(x)), CCBorderRidgeCandPixels);
csize = cellfun(@(x) length(x),CCBorderRidgeCandPixels);
score = (avgIntCandBorderRidge/max(avgIntCandBorderRidge)+csize/max(csize)); % make a linear cost

% keep the one with the largest average intensity
backboneSeed = zeros(size(img));

pixBackbone = CCBorderRidgeCandPixels{score==max(score)} ; % changed 12/11/2013 from original of intensity based only
%pixBackbone= CCBorderRidgeCandPixels{csize == max(csize)};
backboneSeed(pixBackbone) =1 ;
%% Add to Plot III: Border Ridge Candidates and Final Backbone

spy(backboneSeed,'r');
%%
% re-pad the image
%     backboneSeed =  padarray(backboneSeed,[3,3]) ;
%     boundaryMask = padarray(boundaryMask,[3,3]);
%     cleanedRidge = padarray(cleanedRidge,[3,3]);
%

% save backboneSeed


%% Interpolate to the edge if necessary
% test if backbone seed gives you a closed contour.
idxEnterNeurite = find(backboneSeed ==1 & boundaryMask ==1);
[yEnter,xEnter] = ind2sub(size(img),idxEnterNeurite);
if isempty(idxEnterNeurite) % need to interpolate to make sure closed contour
    % interpolate between to nearest point on boundary
    % find endpoint of backboneSeed
    
    EPsBackbone = getEndpoints(pixBackbone,size(img));
    [yBoundary,xBoundary] = find(boundaryMask==1);
    [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary] ,EPsBackbone,p.maxDistanceFromEdge1);
    % quick fix
    %        if isempty(idxBoundaryClose)
    %            [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],EPsBackbone,20);
    %        end
    distAll =  vertcat(dist{:});
    useEPsFlag = 1; % use endpoint flag
    numPixBBTested = length(EPsBackbone(:,1));
    % test if no boundary pixels are within 10 pixels of the end point
    if isempty(distAll)
        % use all pixels along edge
        [yBBAll, xBBAll] = ind2sub(size(img),pixBackbone);
        useEPsFlag = 0 ; % use all the coordinates of the BB and try again
        [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],[xBBAll,yBBAll],p.maxDistanceFromEdge1);
        numPixBBTested = length(xBBAll);
        distAll =vertcat(dist{:});
        if isempty(distAll)
            [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],EPsBackbone,p.maxDistanceFromEdge2); % try a larger search radius
            distAll = vertcat(dist{:});
            numPixBBTested = length(EPsBackbone(:,1));
        end
    end
    if ~isempty(distAll);
        
        
        
        
        
        
        % find the minimum distance
        toSave = find(distAll ==min(distAll));
        if length(toSave) > 1
            toSave = toSave(1);
        end
        
        
        
        E = arrayfun(@(i) [repmat(i, [numel(idxBoundaryClose{i}) 1]) idxBoundaryClose{i}], 1:numPixBBTested, 'UniformOutput', false);
        E = vertcat(E{:});
        
        if useEPsFlag == 1 ;
            
            BBLinkers =  [EPsBackbone(E(toSave,1),1),EPsBackbone(E(toSave,1),2)];
        else
            BBLinkers = [xBBAll(E(toSave)),yBBAll(E(toSave))];
        end
        
        linkCoords = gcaBresenham([BBLinkers(:,1),BBLinkers(:,2)],[xBoundary(E(toSave,2)),yBoundary(E(toSave,2))]);
        backboneSeed(linkCoords(:,2), linkCoords(:,1) ) = 1;
        linked = 1;
        xEnter = xBoundary(E(toSave,2));
        yEnter = yBoundary(E(toSave,2));
        %         % if quick fix
        %         useEPs == 0 ;
        %         % break junctions.
        
        
    else % if it is still empty you have a problem skip the next steps
        warning('GCA:getNeuriteOrient:FailureToCloseContour',['Contour did NOT close for Frame ' ]);
    end
end % if isempty idxNeurite
%% Finish Plot III: Border Ridge Candidates and Final Backbone
if p.plots == 1
    spy(backboneSeed,'r');
    if linked == 1 ;
        scatter(linkCoords(:,1),linkCoords(:,2),'y','filled');
    end
    scatter(xEnter,yEnter,'m','filled');
    saveas(gcf,[figDir3 filesep num2str(iFrame,'%03d') '.tif']);
    gcf
    linked = 0; % reset
end
%% Save the information
backboneInfo.backboneSeedMask = backboneSeed;
backboneInfo.coordsEnterNeurite = [xEnter,yEnter]; % save this to test for consistency
backboneInfo.timeStamp = clock;


end

