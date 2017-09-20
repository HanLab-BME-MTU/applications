function [ backboneInfo,TSFigs] = GCAgetNeuriteOrient(img,varargin)
% GCAgetNeuriteOrient: (STEP I of GCA PACKAGE)
% This function finds the large scale ridge that is most likely to
% correspond to the point of entry of the neurite in the image, hence
% detecting the neurite's orientation in the image. It likewise
% generates and cleans the other large scale ridge detections that will be
% used in GCAVeilStemReconstructMovie.m
%
%% INPUT:
%
%   img: (REQUIRED) : RxC double array
%        of image to analyze where R is the height (ny) and C is the width
%       (nx) of the input image
%
%% PARAMS
%
% %% PARAMS: STEERABLE FILTER: RIDGE FINDING %%
%
%    'BBScale' (PARAM) : Positive scalar or vector
%        Sigmas (standard deviation of the Gaussian kernel) to use for
%        the steerable filter estimation of the backbone (BB). Note if a vector
%        is specified, responses for all scales are calculated and the
%        scale with the largest steerable filter
%        response at each point is chosen for the final backbone estimation.
%        Default: [5,6,7,8,9,10] (in pixels)
%        See gcaMultiscaleSteerableDetector.m
%
%    'FilterOrderBB' (PARAM) : Scalar 2 or 4
%        4 provides better orientation selectivity than 2
%        and is less sensitive to noise, at a small trade-off in
%        computational cost. 
%        Default: 4
%        See gcaMultiscaleSteerableDetector.m
%
% %% PARAMS: RIDGE LINKING %%
%
%    'MaxRadiusLargeScaleLink' (PARAM) :  Positive Scalar
%        Radius for linking neighboring large-scale ridge candidates
%        Default: 10 (in pixels)
%        See gcaConnectLinearRidges.m
%
% %% PARAMS: RIDGE CLEANING %%
%
%     'ThreshNMSResponse'  (PARAM) :  Postive Scalar
%        Percentile for the NMS response threshold cut-off
%        Default: 25
%
%     'MinCCRidgeBeforeConnect' (PARAM) : Positive Scalar
%         Size of the connected components of the ridge NMS to filter
%         BEFORE linking, note for small candidate ridges it is difficult
%         to maintain a good orientation definition for linking.
%         Default: 3 (in pixels)
%
%     'MinCCRidgeAfterConnect' (PARAM) : Positive Scalar
%         Size of the connected components of the ridge NMS to filter
%         AFTER ridge linking.
%         Default: 5 (in pixels)
%
% %% PARAMS: DEFINE NEURITE ENTRANCE CANDIDATES %%
%
%     'MinSizeForEntranceRidge' (PARAM) : Positive Scalar
%         Smallest Size for a ridge connected component to be considered a
%         neurite entrance candidate ridge. Adaptive: if none found will
%         search iteratively for candidates of smaller lengths until it
%         finds a candidate (ie assumes a correct candidate exists)
%         Default: 10 (in pixels)
%
%     'MaxDistBorderFistTry' (PARAM) : Positive Scalar
%         Any ridges with endpoints within this value from the border will
%         be considered ridge candidates. If no candidates are found this
%         cut-off will be iteratively opened until a candidate entrance ridge
%         is found. Therefore, it always assumes a neurite entrance ridge should
%         exist.
%         Default: 10 (in pixels)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%     'backboneInfo': A structure with fields :
%        .backboneSeedMask : an rxc logical mask where r is the height (ny)
%               of the input img and c is the width (nx)
%               marking the ridge that is most likely the neurite entrance.
%               These candidates will be tested for temporal outliers in
%               GCAneuriteOrientConsistencyCheck.m, and these outliers will be
%               corrected.
%               These masks serve as the 'seed' in the first step of
%               GCAveilStemReconstruct.m
%
%       .coordsEnterNeurite : an 1x2 double
%           giving the xy coordinates corresponding to where the neurite
%           entrance point was detected. If there are temporal outliers
%           these values may be refined in
%           GCAneuriteOrientConsistencyCheck.m
%
%       .candSeeds : an rxc logical mask where r is the height (ny) of the
%           input img and c is the width (nx)
%           marking any alterative seed ridges that may be considered as
%           entrance cadidates.
%           These other ridges will be re-considered if the selected ridge
%           was found to be a temporal outlier in the
%           GCAneuriteOrientConsistencyCheck.m
%
%       .linkedRidgesFinal: an rxc logical mask where r is the height (ny)
%           of the input img and c is the width (nx) marking cleaned large
%           scale ridge candidate paths that will be used for the
%           final reconstruction in GCAveilStemReconstruct.m
%
%       .scaleMapLarge
%       .maxNMSLarge
%
%
%   TSFigs (OPTIONAL) : rx1 structure with fields where r is the number of troubleshoot
%                     figures
%                   .h : figure handle
%                   .name : character array containing the name of the figure
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTENDED SUMMMARY:
%
% I.Neurite backbone seed candidates are obtained via a 4th order steerable
%   filtering process followed by non-maximal suppression.
%
% II.As these signals can still be somewhat noisy, some non-elegant
%    ridge cleaning/linking steps are employed:
%
% Steps in Ridge Cleaning:
%    1. Ridges that are associated with a very conservative background
%       estimate based on fluorescence intensity are removed.
%    2. Thresholding of Bottom Percentile of the NMS Response (typically <25 th) is performed
%    3. Any junctions are broken.
%    4. Edge effects along the border are removed
%    5. Small connected components < a user defined pixel value are removed
%       as small candidates tend to not have high fidelity orientation information
%       for linking.
%    6. Ridge linking step.
%    7. Another small filter for connected components post linking. (likely
%    can remove...MB Test)
%
% III. Find Candidate Neurite Entrance Ridges
%    All ridge candidates with at least one point within a user defined radius
%   (default 10) pixels of the image boundary are considered backbone seed
%   candidates, if no candidates are found this search radius is widened
%   until a candidate is found. (MB 20150428 - likely want to switch this
%   so it is looking for the endpoints - just don't want to change it now
%   as did a lot of testing with this framework
%
% A simple length and intensity linear cost for each candidate
% ridge is calculated to determine the probability of being a true backbone
% seed.
%
% Higher average intensities and longer ridges are considered more probable
% true backbone seeds and thus maintained: currently the package
% only considers temporal constraints in the second step
% GCAneuriteOrientConsistencyCheck.m
% 
% The shortest path between the endpoint of the candidate backbone and the
% boundary is interpolated : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 
%REQUIRED
ip.addRequired('img');

% PARAMETERS
ip.addParameter('TSOverlays',true,@(x) islogical(x));

ip.addParameter('BBScale',[5 6 7 8 9 10]);
ip.addParameter('FilterOrderBB',4,@(x) ismember(x,[2 4]));

ip.addParameter('MaxRadiusLargeScaleLink',10) ;
ip.addParameter('NoLinkDistanceFromBorder',10); 
ip.addParameter('MaxDistNoGeoTerm',3); 

ip.addParameter('ThreshNMSResponse',25);
ip.addParameter('MinCCRidgeBeforeConnect',3);
ip.addParameter('MinCCRidgeAfterConnect',5);

ip.addParameter('MinCCEntranceRidgeFirstTry',10);
ip.addParameter('MaxDistBorderFirstTry',10);
ip.parse(img,varargin{:});

%% Initiations

[ny,nx] = size(img);
imSize = size(img);
linked =0 ; % initiate small flag for plotting links in the colinear link step
if ip.Results.TSOverlays == true
    figCount = 1; % initiate counting of figure - makes easier if user wants to add
    % more figures later for troubleshooting if make generic.
    % set up colors
    cMap = brewermap(2,'paired');
    blue = cMap(2,:);
    
    % plot the alignment mask
    cDark = brewermap(2,'Dark2');
    orange = cDark(2,:);
    green = cDark(1,:);
    
else 
    TSFigs = []; 
end

%% START PART I. Detect large scale ridges and perform NMS
[~, ~ ,maxNMSLarge ,scaleMapLarge]= gcaMultiscaleSteerableDetector(img,ip.Results.FilterOrderBB,ip.Results.BBScale);

backboneInfo.scaleMapLarge = scaleMapLarge;

% END PART I

%% START PART II. Ridge Cleaning 

%% 1 Ridge Cleaning: Estimate Background Based On Intensity Values and Remove Ridges
% Associated With This Background
[maskBack,~,~] = gcaEstimateBackgroundArea(img);

maxNMSLarge = maxNMSLarge.*~maskBack ; %
backboneInfo.maxNMSLarge = maxNMSLarge; % save the maxNMS after removing background

%% Optional Troubleshoot (TS) Overlay 1:  Ridge Scale Overlays 
%  (using the current pre-cleaned NMS mask)

if ip.Results.TSOverlays == true
    TSFig1H = setFigure(nx,ny,'off');
    TSFigs(figCount).h  = TSFig1H;
    TSFigs(figCount).name = 'ScalesOfRidges';
    imagesc(scaleMapLarge.*(maxNMSLarge>0));
    colorbar
    figCount = figCount+1; % figure completed
end
%% 2 Ridge Cleaning: Threshold the NonMaximalSuppression Ridge Response Values:
%
% Note this thresholding is a bit arbitrary and one of the weaker steps in
% the algorithm. Response at this stage is typically almost Guassian
% Distributed therefore you might be able to argue that one may keep all of it
% MB: It might make sense to try to skip this step, though all successful
% segmentations as of 20150427 have been performed with this step included.
%
values = maxNMSLarge(maxNMSLarge~=0);
cutoff = prctile(values,ip.Results.ThreshNMSResponse); % note arbitrarily set this cut-off to the 25th percentile- tried to fit the first mode of
%backboneInfo.bodyEst.cutoff = cutoff; % save the information in the output

% get ridge candidates
ridgeCand = maxNMSLarge>cutoff;
ridgeCand(ridgeCand>0)=1;
% Perform a thinning operation to keep tidy
ridgeCand = bwmorph(ridgeCand,'thin','inf');
ridgeCandBeforeClean = ridgeCand; 
%% OPTIONAL Troubleshoot Overlay
if ip.Results.TSOverlays == true
  
    TSFigs(figCount).h  = setAxis('off');
     
    TSFigs(figCount).name = 'HistogramNMSResponse';
    
    [n,b] =  hist(values,100);
    hist(values,100); 
    hold on 
    line([cutoff cutoff],[0,max(n)],'color','r');  
    ylabel('Count'); 
    xlabel('Large Scale NMS Response Values'); 
    
    figCount = figCount + 1; % figure closed
end

%% 2 Ridge Cleaning: Threshold the NonMaximalSuppression Ridge Response Values:
%
% Note this thresholding is a bit arbitrary and one of the weaker steps in
% the algorithm. Response at this stage is typically almost Gaussian
% Distributed. MB: It might make sense to try to skip this step, though all successful
% segmentations as of 20150427 have been performed with this step included.
% 1
values = maxNMSLarge(maxNMSLarge~=0);
cutoff = prctile(values,ip.Results.ThreshNMSResponse); % note arbitrarily set this cut-off to the 25th percentile- tried to fit the first mode of
%backboneInfo.bodyEst.cutoff = cutoff; % save the information in the output
if cutoff <0 
    cutoff = 0 ; 
end 
% get ridge candidates
ridgeCand = maxNMSLarge>cutoff;
ridgeCand(ridgeCand>0)=1;
% Perform a thinning operation to keep tidy
ridgeCand = bwmorph(ridgeCand,'thin','inf');

%% 3 Ridge Cleaning: Break Junctions
% Cut Junctions in the maxNMS
nn = padarrayXT(double(ridgeCand~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn1 = (nn-1) .* (ridgeCand~=0);
junctionMask = nn1>2;
ridgeCand(junctionMask) =0;

%% 4 Ridge Cleaning : Filter Borders

%The Ridges Along the Borders are Tyically Edge Effects So Remove This
%Signal
ridgeCand(:,1:3) =  0;
ridgeCand(:,end-2:end) = 0;
ridgeCand(1:3,:) = 0;
ridgeCand(end-2:end,:) = 0;

%% 5 Ridge Cleaning : Remove small connected components < ip.Results.MinCCRidgeBeforeConnect
%
% Note: might want to reduce to singletons
% NOTE MB: you need to test this setting in the new run reduced to
% singletons. 20150427
CCRidgeBone = bwconncomp(ridgeCand);
% filter out small connected components
csizeRidgeFirst = cellfun(@(x) length(x),CCRidgeBone.PixelIdxList);
CCRidgeBone.PixelIdxList(csizeRidgeFirst < ip.Results.MinCCRidgeBeforeConnect) = [];
CCRidgeBone.NumObjects = CCRidgeBone.NumObjects - sum(csizeRidgeFirst<ip.Results.MinCCRidgeBeforeConnect);
cleanedRidgeLabels = labelmatrix(CCRidgeBone);
cleanedRidge = cleanedRidgeLabels>0;

%% Optional Troubleshoot (TS) Overlay II:: RidgeCand Before After Clean :
if ip.Results.TSOverlays == true
  
    TSFigs(figCount).h  = setFigure(nx,ny,'off');
     
    TSFigs(figCount).name = 'RidgeCandBeforeAfterClean';
    
    imshow(-img,[]);
    hold on
    spy(ridgeCandBeforeClean,'b'); % ridge candidates before connections
    hold on %   
    % finish overlay
    spy(cleanedRidge,'g'); % After Connection
    figCount = figCount + 1; % figure closed
end
%% 6. Ridge Cleaning : Ridge Linking Step:  PREPARE FOR LINKING

CCRidgeBonePreLink = bwconncomp(cleanedRidge);
cleanedRidgeLabelsPreLink = labelmatrix(CCRidgeBonePreLink);

%%% Get endpoints  %%%- note singletons filtered here above so no
%%% need to check for them.  20140922 - fix this for final
EPCandidateSort = cellfun(@(x) getEndpoints(x,size(img),0,1),CCRidgeBonePreLink.PixelIdxList,'uniformoutput',0);
backboneInfo.beforeConnect= cleanedRidge;

%% 6. Ridge Cleaning : Ridge Linking Step:  PERFORM LINKING

[cleanedRidge,linkMask,~,~,madeLinks] = gcaConnectLinearRidges(EPCandidateSort,cleanedRidgeLabelsPreLink,ip.Results);%NEED to make a variable!
% MB CHECK BEFORE RELEASE : part of the problem may be that you are
% introducing junctions at this stage... see how deal with in the next
% step. The two components you ~ out are the EPs of hte new CCs and the
% labelMatof the new CCs- you should really potentially put these new CCs
% into a cellarray so you can pick apart junctions.  
cleanedRidge = bwmorph(cleanedRidge,'thin');  % thinning is very important for subsequent steps

%% 7. Ridge Cleaning :  Remove small connected components < ip.Results.MinCCRidgeAfterConnect
% I believe I had this step here originally as it was after a junction
% break- all data run before 04-27-2015 will have included this step
% STEP II: Clean Pieces Less than 5 pixels
% again this might lead you to lose info- might want to exclude

CCRidgeBone = bwconncomp(cleanedRidge);
csizeRidgeFirst = cellfun(@(x) length(x),CCRidgeBone.PixelIdxList);
CCRidgeBone.PixelIdxList(csizeRidgeFirst < ip.Results.MinCCRidgeAfterConnect) = [];
CCRidgeBone.NumObjects = CCRidgeBone.NumObjects - sum(csizeRidgeFirst< ip.Results.MinCCRidgeAfterConnect);
cleanedRidgeLabels = labelmatrix(CCRidgeBone);
cleanedRidge = cleanedRidgeLabels>0;

%% Optional TroubleShoot Plot III: Before and after connect : AFTER
if ip.Results.TSOverlays == true
     TSFig3H = setFigure(nx,ny,'off');
    TSFigs(figCount).h  = TSFig3H;
    TSFigs(figCount).name = 'BeforeAndAfterConnect';
    
    imshow(-img,[]) ;
    hold on
    % plot ridge before connect
    spy(cleanedRidge)
    hold on
    test = vertcat(EPCandidateSort{:}); % figure open 
    
    
    if madeLinks == 1
        spy(linkMask,'r'); % plot any links made
    end
    scatter(test(:,1),test(:,2),5,'c','filled'); % scatter end points of ridges
    figCount = figCount +1; % figure closed
end % ip.Results.TSOverlays
%  END PART I
%% START PART III. SELECT NEURITE ENTRANCE RIDGE CANDIDATES

% Get connected components of all ridges and make label matrix
CCRidgeBoneLinked = bwconncomp(cleanedRidge);
cleanedRidgeLinkedLabelMat = labelmatrix(CCRidgeBoneLinked);
cleanedRidgeLinked = cleanedRidgeLinkedLabelMat>0; % this  may not be necessary check and remove

% Save for input to GCAveilStemReconstruct.m
backboneInfo.linkedRidgesFinal = cleanedRidgeLinked;

% Define the border of the image
boundaryMask = zeros(imSize);
boundaryMask(1:imSize(1),1) =1;
boundaryMask(1:imSize(1),imSize(2))=1;
boundaryMask(1,1:imSize(2))= 1;
boundaryMask(imSize(1),1:imSize(2)) =1;

% perform distance transfer from boundary
distTransBound =  bwdist(boundaryMask);

% First try to find candidates by distance from border : make the criteria
% less leniant if you do not find any candidates
CCBorderRidgeCandPixels = [];
distFromBorder = ip.Results.MaxDistBorderFirstTry;
while isempty(CCBorderRidgeCandPixels)
    
    % get the indices of ridge pixels ip.Results.MaxDistBorder1stTry
    idxRidgeBorder_Logic = (distTransBound<(distFromBorder) & cleanedRidgeLinked ==1); % remember filtered 3 pixels of the maxNMS around border
    
    % find the connected component ridge labels associated with these pixels
    ridgeLabelsBorder =  unique(cleanedRidgeLinkedLabelMat(idxRidgeBorder_Logic));
    
    CCBorderRidgeCandPixels = CCRidgeBoneLinked.PixelIdxList(ridgeLabelsBorder);
    
    distFromBorder = distFromBorder +1; 
end % while isempty



% filter out very small candidates unless it is the only one found within
% the distance, then keep it. 
x =0;
minEntranceRidgeSize  = ip.Results.MinCCEntranceRidgeFirstTry;
while x == 0
    % test all the candidates sizes
    csizeCand = cellfun(@(x) length(x), CCBorderRidgeCandPixels);
    % test if any of the candidates have a size greater than the cut-off
    % parameter
    x = sum(csizeCand>minEntranceRidgeSize);
    % if not adapt the cut-off to make sure at least one candidate is in
    % region. 
    minEntranceRidgeSize = minEntranceRidgeSize-1;
end % while x == 0

CCBorderRidgeCandPixels = CCBorderRidgeCandPixels(csizeCand>minEntranceRidgeSize);
% save candidateSeeds
candidateSeeds = zeros(size(img));
candidateSeeds(vertcat(CCBorderRidgeCandPixels{:})) = 1;
% save as input to GCAneuriteOrientConsistencyCheckMovie.m
backboneInfo.candSeeds = candidateSeeds; % need to save these as they are possible new seeds

%% Optional TroubleShoot Plot IV: Candidate Ridge Seeds : ALL POTENTIAL CANDS
if ip.Results.TSOverlays == true
    TSFig4H = setFigure(nx,ny,'off');
    TSFigs(figCount).h  = TSFig4H;
    TSFigs(figCount).name = 'CandSeeds';
    imshow(-img,[]);
    hold on
    candidateSeeds = zeros(size(img));
    candidateSeeds(vertcat(CCBorderRidgeCandPixels{:})) = 1;
    [nyCand,nxCand] = ind2sub(size(img),vertcat(CCBorderRidgeCandPixels{:})); 
%     spy(candidateSeeds,'g'); % figure open 
    scatter(nxCand,nyCand,10,blue,'filled'); 
end % ip.Results.TSOverlays
%% Find the Best Entering Ridge Candidate

% Check Average Intensity: Largest intensity likely candidate
avgIntCandBorderRidge = cellfun(@(x) mean(img(x)), CCBorderRidgeCandPixels);
csize = cellfun(@(x) length(x),CCBorderRidgeCandPixels);
score = (avgIntCandBorderRidge/max(avgIntCandBorderRidge)+csize/max(csize)); % make a linear cost

% keep the one seed candidate with the largest score 
backboneSeed = zeros(size(img));

pixBackbone = CCBorderRidgeCandPixels{score==max(score)} ; % changed 12/11/2013 from original of intensity based only

backboneSeed(pixBackbone) =1 ;

%% Interpolate to from the edge of the image to the candidate seed if necessary
% test if backbone seed gives you a closed contour- 
idxEnterNeurite = find(backboneSeed ==1 & boundaryMask ==1);
[yEnter,xEnter] = ind2sub(size(img),idxEnterNeurite);% this should actually
% never be the case anymore as we should have filtered out this information

if isempty(idxEnterNeurite) % need to interpolate to make sure closed contour
    % interpolate between to nearest point on boundary
    % find endpoint of backboneSeed
    
    EPsBackbone = getEndpoints(pixBackbone,size(img));
    [yBoundary,xBoundary] = find(boundaryMask==1);
    
    % find the closest distance between EP and boundary: 
    % MB note 20150428 if you set this up smarter such that you took only
    % the ridges with endpoints in a given region it might have been a bit
    % better as you would know that at least one endpoint is within the
    % distance you found above... anyway lets not shift things up here too
    % much as I tested it mainly with this framework. 
    
    % need to check which is faster KDTree or just a simple distance
    % transform... for now keep as is. MB not 20150428
    
    % First just test for the shortest distance from the backbone endpoints
    % to the boundary within the final open distance. 
    [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary] ,EPsBackbone,distFromBorder-1);
   
    distAll =  vertcat(dist{:});
    useEPsFlag = 1; % use endpoint flag
    numPixBBTested = length(EPsBackbone(:,1));
    % test if no EP boundary pixel is within this range just use the whole
    % thing. 
    if isempty(distAll)
        % use all pixels along edge
        [yBBAll, xBBAll] = ind2sub(size(img),pixBackbone);
        useEPsFlag = 0 ; % use all the coordinates of the BB and try again
        [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],[xBBAll,yBBAll],distFromBorder-1);
        numPixBBTested = length(xBBAll);
        distAll =vertcat(dist{:});     
    end
    
    
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
end % if isempty idxNeurite








%% Optional TroubleShoot Plot IV: Candidate Ridge Seeds : NEURITE ENTRANCE SEED CHOSEN
if ip.Results.TSOverlays == true
    [nyBack,nxBack] = ind2sub(size(img),find(backboneSeed)); 
    scatter(nxBack,nyBack,10,orange,'filled'); 
%     spy(backboneSeed,'r'); % plot the chosen backbone candidate
    if linked == 1 ;
        scatter(linkCoords(:,1),linkCoords(:,2),10,cDark(2,:),'filled');
    end
    scatter(xEnter,yEnter,100,'k','Marker','*'); % plot the neurite entrance point
    %linked = 0; % reset
    %figCount = figCount+1; % figure closed
end % if ip.Results.TSOverlays
%% Save the information

% Input to GCAneuriteOrientConsistencyCheck.m and GCAVeilStemReconstruct.m
backboneInfo.backboneSeedMask = backboneSeed;
backboneInfo.coordsEnterNeurite = [xEnter,yEnter];

backboneInfo.timeStamp = clock; % eventually will be saved in movieData

end % end Function

