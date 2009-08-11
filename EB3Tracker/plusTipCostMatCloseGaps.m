function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
    errFlag] = plusTipCostMatCloseGaps(trackedFeatInfo,...
    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%COSTMATLINEARMOTIONCLOSEGAPS provides a cost matrix for closing gaps and capturing merges/splits using Kalman filter information
%
%SYNOPSIS [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
%    errFlag] = costMatLinearMotionCloseGaps(trackedFeatInfo,...
%    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
%    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesKalman.
%                        Number of rows = number of tracks.
%                        Number of columns = 8*number of frames.
%                        Each row consists of
%                        [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points
%                        where the track does not exist.
%       trackedFeatIndx: Connectivity matrix of features between frames.
%                        Rows indicate continuous tracks, while columns
%                        indicate frames. A track that ends before the
%                        last time point is followed by zeros, and a track
%                        that starts at a time after the first time point
%                        is preceded by zeros.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParam   : Structure containing variables needed for cost
%                        calculation. Contains the fields:
%          .fluctRad             : Size in pixels of tube radius around track
%                                trajectory. The search region in the backward
%                                direction will expand out from the track at
%                                the final point.  This value also determines
%                                the search radius around the final point
%                                wherein any candidates will be considered for 
%                                forward linking, even if they fall slightly 
%                                behind the point.  This ensures that tracks
%                                starting from a fluctuation during a pause
%                                will still be picked up as candidates for
%                                pause.
%          .maxFAngle          : Max angle in degrees allowed between the end 
%                                track's final velocity vector and the
%                                displacement vector between end and start.
%                                Also the max angle between the end and start 
%                                tracks themselves.
%          .maxBAngle          : Angle in degrees used to expand backward
%                                search region, giving a distance-dependent
%                                criterion for how far a start track
%                                can be from the lattice to be considered a
%                                candidate for linking.
%          .backVelMultFactor : Muliplication factor of max growth speed used
%                               to define candidate search area in the
%                               backward direction.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to
%                           it.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%                           For MT tracking, there are no merges/splits, so
%                           this should be 0.
%       kalmanFilterInfo: Structure array with number of entries equal to
%                         number of frames in movie. Contains the fields:
%             .stateVec   : Kalman filter state vector for each
%                           feature in frame.
%             .stateCov   : Kalman filter state covariance matrix
%                           for each feature in frame.
%             .noiseVar   : Variance of state noise for each
%                           feature in frame.
%             .stateNoise : Estimated state noise for each feature in
%                           frame.
%             .scheme     : 1st column: propagation scheme connecting
%                           feature to previous feature. 2nd column:
%                           propagation scheme connecting feature to
%                           next feature.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       movieInfo      : movieInfo as input to trackCloseGapsKalman. Not
%                        really used in this code, but needed for
%                        compatibility with other cost functions.
%
%OUTPUT costMat       : Cost matrix.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%REMARKS
%
%Khuloud Jaqaman, April 2007

%% Output

costMat = [];
nonlinkMarker = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('plusTipCostMatCloseGaps')
    disp('--plusTipCostMatCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

doTest=1
doPlot=0;

fluctRad=costMatParam.fluctRad;
maxFAngle = costMatParam.maxFAngle*pi/180;
costMatParam.maxBAngle=10;
maxBAngle = costMatParam.maxBAngle*pi/180;
backVelMultFactor = costMatParam.backVelMultFactor;


%find the number of tracks to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
for iFrame = 1 : numFrames
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Gap closing

% extract feature positions and velocity components
px=trackedFeatInfo(:,1:8:end);
py=trackedFeatInfo(:,2:8:end);
vx=diff(px,1,2);
vy=diff(py,1,2);
%vMag=sqrt(vx.^2+vy.^2);

%trackLengths=nansum(sqrt(vx.^2+vy.^2),2);

% ax=vx(:,1:end-1); [row col]=size(ax); ax=ax(:)';
% ay=vy(:,1:end-1); ay=ay(:)';
%
% bx=vx(:,2:end); bx=bx(:)';
% by=vy(:,2:end); by=by(:)';
%
% c=cross([ax;ay;zeros(1,length(ax))],[bx;by;zeros(1,length(bx))]);
% magC=sqrt(c(1,:).^2+c(2,:).^2+c(3,:).^2);
% magC=reshape(magC,row,col);



% TRACK STARTS
trackStartPxyVxy = zeros(numTracks,4);
% x and y coordinates of the track's first point
trackStartPxyVxy(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'first')),[1:numTracks]','UniformOutput',0));
trackStartPxyVxy(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'first')),[1:numTracks]','UniformOutput',0));
% average of first three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackStartPxyVxy(:,3)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'first'))),[1:numTracks]','UniformOutput',0));
trackStartPxyVxy(:,4)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'first'))),[1:numTracks]','UniformOutput',0));

% TRACK ENDS
trackEndPxyVxy = zeros(numTracks,4);
% x and y coordinates of the track's last point
trackEndPxyVxy(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'last')),[1:numTracks]','UniformOutput',0));
trackEndPxyVxy(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'last')),[1:numTracks]','UniformOutput',0));
% average of last three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackEndPxyVxy(:,3)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'last'))),[1:numTracks]','UniformOutput',0));
trackEndPxyVxy(:,4)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'last'))),[1:numTracks]','UniformOutput',0));

% get velocity components from kalman filter (very similar to trackEndVxy)
xyzVel=cell2mat(arrayfun(@(iTrack) kalmanFilterInfo(trackEndTime(iTrack))...
    .stateVec(trackedFeatIndx(iTrack,trackEndTime(iTrack)),probDim+1:2*probDim),...
    [1:numTracks]','UniformOutput',0));


trackEndSpeed=sqrt(sum(xyzVel.^2,2));
vMax=prctile(trackEndSpeed,95);
vMed=median(trackEndSpeed);

tMax=gapCloseParam.timeWindow;

indx1 = zeros(10*numTracks,1);
indx2 = zeros(10*numTracks,1);
costComponents  = zeros(10*numTracks,5);

% get start and end frames for each track
sFrameAll=zeros(numTracks,1);
eFrameAll=zeros(numTracks,1);
for iFrame=1:numFrames
    sFrameAll(tracksPerFrame(iFrame).starts)=iFrame;
    eFrameAll(tracksPerFrame(iFrame).ends)=iFrame;
end

linkCount = 1;
for iFrame = 1:numFrames-1
    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;

    if isempty(endsToConsider)
        continue
    end

    % these are the frames to consider for finding starts
    jFrame=iFrame+1:min(iFrame+tMax,numFrames);

    % find tracks that start in possible frame range
    startsToConsider=arrayfun(@(x) tracksPerFrame(x).starts,jFrame,'uniformoutput',0);
    % get number of starts in each frame
    nStarts=cellfun(@(x) length(x),startsToConsider,'uniformoutput',0);

    if isempty(startsToConsider)
        continue
    end

    % n-vector of time gaps allowable
    tGap = jFrame - iFrame;

    % forward and backward cutoff distances - based on velocity and the
    % current time gap
    cutoffDistFwd = vMax*min(sqrt(tMax),tGap);
    cutoffDistBwd = min(vMed*tMax,backVelMultFactor*vMax*tGap);

    % plot the cutoff distances as a function of time gap
    if doPlot==1
        figure; plot(cutoffDistFwd); hold on; plot(cutoffDistBwd,'r');
    end

    % make vectors containing forward/backward cutoffs for each start
    cutFwdPerVec=cell2mat(arrayfun(@(i,j) repmat(i,[j,1]),cutoffDistFwd,cell2mat(nStarts),'uniformoutput',0)');
    cutBwdPerVec=cell2mat(arrayfun(@(i,j) repmat(i,[j,1]),cutoffDistBwd,cell2mat(nStarts),'uniformoutput',0)');

    startsToConsider=cell2mat(startsToConsider');

    nStarts = length(startsToConsider);
    nEnds = length(endsToConsider);

    % end track last pt
    epX=repmat(trackEndPxyVxy(endsToConsider,1),[1 nStarts]);
    epY=repmat(trackEndPxyVxy(endsToConsider,2),[1 nStarts]);
    % start track first pt
    spX=repmat(trackStartPxyVxy(startsToConsider,1)',[nEnds 1]);
    spY=repmat(trackStartPxyVxy(startsToConsider,2)',[nEnds 1]);
    % nEnds x nStarts distance matrix containing displacement vector
    % magnittude
    dispMag=sqrt((epX-spX).^2+(epY-spY).^2);

    % we will only consider end/start pairs where the distance from end to start is less than
    % the max(forwardCutoff,backwardCutoff)
    maxCut=max(repmat(cutFwdPerVec',[nEnds 1]),repmat(cutBwdPerVec',[nEnds 1]));

    dPerp=zeros(nStarts*nEnds,1);
    dPara=zeros(nStarts*nEnds,1);
    evYc=zeros(nStarts*nEnds,1);
    evXc=zeros(nStarts*nEnds,1);
    endLinkIdx=zeros(nStarts*nEnds,1);
    startLinkIdx=zeros(nStarts*nEnds,1);
    sAll=zeros(nStarts*nEnds,1);

    %%%%%%%%%%%%%%% load b4gapclosing from NIH_Waterman
    % need to find direction of velocity in bw direction and backproject by
    % large amt (maxCut?)
    % plot this to check, then see how this affects cand selection - is
    % criterion based on whether pt is trackLength-distance away?  if so, need
    % to redefine b/c this will change when we add a phantom point...


    count=1;
    if doTest==1
        endRange=1:5;
    else
        endRange=1:nEnds;
    end


    for iEnd=endRange
        % indices correspoinding to current set of starts and ends,
        % respectively
        startCandidateIdx=find(dispMag(iEnd,:)<maxCut(iEnd,:))';
        endCandidateIdx=repmat(iEnd,[length(startCandidateIdx) 1]);

        if isempty(startCandidateIdx)
            continue
        end

        % coordinates of the first point in each startsToConsider track
        sX=trackStartPxyVxy(startsToConsider(startCandidateIdx),1);
        sY=trackStartPxyVxy(startsToConsider(startCandidateIdx),2);

        if doTest==1
            %plot the track and the potential starts...
            % figure
            xtemp=px(endsToConsider(iEnd),:)'; ytemp=py(endsToConsider(iEnd),:)';
            xtemp(isnan(xtemp))=[]; ytemp(isnan(ytemp))=[];
            
            
           


            % add in some fake pts...
%             r=150; nPts=10000;
%             sX=r.*rand(nPts,1)+xtemp(end)-r/2;
%             sY=r.*rand(nPts,1)+ytemp(end)-r/2;
%             cutFwdPerVec=randsample(cutoffDistFwd,nPts,'true')';
%             cutBwdPerVec=randsample(cutoffDistBwd,nPts,'true')';

            [sX sY]= meshgrid([round(xtemp(end))-100:round(xtemp(end))+100],[round(ytemp(end))-100:round(ytemp(end))+100]);
            sX=sX(:);
            sY=sY(:);
            negIdx=find(sX<=0 | sY<=0);
            sX(negIdx)=[];
            sY(negIdx)=[];
            
            img=zeros(max(sY),max(sX));
            
            
            cutFwdPerVec=repmat(max(cutoffDistFwd),size(sX));
            cutBwdPerVec=repmat(max(cutoffDistBwd),size(sX));
            
            %scatter(sX(:),sY(:))
            
            
            %scatter(sX,sY,'r.')

            
        end

        % call subfunction to calculate magnitude of the components of
        % the vector pointing from end to start, as well as the components
        % of the local velocity along the end track at its closest point to
        % the start track
        [dPerpTemp,dParaTemp,evYcTemp,evXcTemp]=pt2segDist([py(endsToConsider(iEnd),:)',px(endsToConsider(iEnd),:)'],trackStartPxyVxy(endsToConsider(iEnd),:),[sY,sX],cutoffDistBwd,0);

        % dPerp is the component perpendicular to the end track
        dPerp(count:count+length(sX)-1)=dPerpTemp;
        % dPara is the component parallel to the end track
        dPara(count:count+length(sX)-1)=dParaTemp;

        % evX/Yc are the velocity components of the end track at
        % the point closest (c) each startsToConsider track starts
        evXc(count:count+length(sX)-1)=evXcTemp;
        evYc(count:count+length(sX)-1)=evYcTemp;



        if doTest==1
            % redefine these based on test points
            endCandidateIdx=repmat(iEnd,[length(sX) 1]);
            endLinkIdx  (count:count+length(sX)-1) = endsToConsider(endCandidateIdx);
            startLinkIdx(count:count+length(sX)-1)=1:length(sX);
            sAll(count:count+length(sX)-1)=1:length(sX);
        else
            % starts/endsToConsider indices of the ones checked
            endLinkIdx  (count:count+length(sX)-1) =   endsToConsider(endCandidateIdx);
            startLinkIdx(count:count+length(sX)-1) = startsToConsider(startCandidateIdx);
            sAll(count:count+length(sX)-1)=startCandidateIdx;
        end



        count=count+length(sX);
    end
    if count==1
        continue
    end

    dPerp(count:end)=[];
    dPara(count:end)=[];
    evYc(count:end)=[];
    evXc(count:end)=[];
    endLinkIdx(count:end)=[];
    startLinkIdx(count:end)=[];
    sAll(count:end)=[];

    % velocity at starts of startsToConsider tracks
    if doTest==1
        %svX = 1*rand(length(sX),1).*randsample([-1;1],length(sX),'true');
        %svY = 1*rand(length(sX),1).*randsample([-1;1],length(sX),'true');
        svX = repmat(trackEndPxyVxy(endLinkIdx(1),3),size(sX));
        svY = repmat(trackEndPxyVxy(endLinkIdx(1),4),size(sX));
    else
        svX = trackStartPxyVxy(startLinkIdx,3);
        svY = trackStartPxyVxy(startLinkIdx,4);
    end
    svMag = sqrt(svX.^2 + svY.^2);

    % cos of angle between start track beginning and direction of end
    % track at closest point to start
    evMagC=sqrt(evXc.^2+evYc.^2);
    cosTheta = (evXc.*svX + evYc.*svY)./(evMagC.*svMag);

    % velocity at final point (f) of endsToConsider tracks
    evXf = trackEndPxyVxy(endLinkIdx,3);
    evYf = trackEndPxyVxy(endLinkIdx,4);
    evMagF = sqrt(evXf.^2 + evYf.^2);


    % displacement vector (start minus end)
    if doTest==1
        dispX = sX-trackEndPxyVxy(endLinkIdx,1);
        dispY = sY-trackEndPxyVxy(endLinkIdx,2);
    else
        dispX = trackStartPxyVxy(startLinkIdx,1)-trackEndPxyVxy(endLinkIdx,1);
        dispY = trackStartPxyVxy(startLinkIdx,2)-trackEndPxyVxy(endLinkIdx,2);
    end
    dispMag = sqrt(dispX.^2 + dispY.^2);

    % cos angle between track 1 end and track 2 start
    cosEF_SF = (evXf.*svX + evYf.*svY)./(evMagF.*svMag); % cos(alpha)

    % cos angle between track 1 end and displacement vector
    cosEF_D  = (evXf.*dispX + evYf.*dispY)./(evMagF.*dispMag); % cos(beta)

    % criteria for backward linking:
    % perp dist (dPerp) must be smaller than user-set fluctRad
    % nearest pt needs to not be the end of the endTrack and parallel dist
    % should be smaller than backward cutoff
    % angle between tracks should be less than max forward angle
    bwdIdx=find(dPerp<=(fluctRad+dPara*tan(maxBAngle)) & (dPara>fluctRad & dPara<=cutBwdPerVec(sAll)) & cosTheta>=cos(maxFAngle));

    bwdIdx1=find(dPerp<=(fluctRad) & (dPara>fluctRad & dPara<=cutBwdPerVec(sAll)) & cosTheta>=cos(maxFAngle));
    
    if ~isempty(bwdIdx)
        % record indices and parts of cost for forward links
        indx1(linkCount:linkCount+length(bwdIdx)-1) = endLinkIdx(bwdIdx);
        indx2(linkCount:linkCount+length(bwdIdx)-1) = startLinkIdx(bwdIdx);

        % cost - keep several pieces of data here for now
        % [dPerp dPara cosTheta 2 (for backward)]
        costComponents(linkCount:linkCount+length(bwdIdx)-1,1:4) = [dPerp(bwdIdx) dPara(bwdIdx) cosTheta(bwdIdx) 2*ones(length(bwdIdx),1)];
        linkCount = linkCount+length(bwdIdx);
    end

    % criteria for forward linking:
    % parallel dist (dPara) must be 0 (indicates closest pt is the end pt)
    % end-start dist must be smaller than forward cutoff
    % end-displacement angle must be smaller than max forward angle
    % angle between tracks should be less than max forward angle 
    fwdIdx1=find(dPerp<=cutFwdPerVec(sAll) & dPara==0 & cosEF_D>=cos(maxFAngle) & cosTheta>=cos(maxFAngle));
    % but also count those tracks falling within fluctRad of the track
    % end in the backward direction as forward links
    fwdIdx2=setdiff(find(dispMag<=fluctRad),fwdIdx1);
    
    if ~isempty(unique([fwdIdx1; fwdIdx2]))
        % for forward links, currently cosTheta=cosEF_SF and dPerp=dispMag
        % reassign dPerp and dPara with components of displacement vector and
        % cosTheta with more accurate measurement of cos(alpha)
        dPara(fwdIdx1)=dPerp(fwdIdx1).*cosEF_D(fwdIdx1);
        dPerp(fwdIdx1)=sqrt(dPerp(fwdIdx1).^2-dPara(fwdIdx1).^2);
        cosTheta(fwdIdx1)=cosEF_SF(fwdIdx1);
        
        fwdIdx=[fwdIdx1; fwdIdx2];

        % record indices and parts of cost for forward links
        indx1(linkCount:linkCount+length(fwdIdx)-1) = endLinkIdx(fwdIdx);
        indx2(linkCount:linkCount+length(fwdIdx)-1) = startLinkIdx(fwdIdx);

        % cost - keep several pieces of data here for now
        % [dPerp dPara cosTheta 1 (for forward)]
        costComponents(linkCount:linkCount+length(fwdIdx)-1,1:4) = [dPerp(fwdIdx) dPara(fwdIdx) cosTheta(fwdIdx) ones(length(fwdIdx),1)];
        linkCount = linkCount+length(fwdIdx);
    end
    
    if doTest==1
        
        img(sub2ind(size(img),sY(fwdIdx),sX(fwdIdx)))=1;
        img(sub2ind(size(img),sY(bwdIdx),sX(bwdIdx)))=2;
        img(sub2ind(size(img),sY(bwdIdx1),sX(bwdIdx1)))=3;

        
        figure(gcf)
        
            hold on
            imagesc(img); 
            plot(xtemp,ytemp,'linewidth',2)            
            axis equal
            scatter(xtemp,ytemp,'b.')
            scatter(xtemp(end),ytemp(end),'b*')
        % show velocity vectors for those simulated start points not
        % selected as either forward or backward links
        %quiver(sX(setdiff(sAll,[fwdIdx; bwdIdx])),sY(setdiff(sAll,[fwdIdx; bwdIdx])),...
        %   svX(setdiff(sAll,[fwdIdx; bwdIdx])),svY(setdiff(sAll,[fwdIdx; bwdIdx])),'k')
        
%         % green are all backward links
%         quiver(sX(bwdIdx),sY(bwdIdx),svX(bwdIdx),svY(bwdIdx),'g')
%         % red are backward links within sausage only - not used right now
%         quiver(sX(bwdIdx1),sY(bwdIdx1),svX(bwdIdx1),svY(bwdIdx1),'r')
%         % blue are forward links
%         quiver(sX(fwdIdx),sY(fwdIdx),svX(fwdIdx),svY(fwdIdx),'b')
        
        %plusTipSearchRadPlot(trackEndPxyVxy(endLinkIdx(1),:),cutoffDistFwd,maxFAngle)
        %plusTipSearchRadPlot([1 1 -1 -1].*trackEndPxyVxy(endLinkIdx(1),:),cutoffDistBwd,maxBAngle)
    end
end


% figure 
% hold on
% for i=5:5:20
% plot(1:20,1+tan(i*pi/180).*[1:20])
% plot(1:20,-(1+tan(i*pi/180).*[1:20]))
% end

indx1(linkCount:end) =[];
indx2(linkCount:end) =[];
costComponents(linkCount:end,:)=[];
costComponents(:,5)=sFrameAll(indx2)-eFrameAll(indx1);

% type is 1 for forward, 2 for backward
type=costComponents(:,4);

% calculate the cost
d1NormFactor=prctile(costComponents(:,1),99);
cost=1.5.^costComponents(:,5).*(costComponents(:,1)./d1NormFactor + (1-costComponents(:,3)));


% plot histograms of costs for forward and backward
doPlot=0;
if doPlot==1
    % to make a stacked plot we need equal sample sizes for both forward and
    % backward populations.  here we find the max sample size which is the
    % min of nForwardEvents or nBackwardEvents
    m=min(length(find(type==1)),length(find(type==2)));
    sub1=randsample(find(type==1),m); % forward indices
    sub2=randsample(find(type==2),m); % backward indices

    pop1=cost(sub1); % sampled forward costs
    pop2=cost(sub2); % sampled backward costs

    % create x-axis bins spanning all costs in sample
    n=linspace(min([pop1;pop2]),max([pop1;pop2]),25);

    % bin the samples
    [x1,nbins1] = histc(pop1,n); % forward
    [x2,nbins2] = histc(pop2,n); % backward

    % make the plot
    subplot(2,1,1)
    bar(n,[x1 x2],'stack')
    colormap([1 0 0;0 0 1])
    legend('Forward Costs','Shrinkage Costs')
    %title('Cost for all tracks')
    hold on
    deathCost=prctile(cost,90);
    plot([deathCost;deathCost],[0,max([x1+x2])])

    % get histogram of costs for those ends which only link to one start
    [u,num]=countEntries(indx1);
    only1=u(num==1); % indices from endsToConsider with only 1 potential start link

    costSingles=cost(cell2mat(arrayfun(@(x) find(indx1==x),only1,'uniformoutput',0)));
    typeSingles=type(cell2mat(arrayfun(@(x) find(indx1==x),only1,'uniformoutput',0)));

    m=min(length(find(typeSingles==1)),length(find(typeSingles==2)));
    sub1=randsample(find(typeSingles==1),m); % forward indices
    sub2=randsample(find(typeSingles==2),m); % backward indices

    pop1=costSingles(sub1); % sampled forward costs
    pop2=costSingles(sub2); % sampled backward costs

    % bin the samples
    [x1,nbins1] = histc(pop1,n); % forward
    [x2,nbins2] = histc(pop2,n); % backward

    % make the plot
    subplot(2,1,2)
    bar(n,[x1 x2],'stack')
    colormap([1 0 0;0 0 1])
    %     legend('Forward Costs','Shrinkage Costs')
    %     title('Costs for track ends with only one potential link')
    hold on
    deathCost=prctile(cost,90);
    plot([deathCost;deathCost],[0,max([x1+x2])])

end

% plot those tracks with more than 3 potential connections and their costs
if doPlot==1
    figure
    imagesc(.75*zeros(round(max(py(:))),round(max(px(:))))); colormap gray
    hold on

    [u,num]=countEntries(indx1);
    c=u(num>3);
    for j=1:length(c);
        b=find(indx1==c(j));
        [indx1(b) indx2(b) type(b) cost(b)]
        for i=1:length(b)
            idx=b(i);

            %get current end track's coordinates
            currentTrackE = [px(indx1(idx),:); py(indx1(idx),:)]';
            currentTrackE = currentTrackE(trackStartTime(indx1(idx)):trackEndTime(indx1(idx)),:);

            %get current start track's coordinates
            currentTrackS = [px(indx2(idx),:); py(indx2(idx),:)]';
            currentTrackS = currentTrackS(trackStartTime(indx2(idx)):trackEndTime(indx2(idx)),:);

            % plot the tracks in blue
            plot(currentTrackE(:,1),currentTrackE(:,2),'g')
            plot(currentTrackS(:,1),currentTrackS(:,2),'r')

            % plot points along tracks in red (ends) or green (starts)
            scatter(currentTrackE(:,1),currentTrackE(:,2),'b.')
            scatter(currentTrackS(:,1),currentTrackS(:,2),'b.')

            % plot possible connections
            if type(iEnd)==1
                x=[currentTrackE(end,1);currentTrackS(1,1)];
                y=[currentTrackE(end,2);currentTrackS(1,2)];
                plot(x,y,'c','LineWidth',2)
                text(mean(x),mean(y),[' \leftarrow ' sprintf('%3.2f',cost(idx))],'color','c');
            else
                x=[currentTrackE(end,1);currentTrackS(1,1)];
                y=[currentTrackE(end,2);currentTrackS(1,2)];
                plot(x,y,'y','LineWidth',2)
                text(mean(x),mean(y),[sprintf('%3.2f',cost(idx)) '\rightarrow '],'color','y','horizontalAlignment','right');
            end

            % end track end vectors
            %quiver(currentTrackE(end,1),currentTrackE(end,2),xyzVel(indx1(iEnd),1),xyzVel(indx1(iEnd),2),'r')
            quiver(trackEndPxyVxy(indx1(idx),1),trackEndPxyVxy(indx1(idx),2),trackEndPxyVxy(indx1(idx),3),trackEndPxyVxy(indx1(idx),4),'g')
            % start track end vectors
            %quiver(currentTrackS(end,1),currentTrackS(end,2),xyzVel(indx2(iStart),1),xyzVel(indx2(iStart),2),'r')
            quiver(trackEndPxyVxy(indx2(idx),1),trackEndPxyVxy(indx2(idx),2),trackEndPxyVxy(indx2(idx),3),trackEndPxyVxy(indx2(idx),4),'b')

        end

    end

    axis equal
end


%% Merging and splitting

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%create cost matrix without births and deaths
numEndSplit = numTracks;
numStartMerge = numTracks;
costMat = sparse(indx1,indx2,cost,numEndSplit,numStartMerge);

%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
costBD = prctile(cost,90);

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);

% % %create cost matrix that allows for births and deaths
costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags([costBD*ones(numTracks,1); altCostSplit],0,numEndSplit,numEndSplit); ... %costs for death
    spdiags([costBD*ones(numTracks,1); altCostMerge],0,numStartMerge,numStartMerge) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);


%% ~~~ the end ~~~


function [dPerp,dPara,evY,evX]=pt2segDist(segYX,endTrackStartPxyVxy,ptYX,maxCutoff,doPlot)
% find nearest point on the directed segment B-C to point P
% see http://www.geometrictools.com/Documentation/DistancePointLine.pdf for
% inspiration

% the point of this subfunction is to find the perpendicular and parallel
% distances of one or more points (stored in ptYX), which represent a
% candidate track's starting position, to a line segment (segYX), which
% represents the end track itself.  we also want to know the
% velocity components of the point on segment nearest the
% track start.  because the MT tarck may not be straight, dPara, the
% component parallel to the track, is here the actual distance along the
% lattice, not the euclidean distance from track end to the candidate start.
% dPerp is perpendicular to the line segment unless of course
% the point of interest is nearest one of the two segment ends; then it is
% defined as the Euclidean distance.

% to extend the track in the negative direction, I suggest adding more
% points here to segYX at the start....

if doPlot==1
    figure
end

segYXorig=segYX;

% find how many phantom pts are needed to extend bwd vector past max bwd
% cutoff
endMag=sqrt(sum(endTrackStartPxyVxy(1,3:4).^2));
nRepPhantom=ceil(max(maxCutoff)/endMag);

phantomYX=repmat(endTrackStartPxyVxy(1,2:-1:1),[nRepPhantom,1])-repmat([nRepPhantom:-1:1]',[1,2]).*repmat(endTrackStartPxyVxy(1,4:-1:3),[nRepPhantom,1]);
segYX=[phantomYX; segYX];

% treat every consecutive pair of points in segYX as a line segment BC
bYX=segYX(1:end-1,:); % here's a vector containing the first pt of each line segment (B)
bcYX=diff(segYX); % velocity components
BC=sqrt(sum(bcYX.^2,2)); % BC length

% keep track of which entries don't exist
nanEntries=double(isnan(bcYX(:,1)));
nanEntries=swapMaskValues(nanEntries,[0,1],[1,nan]);
nPts=size(ptYX,1);
dPerp=zeros(nPts,1);
dPara=zeros(nPts,1);
evY=zeros(nPts,1);
evX=zeros(nPts,1);
for i=1:size(ptYX,1)
    % velocity components for vectors from all points on seg BC to P
    temp=repmat(ptYX(i,:),[size(segYX,1),1])-segYX;

    bpYX=temp(1:end-1,:); % velocity components for P-B
    BP=sqrt(sum(bpYX.^2,2)); % distance from P to B

    cpYX=temp(2:end,:); % velocity components for P-C
    CP=sqrt(sum(cpYX.^2,2)); % distance from P to C

    % get fraction mag(vector from B to point on line closest to P)/mag(BC)
    % if t0<0, closest to b; if t0>0 then closet to c
    t0=(bcYX(:,1).*bpYX(:,1)+bcYX(:,2).*bpYX(:,2))./(BC.^2);

    D=zeros(length(t0),1);
    extraPt=zeros(length(t0),2);

    % P falls outside segment and is closest to B
    idx=find(t0<=0);
    D(idx)=BP(idx);
    extraPt(idx,:)=segYX(idx,:); % just duplicate the first point

    % P falls outside segment and is closest to C
    idx=find(t0>=1);
    D(idx)=CP(idx);
    extraPt(idx,:)=segYX(idx+1,:); % duplicate the last point

    % P falls within BC segment
    idx=find(t0>0 & t0<1);
    pYX=repmat(ptYX(i,:),[length(idx),1]);
    b_plus_t0M=bYX(idx,:)+repmat(t0(idx),[1 2]).*bcYX(idx,:); % location of perp point along segment
    nearVec=pYX-b_plus_t0M;
    D(idx)=sqrt(sum(nearVec.^2,2));
    extraPt(idx,:)=b_plus_t0M; % we'll have to insert this point into the list of pts in segYX

    % don't consider where track didn't exist
    D=D.*nanEntries;

    % this is the segment index with the lowest distance
    d1Idx=find(D==nanmin(D),1,'first'); % this is where we will insert the extra point
    dPerp(i)=D(d1Idx); % distance from P to nearest point on BC

    % add in the extra point corresponding to where the dPerp vector falls on
    % the BC line
    temp=[segYX; 0 0];
    temp(d1Idx+2:end,:)=temp(d1Idx+1:end-1,:);
    temp(d1Idx+1,:)=extraPt(d1Idx,:);

    % get local segment velocity at the extra point. do this by finding the
    % pt behind up and the three pts ahead of the extra pt (closest pt on
    % the line to the start pt of interest). since there is a duplicated
    % pt, this could throw off the velocity calculation if we just did a
    % mean of the differences between pts. instead, sum the differences and
    % divide by the number of *segments* (nPts-1), minus 1 for the zero
    % difference from one point to its duplicate. this should give a good
    % local estimate of the instantaneous velocity
    pts2getLocalVel = temp(max(1,d1Idx-1):min(size(temp,1),d1Idx+3),:);
    pts2getLocalVel(isnan(pts2getLocalVel(:,1)),:)=[];
    velYX=sum(diff(pts2getLocalVel))./(size(pts2getLocalVel,1)-2);
    evY(i)=velYX(1);
    evX(i)=velYX(2);

    % calculate pt-to-pt displacements towards the track end and sum them
    % this is the total shrinkage distance
    % if dPara=0, P is nearest the track end (no shrinkage)
    % if dPara=trackLength, P is nearest the track start (complete shrinkage)
    dPara(i)=nansum(sqrt(sum(diff(temp(d1Idx+1:end,:)).^2,2)));

    if doPlot==1
        % plot end track as a blue line
        plot(segYX(:,2),segYX(:,1))
        hold on;
        % add blue dots for detection events
        scatter(temp(:,2),temp(:,1),'b.')
        % add magenta line showing start velocity vector
        quiver(endTrackStartPxyVxy(1),endTrackStartPxyVxy(2),endTrackStartPxyVxy(3),endTrackStartPxyVxy(4),0,'m')
        % plot candidate start point in red
        scatter(ptYX(i,2),ptYX(i,1),'r.')
        % show dPerp,dPara distances next to point
        text(ptYX(i,2),ptYX(i,1),['\leftarrow ' sprintf('%3.2f',dPerp(i)) ', ' sprintf('%3.2f',dPara(i))])
        % show newly-created extra pt as green circle
        scatter(extraPt(d1Idx,2),extraPt(d1Idx,1),'g')
        % show end-track calculated instantaneous velocity at new pt
        quiver(extraPt(d1Idx,2),extraPt(d1Idx,1),velYX(2)+.1,velYX(1)+.1,0,'r')
        % connect new pt to candidate start pt with green line
        plot([ptYX(i,2); extraPt(d1Idx,2)],[ptYX(i,1); extraPt(d1Idx,1)],'g')
        axis equal
    end
end