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
%             .minSearchRadiusCG:Minimum allowed search radius (in pixels).
%             .maxSearchRadiusCG:Maximum allowed search radius (in pixels).
%                               This value is the maximum search radius
%                               between two consecutive frames as when
%                               linking between consecutive frames. It will
%                               be calcualted for different time gaps
%                               based on the scaling factor of Brownian
%                               motion (expanding it will make use of the
%                               field .timeReachConfB).
%             .brownStdMultCG : Factor multiplying Brownian
%                               displacement std to get search radius.
%                               Vector with number of entries equal to
%                               gapCloseParam.timeWindow (defined below).
%             .linStdMultCG   : Factor multiplying linear motion std to get
%                               search radius. Vector with number of entries
%                               equal to gapCloseParam.timeWindow (defined
%                               below).
%             .timeReachConfB : Time gap for reaching confinement for
%                               2D Brownian motion. For smaller time gaps,
%                               expected displacement increases with
%                               sqrt(time gap). For larger time gaps,
%                               expected displacement increases slowly with
%                               (time gap)^0.1.
%             .timeReachConfL : Time gap for reaching confinement for
%                               linear motion. Time scaling similar to
%                               timeReachConfB above.
%             .lenForClassify : Minimum length of a track to classify it as
%                               directed or Brownian.
%             .maxAngleES     : Maximum allowed angle between two
%                               directions of motion for potential linking
%                               (in degrees).
%             .maxAngleVD     : Max allowed angle between each linear track
%                               and the vector connecting the centers of
%                               the tracks (tests for colinearity)
%             .closestDistScaleCG:Scaling factor of nearest neighbor
%                                 distance.
%             .maxStdMultCG   : Maximum value of factor multiplying
%                               std to get search radius.
%             .ampRatioLimitCG: Minimum and maximum allowed ratio between
%                               the amplitude of a merged feature and the
%                               sum of the amplitude of the two features
%                               making it.
%             .lftCdf         : Lifetime cumulative density function.
%                               Column vector, specifying cdf for
%                               lifetime = 0 to movie length.
%                               Enter [] if cdf is not to be used.
%                               Optional. Default: [].
%             .useLocalDensity: 1 if local density of features is used to expand
%                               their search radius if possible, 0 otherwise.
%             .nnWindow       : Time window to be used in estimating the
%                               nearest neighbor distance of a track at its start
%                               and end.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%             .tolerance  : Relative change in number of tracks in two
%                           consecutive gap closing steps below which
%                           iteration stops.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%       kalmanFilterInfo:Structure array with number of entries equal to
%                        number of frames in movie. Contains the fields:
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


doPlot=0;

d1Max=costMatParam.d1Max;
maxFAngle = costMatParam.maxFAngle*pi/180;
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
px=trackedFeatInfo(:,1:8:end);
py=trackedFeatInfo(:,2:8:end);
vx=diff(px,1,2);
vy=diff(py,1,2);
vMag=sqrt(vx.^2+vy.^2);

trackLengths=nansum(sqrt(vx.^2+vy.^2),2);

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


if doPlot==1
    plot forward cutoff
    x=[1:tMax]';
    y=vMax*x;
    figure; plot(x,y,'r')
    hold on
    plot([x(1); x(end)],[vMax*sqrt(tMax); vMax*sqrt(tMax)],'r')
    title({'cutoffFwd (red): vMax*min(sqrt(tMax),tGap)';'cutoffBwd (blue): min(vMed*tMax,backVelMultFactor*vMax*tGap)'})

    plot backward cutoff
    y=backVelMultFactor*vMax*x;
    plot(x,y)
    plot([x(1); x(end)],[vMed*tMax; vMed*tMax])
end

tMax=gapCloseParam.timeWindow;

indx1 = zeros(10*numTracks,1);
indx2 = zeros(10*numTracks,1);
cost  = zeros(10*numTracks,4);

linkCount = 1;
for iFrame = 1:numFrames-1
    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;

    if isempty(endsToConsider)
        continue
    end

    for jFrame = iFrame + 1 : min(iFrame+tMax,numFrames)

        %find tracks that start in this frame
        startsToConsider = tracksPerFrame(jFrame).starts;

        if isempty(startsToConsider)
            continue
        end

        nStarts = length(startsToConsider);
        nEnds = length(endsToConsider);

        % time gap
        tGap = jFrame - iFrame;

        % forward and backward cutoff distances - based on velocity and the
        % current time gap
        cutoffDistFwd = vMax*min(sqrt(tMax),tGap);
        cutoffDistBwd = min(vMed*tMax,backVelMultFactor*vMax*tGap);

        % coordinates of the first point in each startsToConsider track
        spX=mat2cell(repmat(trackStartPxyVxy(startsToConsider,1),[nEnds 1]),repmat(nStarts,[nEnds 1]),1);
        spY=mat2cell(repmat(trackStartPxyVxy(startsToConsider,2),[nEnds 1]),repmat(nStarts,[nEnds 1]),1);

        [d1,d2,evYc,evXc]=arrayfun(@(i) pt2segDist([py(endsToConsider(i),:)'...
            px(endsToConsider(i),:)'],[spY{i,1},spX{i,1}],0),[1:nEnds]','UniformOutput',0);

        % for each endsToConsider track, find the distance from
        % the point along the track nearest each startsToConsider track
        % to the starting point of the startsToConsider track
        d1=cell2mat(d1')';
        % for each endsToConsider track, find the distance from the
        % track end along the lattice towards the point along the track
        % nearest each startsToConsider track
        d2=cell2mat(d2')';

        % for each endsToConsider track, find the instantaneous velocity at
        % the point along the track nearest each startsToConsider track
        % starts
        evXc=cell2mat(evXc')';
        evYc=cell2mat(evYc')';
        evMagC=sqrt(evXc.^2+evYc.^2);

        % velocity at starts of startsToConsider tracks
        svX = repmat(trackStartPxyVxy(startsToConsider,3)',[nEnds 1]);
        svY = repmat(trackStartPxyVxy(startsToConsider,4)',[nEnds 1]);
        svMag = sqrt(svX.^2 + svY.^2);

        % cos of angle between start track beginning and direction of end
        % track at closet point to start
        cosTheta = (evXc.*svX + evYc.*svY)./(evMagC.*svMag);

        % velocity at final point (f) of endsToConsider tracks
        evXf = repmat(trackEndPxyVxy(endsToConsider,3),[1 nStarts]);
        evYf = repmat(trackEndPxyVxy(endsToConsider,4),[1 nStarts]);
        evMagF = sqrt(evXf.^2 + evYf.^2);

        % displacement vector (start minus end)
        dispX = repmat(trackStartPxyVxy(startsToConsider,1)',[nEnds 1]) - repmat(trackEndPxyVxy(endsToConsider,1),[1 nStarts]);
        dispY = repmat(trackStartPxyVxy(startsToConsider,2)',[nEnds 1]) - repmat(trackEndPxyVxy(endsToConsider,2),[1 nStarts]);
        dispMag = sqrt(dispX.^2 + dispY.^2);

        % cos angle between track 1 end and track 2 start
        cosEF_SF = (evXf.*svX + evYf.*svY)./(evMagF.*svMag); % cos(alpha)
        % cos angle between track 1 end and displacement vector
        cosEF_D  = (evXf.*dispX + evYf.*dispY)./(evMagF.*dispMag); % cos(beta)

        % find candidates for forward linking
        fwdIdx=find(d1<=cutoffDistFwd & d2==0 & cosEF_D>=cos(maxFAngle) & cosTheta>=cos(maxFAngle));

        if ~isempty(fwdIdx)
            % for forward links, cosTheta==cosEF_SF and d1==dispMag
            % reassign d1 and d2 with components of displacement vector and
            % cosTheta with more accurate measurement of cos(alpha)
            d2(fwdIdx)=d1(fwdIdx).*cosEF_D(fwdIdx);
            d1(fwdIdx)=sqrt(d1(fwdIdx).^2-d2(fwdIdx).^2);
            cosTheta(fwdIdx)=cosEF_SF(fwdIdx);

            % record indices and parts of cost for forward links
            [r c]=ind2sub(size(d1),fwdIdx);
            indx1(linkCount:linkCount+length(fwdIdx)-1) = endsToConsider(r);
            indx2(linkCount:linkCount+length(fwdIdx)-1) = startsToConsider(c);
            cost(linkCount:linkCount+length(fwdIdx)-1,:) = [d1(fwdIdx) d2(fwdIdx) cosTheta(fwdIdx) ones(length(fwdIdx),1)];
            linkCount = linkCount+length(fwdIdx);
        end

        % find candidates for backward linking
        bwdIdx=find(d1<=d1Max & (d2>0 & d2<=cutoffDistBwd) & cosTheta>=cos(maxFAngle));

        if ~isempty(bwdIdx)
            % record indices and parts of cost for forward links
            [r c]=ind2sub(size(d1),bwdIdx);
            indx1(linkCount:linkCount+length(bwdIdx)-1) = endsToConsider(r);
            indx2(linkCount:linkCount+length(bwdIdx)-1) = startsToConsider(c);
            cost(linkCount:linkCount+length(bwdIdx)-1,:) = [d1(bwdIdx) d2(bwdIdx) cosTheta(bwdIdx) 2*ones(length(bwdIdx),1)];
            linkCount = linkCount+length(bwdIdx);
        end

        if doPlot==1
            figure(1);
            % plot ends
            for iTrack=1:3 %nEnds %numTracks
                figure

                currentTrack = [px(endsToConsider(iTrack),:)' py(endsToConsider(iTrack),:)'];
                currentTrack = currentTrack(trackStartTime(endsToConsider(iTrack)):trackEndTime(endsToConsider(iTrack)),:);

                plot(currentTrack(:,1),currentTrack(:,2),'b')
                hold on;
                scatter(currentTrack(:,1),currentTrack(:,2),'b.')
                quiver(trackEndPxyVxy(endsToConsider(iTrack),1),trackEndPxyVxy(endsToConsider(iTrack),2),trackEndPxyVxy(endsToConsider(iTrack),3),trackEndPxyVxy(endsToConsider(iTrack),4),'b')

                fwd=fwdIdx{iTrack};

                % plot corresponding starts from fwd
                for jTrack=1:length(fwd)
                    %get current track's coordinates
                    currentTrack = [px(startsToConsider(fwd(jTrack)),:)' py(startsToConsider(fwd(jTrack)),:)'];
                    currentTrack = currentTrack(trackStartTime(startsToConsider(fwd(jTrack))):trackEndTime(startsToConsider(fwd(jTrack))),:);

                    plot(currentTrack(:,1),currentTrack(:,2),'r')
                    hold on;
                    scatter(currentTrack(:,1),currentTrack(:,2),'r.')
                    quiver(trackEndPxyVxy(startsToConsider(fwd(jTrack)),1),trackEndPxyVxy(startsToConsider(fwd(jTrack)),2),trackEndPxyVxy(startsToConsider(fwd(jTrack)),3),trackEndPxyVxy(startsToConsider(fwd(jTrack)),4),'b')

                end
            end
            axis equal
        end
    end
end

indx1(linkCount:end) =[];
indx2(linkCount:end) =[];
cost(linkCount:end,:)=[];


type=cost(:,4);
d1Max=prctile(cost(:,1),99);
cost=cost(:,1)./d1Max + (1-cost(:,3));



% plot histograms of costs for forward and backward
if doPlot==1
    m=min(length(find(type==1)),length(find(type==2)));
    sub1=randsample(find(type==1),m);
    sub2=randsample(find(type==2),m);

    pop1=cost(sub1); t1=type(sub1);
    pop2=cost(sub2); t2=type(sub2);

    n=linspace(min([pop1;pop2]),max([pop1;pop2]),25);

    [x1,nbins1] = histc(pop1,n); %growth
    [x2,nbins2] = histc(pop2,n); %shrinkage

    figure
    bar(n,[x1 x2],'stack')
    colormap([1 0 0;0 0 1])
    legend('growth costs','shrinkage costs')
    title('d1/dmax + (1-cosTheta)')
    hold on
    plot([prctile(cost,90);prctile(cost,90)],[0,max([x1+x2])])

    % get histogram of costs for those ends which only link to one start
    [u,num]=countEntries(indx1);
    only1=u(num==1);

    costSingles=cost(cell2mat(arrayfun(@(x) find(indx1==x),only1,'uniformoutput',0)));
    [c,nbins1] = histc(costSingles,n);
    figure
    bar(n,c)
    title('costs for track ends with only one potential link')


end

% plot those tracks with more than 3 potential connections and their costs
if doPlot==1
    figure
    imagesc(.75*zeros(round(max(py(:))),round(max(px(:))))); colormap gray
    hold on

    [u,num]=countEntries(indx1);
    c=u(num>3);
    for j=1:length(c);
        b=find(indx1==c(j) | indx2==c(j));

        for i=1:length(b)
            a=find(indx1==indx1(b(i)) | indx2==indx2(b(i)));
            for k=1:length(a)
                iEnd=a(k);
                iStart=a(k);

                %get current end track's coordinates
                currentTrackE = [px(indx1(iEnd),:); py(indx1(iEnd),:)]';
                currentTrackE = currentTrackE(trackStartTime(indx1(iEnd)):trackEndTime(indx1(iEnd)),:);

                %get current start track's coordinates
                currentTrackS = [px(indx2(iStart),:); py(indx2(iStart),:)]';
                currentTrackS = currentTrackS(trackStartTime(indx2(iStart)):trackEndTime(indx2(iStart)),:);

                % plot the tracks in blue
                plot(currentTrackE(:,1),currentTrackE(:,2),'r')
                plot(currentTrackS(:,1),currentTrackS(:,2),'r')

                % plot points along tracks in red (ends) or green (starts)
                scatter(currentTrackE(:,1),currentTrackE(:,2),'b.')
                scatter(currentTrackS(:,1),currentTrackS(:,2),'g.')

                % plot possible connections
                if type(iEnd)==1
                    x=[currentTrackE(end,1);currentTrackS(1,1)];
                    y=[currentTrackE(end,2);currentTrackS(1,2)];
                    plot(x,y,'c','LineWidth',2)
                    text(mean(x),mean(y),[' \leftarrow ' sprintf('%3.2f',cost(iEnd))],'color','c');
                else
                    x=[currentTrackE(end,1);currentTrackS(1,1)];
                    y=[currentTrackE(end,2);currentTrackS(1,2)];
                    plot(x,y,'y','LineWidth',2)
                    text(mean(x),mean(y),[sprintf('%3.2f',cost(iEnd)) '\rightarrow '],'color','y','horizontalAlignment','right');
                end

                % end track end vectors
                %quiver(currentTrackE(end,1),currentTrackE(end,2),xyzVel(indx1(iEnd),1),xyzVel(indx1(iEnd),2),'r')
                quiver(trackEndPxyVxy(indx1(iEnd),1),trackEndPxyVxy(indx1(iEnd),2),trackEndPxyVxy(indx1(iEnd),3),trackEndPxyVxy(indx1(iEnd),4),'b')
                % start track end vectors
                %quiver(currentTrackS(end,1),currentTrackS(end,2),xyzVel(indx2(iStart),1),xyzVel(indx2(iStart),2),'r')
                quiver(trackEndPxyVxy(indx2(iStart),1),trackEndPxyVxy(indx2(iStart),2),trackEndPxyVxy(indx2(iStart),3),trackEndPxyVxy(indx2(iStart),4),'b')
            end
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
function [d1,d2,evY,evX]=pt2segDist(segYX,ptYX,doPlot)
% find nearest point on the directed segment B-C to point P

if doPlot==1
    figure
end

bYX=segYX(1:end-1,:); % first pt of each line segment (B)

bcYX=diff(segYX); % velocity components for C-B
BC=sqrt(sum(bcYX.^2,2)); % distance from B to C

% keep track of which entries don't exist
nanEntries=double(isnan(bcYX(:,1)));
nanEntries=swapMaskValues(nanEntries,[0,1],[1,nan]);
nPts=size(ptYX,1);
d1=zeros(nPts,1);
d2=zeros(nPts,1);
evY=zeros(nPts,1);
evX=zeros(nPts,1);
for i=1:size(ptYX,1)
    % velocity components for vectors from all points on seg to P
    temp=repmat(ptYX(i,:),[size(segYX,1),1])-segYX;

    bpYX=temp(1:end-1,:); % velocity components for P-B
    BP=sqrt(sum(bpYX.^2,2)); % distance from P to B

    cpYX=temp(2:end,:); % velocity components for P-C
    CP=sqrt(sum(cpYX.^2,2)); % distance from P to C

    % get fraction mag(vector from B to point on line closest to P)/mag(BC)
    % ref: http://www.geometrictools.com/Documentation/DistancePointLine.pdf
    t0=(bcYX(:,1).*bpYX(:,1)+bcYX(:,2).*bpYX(:,2))./(BC.^2);

    D=zeros(length(t0),1);
    extraPt=zeros(length(t0),2);

    % P falls outside segment and is closest to B
    idx=find(t0<=0);
    D(idx)=BP(idx);
    extraPt(idx,:)=segYX(idx,:);

    % P falls outside segment and is closest to C
    idx=find(t0>=0);
    D(idx)=CP(idx);
    extraPt(idx,:)=segYX(idx+1,:);

    % P falls within BC segment
    idx=find(t0>0 & t0<1);
    pYX=repmat(ptYX(i,:),[length(idx),1]);
    b_plus_t0M=bYX(idx,:)+repmat(t0(idx),[1 2]).*bcYX(idx,:);
    nearVec=pYX-b_plus_t0M;
    D(idx)=sqrt(sum(nearVec.^2,2));
    extraPt(idx,:)=b_plus_t0M;

    % don't consider where track didn't exist
    D=D.*nanEntries;

    % this is the segment index with the lowest distance
    d1Idx=find(D==nanmin(D),1,'first');
    d1(i)=D(d1Idx); % distance from P to nearest point on BC

    % add in the extra point corresponding to where the d1 vector falls on
    % the BC line
    temp=[segYX; 0 0];
    temp(d1Idx+2:end,:)=temp(d1Idx+1:end-1,:);
    temp(d1Idx+1,:)=extraPt(d1Idx,:);

    % get local segment velocity at the extra point
    velYX=nanmean(diff(temp(max(1,d1Idx-1):min(size(temp,1),d1Idx+3),:)));
    evY(i)=velYX(1);
    evX(i)=velYX(2);

    % calculate pt-to-pt displacements towards the track end and sum them
    % this is the total shrinkage distance
    % if d2=0, P is nearest the track end (no shrinkage)
    % if d2=trackLength, P is nearest the track start (complete shrinkage)
    d2(i)=nansum(sqrt(sum(diff(temp(d1Idx+1:end,:)).^2,2)));

    if doPlot==1
        plot(segYX(:,2),segYX(:,1))
        hold on;
        scatter(temp(:,2),temp(:,1),'b.')
        scatter(ptYX(i,2),ptYX(i,1),'r.')
        text(ptYX(i,2),ptYX(i,1),['\leftarrow ' sprintf('%3.2f',d1(i)) ', ' sprintf('%3.2f',d2(i))])
        scatter(extraPt(d1Idx,2),extraPt(d1Idx,1),'g')
        quiver(extraPt(d1Idx,2),extraPt(d1Idx,1),vYX(2)+.1,vYX(1)+.1,0,'r')
        plot([ptYX(i,2); extraPt(d1Idx,2)],[ptYX(i,1); extraPt(d1Idx,1)],'g')
        axis equal
    end
end