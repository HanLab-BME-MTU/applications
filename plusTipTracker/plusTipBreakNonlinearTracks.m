function [newTrackedFeatureInfo,newTrackedFeatureIndx,newNnDistFeatures]=plusTipBreakNonlinearTracks(trackedFeatureInfo,trackedFeatureIndx,nnDistFeatures)
% plusTipBreakNonlinearTracks splits up tracks not following uni-directional behavior

% this function splits up individial tracks where the displacement vectors
% created by consecutive frame-frame pairs in a track differ in direction
% by more than 45 degrees. for each instance like this, one new track is
% born. 
%
% some instances where the angle is greater than 45 degrees are retained if
% one or both of the vectors is very short (<3rd percentile), which may
% simply reflect the uncertainty in the detected position if the features
% are not moving much per frame.
%
% output has the same format as input. NOTE: nnDistFeatures will now be
% nonsensical, but these are currently not used in the gap closing cost
% matrix.  many links will now shorter than minLength, but those can be
% filtered out by post-processing step.

% this function is currently called in conditional statement during
% trackCloseGapsKalman


%get total number of tracks
[nTracks, nFrames] = size(trackedFeatureIndx);

% find track start and end frames
trackSEL = getTrackSEL(trackedFeatureInfo);

% extract coordinates
px=trackedFeatureInfo(:,1:8:end);
py=trackedFeatureInfo(:,2:8:end);

% get vector coordinates between linked features
vx=diff(px,1,2);
vy=diff(py,1,2);
vmag=sqrt(vx.^2+vy.^2);

% first vector matrix
v1x=vx(:,1:end-1);
v1y=vy(:,1:end-1);
v1mag=sqrt(v1x.^2+v1y.^2);

% second vector matrix
v2x=vx(:,2:end);
v2y=vy(:,2:end);
v2mag=sqrt(v2x.^2+v2y.^2);

% cos of angle between consecutive vectors (displacements)
cosV12=(v1x.*v2x+v1y.*v2y)./(v1mag.*v2mag);

% assume max angle is 45 degrees
cosMax=cos(45*pi/180);

% lower bound displacement - if smaller than this, may just be jitter
lb=prctile(vmag(:),3); 

% keep track of where cos or displacement vectors are NaNs
nanMat=swapMaskValues(isnan(cosV12) | isnan(v1mag) | isnan(v2mag),[0 1],[1 NaN]);

% these are within forward cone, so they're ok
okAngles=cosV12>cosMax;

% these are not in the forward cone but one or both of the vectors is shorter 
% than the nth percentile of all vectors, so maybe just jitter
okLength=cosV12<cosMax & (v1mag<lb | v2mag<lb); 

% if either the angle or the length criterion isn't met, then it's a bad
% link which we will break
badLinks=swapMaskValues(nanMat.*(okAngles | okLength),[0 1],[1 0]);


[badTrackIdx badTrackVecPair]=find(badLinks==1);
badTrackVecHead=badTrackVecPair+1;

doPlot=0;
if doPlot==1
    % plot the first 50 bad tracks
    b=badTrackIdx(1:50);
    figure
    plot(px(b,:)',py(b,:)')
    hold on
    x=px(b,:)'; x=x(:);
    y=py(b,:)'; y=y(:);
    scatter(x,y,'b')
    % plot break points in red
    for i=1:length(b)
        a=badTrackIdx(badTrackIdx==b(i));
        d=badTrackVecPair(badTrackIdx==b(i));
        x=px(sub2ind(size(badLinks),a,d)+length(px));
        y=py(sub2ind(size(badLinks),a,d)+length(px));

        scatter(x,y,'r')
    end
end


badLinkIdx=[badTrackIdx badTrackVecHead];
badLinkIdx=sortrows(badLinkIdx,1); % sorted indices [trackNumber headPosition]

[trackIdxWithBadLink,nBadLinks,whereIdx] = countEntries(badLinkIdx(:,1));
% n links to break creates n+1 segments. but, since we retain the original row
% for the first segment, we only need to add n rows
nRows2add = sum(nBadLinks);

% initialize new matrices to contain new rows
newTrackedFeatureIndx = [trackedFeatureIndx; zeros(nRows2add,nFrames)];
newNnDistFeatures = [nnDistFeatures; nan(nRows2add,nFrames)];
newTrackedFeatureInfo = [trackedFeatureInfo; nan(nRows2add,8*nFrames)];

counter=nTracks+1;
for i=1:length(trackIdxWithBadLink);
    idx = badLinkIdx(badLinkIdx(:,1)==trackIdxWithBadLink(i),2);
    newSegS = [trackSEL(trackIdxWithBadLink(i),1); idx+1];  % new segment start
    newSegE = [idx; trackSEL(trackIdxWithBadLink(i),2)];    % new segment end

    for j=2:length(newSegS) % leave the first one as it is - will contain first segment
        % assign original values to new rows
        newTrackedFeatureIndx(counter,newSegS(j):newSegE(j)) = newTrackedFeatureIndx(trackIdxWithBadLink(i),newSegS(j):newSegE(j));
        newNnDistFeatures(counter,newSegS(j):newSegE(j)) = newNnDistFeatures(trackIdxWithBadLink(i),newSegS(j):newSegE(j));
        newTrackedFeatureInfo(counter,8*(newSegS(j)-1)+1:8*newSegE(j)) = newTrackedFeatureInfo(trackIdxWithBadLink(i),8*(newSegS(j)-1)+1:8*newSegE(j));

        % erase original values from original rows
        newTrackedFeatureIndx(trackIdxWithBadLink(i),newSegS(j):newSegE(j)) = zeros(size(newSegS(j):newSegE(j)));
        newNnDistFeatures(trackIdxWithBadLink(i),newSegS(j):newSegE(j)) = nan(size(newSegS(j):newSegE(j)));
        newTrackedFeatureInfo(trackIdxWithBadLink(i),8*(newSegS(j)-1)+1:8*newSegE(j)) = nan(size(8*(newSegS(j)-1)+1:8*newSegE(j)));

        counter=counter+1;
    end
end

doPlot=0;
if doPlot==1
    
    % plot broken links in red, retained tracks in blue
    px=trackedFeatureInfo(:,1:8:end)';
    py=trackedFeatureInfo(:,2:8:end)';
    
    % limit t
    px=px(5:25,:); py=py(5:25,:);   
    
    x=px(:); x(isnan(x))=[];
	y=py(:); y(isnan(y))=[];
    
    figure
    plot(px(:,:),py(:,:),'r') % original
    hold on
    scatter(x,y,'.b')
    
    px=newTrackedFeatureInfo(:,1:8:end)';
    py=newTrackedFeatureInfo(:,2:8:end)';
    px=px(5:25,:); py=py(5:25,:);
    
    plot(px(:,:),py(:,:),'b') % new

    axis equal
end

%rearrange "newTrackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the
%same time in descending order from longest to shortest.

trackSEL = getTrackSEL(newTrackedFeatureInfo);
[list,indx] = sortrows(trackSEL,[1 -3]);

% rearrange data
newTrackedFeatureIndx = newTrackedFeatureIndx(indx,:);
newNnDistFeatures = newNnDistFeatures(indx,:);
newTrackedFeatureInfo = newTrackedFeatureInfo(indx,:);



