function [newTrackedFeatureInfo,newTrackedFeatureIndx,newNnDistFeatures]=breakNonlinearTracks(trackedFeatureInfo,trackedFeatureIndx,nnDistFeatures)
%breaks up tracks if not following roughly linear behavior

% this function looks for instances where the vectors created by consecutive 
% features in a track have a negative dot product, indicating that they do
% not follow a unidirectional path. for each broken track two new tracks
% are born. output has the same format as input. note that nnDistFeatures
% may be nonsensical...need to check with Khuloud whether these should be
% recalculated.  many links are now shorter than minLength, but those can
% be filtered out by post-processing step.

% currently called in conditional statement during trackCloseGapsKalman


%get total number of tracks
[numTracks, numFrames] = size(trackedFeatureIndx);

% find track start and end frames
trackSEL = getTrackSEL(trackedFeatureInfo);

% extract coordinates
px=trackedFeatureInfo(:,1:8:end);
py=trackedFeatureInfo(:,2:8:end);

% get vector coordinates between linked features
vx=diff(px,1,2);
vy=diff(py,1,2);

% take dot product of each pair of consecutive vectors along all the tracks
dotProd = vx(:,1:end-1).*vx(:,2:end) + vy(:,1:end-1).*vy(:,2:end);

% dot < 0 when vectors off by 90 degrees or more - indicates a "bad" (ie
% nonlinear) link
[r c]=find(dotProd<0); c=c+1; % add one to get tail position of vector to break
badLinkIdx=[r c];
badLinkIdx=sortrows(badLinkIdx,1); % sorted indices [trackNumber tailPosition]

[trackIdxWithBadLink,nBadLinks,whereIdx] = countEntries(badLinkIdx(:,1));
% n links to break creates n+1 segments. but, since we retain the original row
% for the first segment, we only need to add n rows
nRows2add = sum(nBadLinks);

% initialize new matrices to contain new rows
newTrackedFeatureIndx = [trackedFeatureIndx; zeros(nRows2add,numFrames)];
newNnDistFeatures = [nnDistFeatures; nan(nRows2add,numFrames)];
newTrackedFeatureInfo = [trackedFeatureInfo; nan(nRows2add,8*numFrames)];

counter=numTracks+1;
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
    % plot broken links in red, kept tracks in blue
    px=trackedFeatureInfo(:,1:8:end)';
    py=trackedFeatureInfo(:,2:8:end)';
    plot(px(:,:),py(:,:),'r') % original

    hold on
    px=newTrackedFeatureInfo(:,1:8:end)';
    py=newTrackedFeatureInfo(:,2:8:end)';
    plot(px(:,:),py(:,:),'b') % new
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



