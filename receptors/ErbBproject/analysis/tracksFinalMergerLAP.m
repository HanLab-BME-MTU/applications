function [tracksFinal]= tracksFinalMergerLAP(track1,track2,sep,costMatrices,gapCloseParam)
%tracksFinalMerger mergers tracking results of movies that have been split
%in two pieces. The split should occur at track1(end,:) and track2(1,:)
%The result is a fully merged tracksFinal that can be used in subsequent
%data analysis.
%
% Uses LAP to merge tracks withing GapLen of the end of track1 and the
% begining of track2 and then merges them into one track and returns the
% merged track
%
% Only handles 2D data
% 
% sep is the frame number where the movies are split track1 ends at frame
% sep and track2 starts at sep+1
%
%Note track1 should come before track2 in time

Len1 = size(track1,1);
Len2 = size(track2,1);
GapLen = gapCloseParam.timeWindow;

%Adjust frame time information in track2

ind2 = false(size(track2));
Frame2 = [];

for i = 1: numel(track2)
    
    if track2(i).seqOfEvents(1,1) < GapLen +1
        %finds tracks that may need to be merged
        ind2(i)=true;
        %separates the location at the first point of the track
        Frame2 = vertcat(Frame2, track2(i).tracksCoordAmpCG(1:2));
    end
    
    track2(i).seqOfEvents(1:2,1) = track2(i).seqOfEvents(1:2,1) + sep;
end
 
% lastFrame = vertcat(track2.seqOfEvents);
% lastFrame = max(lastFrame(:,1));

% Identify tracks in both "movies" that should be considered for merging
% and remove these tracks from track1 and track2
ind1 = vertcat(track1.seqOfEvents);
ind1 = ind1(2:2:end,1) > sep - GapLen;

% create a list of points (last known loc of track1 and first known of
% track2) Note this is already handled for track2
Save1 = track1(ind1);
Save2 = track2(ind2);
track1 = track1(~ind1);
track2 = track2(~ind2);


%end-7 is x pos end-6 is y position
Frame1 = arrayfun(@(x) x.tracksCoordAmpCG(end-7:end-6),Save1,'UniformOutput',false);
Frame1 = vertcat(Frame1{:});

%numFeatures
numFeaturesFrame1 = size(Frame1,1);
numFeaturesFrame2 = size(Frame2,1);

%Make a costmatrix for using LAP (adapted from costMatStationaryLink)
searchRadius = costMatrices(1).parameters.searchRadius;

costMat = createDistanceMatrix(Frame1,Frame2);
CMind = costMat > searchRadius; % identifies non allowed event for LAP

%cost is the distance squared
costMat = costMat.^2;

costMat(CMind) = -1; % sets to non allowed for LAP

% Birth and death

maxCost = 2*max(max(costMat(:)),eps);

deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = -1;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = -1;

%get the cost for the lower right block
% costLR = min(min(min(costMat))-1,-1);
costLR = maxCost;
lrBlock = costMat';
lrBlock(lrBlock>0) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%Find matching tracks
[x, y] = lap(costMat, -1, true);

indLap = find(x(1:numFeaturesFrame1)<numFeaturesFrame2);

%to be removed
tbr = [];

for i = indLap'
    temp1 = Save1(i);
    temp2 = Save2(x(i));
    
    seqOfEvents = vertcat(temp1.seqOfEvents,temp2.seqOfEvents);
    if seqOfEvents(3,1)-seqOfEvents(2,1) <= GapLen  %if gap is to large tracks are not merged
        
        temp1.tracksFeatIndxCG = [temp1.tracksFeatIndxCG,temp2.tracksFeatIndxCG];
        
        if seqOfEvents(2,1) <  seqOfEvents(3,1)-1
            gapsize = seqOfEvents(3,1)-1-seqOfEvents(2,1);
            gapfill = NaN(1,gapsize*8);
            temp1.tracksCoordAmpCG = [temp1.tracksCoordAmpCG,gapfill,temp2.tracksCoordAmpCG];
        else
            temp1.tracksCoordAmpCG = [temp1.tracksCoordAmpCG,temp2.tracksCoordAmpCG];
        end
        
        temp1.seqOfEvents = seqOfEvents([1,4],:);
        
        Save1(i) = temp1;
        tbr = [tbr,x(i)];
        
    else
        %do nothing        
    end
    
end

Save2(tbr)=[];


%recombine with tracks 1 and 2

temp = vertcat(track1,track2,Save1,Save2);

ind = vertcat(temp.seqOfEvents);
ind = ind(1:2:end-1,1);
[sorted, ix]= sort(ind);
tracksFinal = temp(ix);
end