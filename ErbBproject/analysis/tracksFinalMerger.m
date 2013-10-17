function [tracksFinal]= tracksFinalMerger(track1,track2,sep,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose)
%tracksFinalMerger mergers tracking results of movies that have been split
%in two pieces. The split should occur at track1(end,:) and track2(1,:)
%The result is a fully merged tracksFinal that can be used in subsequent
%data analysis.
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

for i = 1: numel(track2)
    
    if track2(i).seqOfEvents(1,1) < GapLen +1
        ind2(i)=true;
    end
    
    track2(i).seqOfEvents(1:2,1) = track2(i).seqOfEvents(1:2,1) + sep;
end

lastFrame = vertcat(track2(i).seqOfEvents);
lastFrame = max(lastFrame(:,1));

% Identify tracks in both "movies" that should be considered for merging
% and remove these tracks from track1 and track2
ind1 = vertcat(track1.seqOfEvents);
ind1 = ind1(2:2:end,1) > sep - GapLen;

savedTracks = vertcat(track1(ind1),track2(ind2));

track1 = track1(~ind1);
track2 = track2(~ind2);

% unravel tracks and make a new point list into movieInfo format

movieInfo = struct('xCoord',[],'yCoord',[],'amp',[]);
movieInfo = repmat(movieInfo,[LastFrame,1]);

for i = 1:numel(savedTracks)
    temp = savedTracks(i);
    index = [temp.seqOfEvents(1,1):temp.seqOfEvents(1,2)]; %all the frames in which this track appears
    for j = 1: numel(index)
        movieInfo(index(j)).xCoord =  vertcat(movieInfo(index(j)).xCoord, ...
            [temp.tracksCoordAmpCG(8*(j-1)+1), temp.tracksCoordAmpCG(8*(j-1)+5)]);
        
        movieInfo(index(j)).yCoord =  vertcat(movieInfo(index(j)).yCoord, ...
            [temp.tracksCoordAmpCG(8*(j-1)+2), temp.tracksCoordAmpCG(8*(j-1)+6)]);
        
        movieInfo(index(j)).amp =  vertcat(movieInfo(index(j)).amp, ...
            [temp.tracksCoordAmpCG(8*(j-1)+4), temp.tracksCoordAmpCG(8*(j-1)+8)]);
    end
    
end



% use tracker to track, merge and gapclose this new point list
    [tracksMerge,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

%recombine with tracks 1 and 2

temp = vertcat(track1,track2,tracksMerge);

ind = vertcat(temp.seqOfEvents);
ind = ind(1:2:end-1,1);
[sorted, ix]= sort(ind);
tracksFinal = temp(ix);
end