function [ success] = unitTest_tracksFinalMerger(nFrames)
% unitTest_tracksFinalMerger test the functionality of the function
% tracksFinalMerger
%
% Written by
% Jeffrey Werbin 2013/10/17
%

success =false;

%set tracking parameters
%% Creates parameter structures for tracking
%% general gap closing parameters
gapCloseParam.timeWindow = 3; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatStationaryLink';

%parameters
parameters.searchRadius = 2;

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatStationaryCloseGaps';

%parameters
parameters.searchRadius = 2; 
parameters.gapPenalty = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions = [];

%% additional input
saveResults = 0;

%verbose
verbose = 1;

%problem dimension
probDim = 2;



% Generate points centers at all unit intersections on a 100x100 latice
[x,y] = ind2sub([40,40],find(ones(50)));
pnts = [x*2.5,y*2.5];

%determine when each track appears and disapears
limits = sort(random('unid',nFrames,[numel(x),2]),2);

fake = sparse(nFrames,numel(x));

for i=1:numel(x)
    fake(limits(i,1):limits(i,2),i)=1;
end

%add drift to appease Khuloud's tracker
drift = 0:0.01:nFrames*0.01;

%create movieTemplate
movieTemplate = struct('xCoord',[],'yCoord',[],'amp',[]);
movieInfo = repmat(movieTemplate,[nFrames,1]);


for i=1:nFrames
    ind = find(fake(i,:));
    movieInfo(i).amp = [50*ones(size(ind))', 5*ones(size(ind))'];
    movieInfo(i).xCoord = [(pnts(ind,1)+drift(i)), 0.2*ones(size(ind))'];
    movieInfo(i).yCoord = [(pnts(ind,2)+drift(i)), 0.2*ones(size(ind))'];
end

%Track whole thing first
tic
[tWhole,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
'Whole'
numel(tWhole)
toc


%split in the middle
sep = floor(nFrames/3);

tic
[t1,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1:sep),costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
[t2,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(sep+1:end),costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
'two halves'
numel(t1)+numel(t2)
toc

tic
[tMerge] = tracksFinalMerger(t1,t2,sep,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
'Merge'
numel(tMerge)
toc

%plots the two
figure,hold;
for i=1:numel(tWhole)
    plot(tWhole(i).tracksCoordAmpCG(1:8:end),tWhole(i).tracksCoordAmpCG(2:8:end),'r');
end

for i = 1:numel(tMerge)
    plot(tMerge(i).tracksCoordAmpCG(1:8:end),tMerge(i).tracksCoordAmpCG(2:8:end)+0.5,'b');
end


if numel(tWhole) == numel(tMerge)
    success = true;
end


end

