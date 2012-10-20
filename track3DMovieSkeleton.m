function movieData = track3DMovieSkeleton(movieData,paramsIn)

%Somebody else will clean this up and add documentation later, right? ;)



%% ------------- Parameters -------------- %%

fileName = 'pruned skeleton tracking';%String for naming output file

iProcChan = 1;%Channel with wich non-specific processing is associated
    
%% ------------- Input ------------ %%

if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object for a 3D movie!')
end

if nargin < 2
    paramsIn = [];
end

iProc = movieData.getProcessIndex('SkeletonTrackingProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SkeletonTrackingProcess(movieData));                                                                                                 
end

%Parse input, store in parameter structure. We also need to get the
%additional options for passing to track3DSkeletonGraphs.m
[p,trackParam] = parseProcessParams(movieData.processes_{iProc},paramsIn,true);

p.ChannelIndex = iProcChan;

%% ------------ Init --------------- %%

%Make sure the movie has had skeleton pruning performed
iSkelProc = movieData.getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);

if isempty(iSkelProc) || ~movieData.processes_{iSkelProc}.checkChannelOutput(iProcChan)
    error('The input movie has no valid skeleton pruning! Please run skeleton pruning prior to trackign!')
end


%Load the skeleton graphs from each frame
nFrames = movieData.nFrames_;
allSkel(1:nFrames,1) = struct('vertices',[],'edges',[],'edgePaths',[],'edgeLabels',[]);
for iFrame = 1:nFrames
    
    allSkel(iFrame) = movieData.processes_{iSkelProc}.loadChannelOutput(p.ChannelIndex,iFrame);
        
end

%Use movie units to convert tracking params to pixel units
%maxDisp = p.MaxDisp / movieData.pixelSize_ * movieData.timeInterval_;%Z-size has already been corrected, so multiply by XY size only
%TEMP
maxDisp = 30;


%% ---------- Tracking ---------- %%

%Just call the tracking function.

skelTracks = track3DSkeletonGraphs(allSkel,'MaxDisp',maxDisp);



%% -------------- Output ------------ %%
%Write that shit to disk.
mkClrDir(p.OutputDirectory);
outPath = [p.OutputDirectory filesep fileName '.mat'];
save(outPath,'skelTracks');
movieData.processes_{iProc}.setPara(p);%Store any tracking parameters
movieData.processes_{iProc}.setOutFilePath(p.ChannelIndex,outPath);
movieData.processes_{iProc}.setDateTime;
movieData.save;

