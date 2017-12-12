% pcTracking
% Input: detected cells
% Output: trajectories

% Assaf Zaritsky, July 2015
function [doneTracking] = pcTracking(params,dirs)

if params.always
    unix(sprintf('rm %strackingOutput.mat',[dirs.tracking filesep]));    
    fprintf('tracking (always): clean output directories\n'); 
end

doneTracking = false;

% params.always = true;

matchingThresholdInPixels = 20/params.pixelSize; % 20 um

%% Prepare for tracking

outFname = [dirs.tracking 'trackingOutput.mat'];
outFnameLong = [dirs.tracking 'trackingOutputLong.mat'];

if exist(outFname,'file') && exist(outFnameLong,'file') && ~params.always
    return;
end

doneTracking = true;

% Tracking parameters
[gapCloseParam,costMatrices, ...
    kalmanFunctions, ...
    probDim,verbose ...
    ] = tracker_paramPhillipe(matchingThresholdInPixels);

saveResults.dir =  dirs.tracking; %directory where to save input and output
saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved

% prepare detection data
ntime = params.nTime-2;
% movieInfoPSD = nan(1,ntime);

if ~isfield(params,'sTime')
    params.sTime = 1;
end

fprintf('preparing for tracking\n');  
for t = params.sTime : ntime    
    PPBacktrackFname = [dirs.detectPPData sprintf('%03d',t) '_backtrack.mat'];
    load(PPBacktrackFname); % detections,combinedImage,pstruct,filterMask
    
    nDetections = length(detections.x) + length(detections.xBacktrack);
    
    % TODO: use cell segmentation to resolve multiple detections!!
    %     assert(multipleDetections([detections.x detections.xBacktrack],[detections.y detections.yBacktrack],matchingThresholdInPixels));
    
    movieInfoPSD(t).xCoord = [[detections.x detections.xBacktrack];params.patchSize .* ones(1,nDetections)]';
    movieInfoPSD(t).yCoord = [[detections.y detections.yBacktrack];params.patchSize .* ones(1,nDetections)]';
    
    movieInfoPSD(t).amp = ones(nDetections,2);
    movieInfoPSD(t).int = ones(nDetections,2);    
end    
    

%% Actual tracking

fprintf('actual tracking\n');     
[tracksFinal,kalmanInfoLink,errFlag] = ...
    trackCloseGapsKalmanSparse(movieInfoPSD, ...
    costMatrices,gapCloseParam,kalmanFunctions,...
    probDim,saveResults,verbose);

allTrajectories = TracksHandle(tracksFinal);

[allTrajectoriesByStartTime,trajStats] = getTrajectoriesDS(allTrajectories,params.minLengthTH);
[allTrajectoriesByStartTimeLong,trajStatsLong] = getTrajectoriesDS(allTrajectories,params.minLengthTHLong);

save(outFname,'allTrajectories','allTrajectoriesByStartTime','trajStats');
save(outFnameLong,'allTrajectories','allTrajectoriesByStartTimeLong','trajStatsLong');

fprintf('done tracking\n');
end

%% 
function [allTrajectoriesByStartTime,trajStats] = getTrajectoriesDS(allTrajectories,minLengthTH)
nTraj = length(allTrajectories);
trajStats.nTrajFinal = 0;
trajStats.accTrajLength = 0;
trajStats.nSubTraj = 0;
% trajStats.subTrajLength = [];
trajStats.trajLength = [];
allTrajectoriesByStartTime = cell(1,allTrajectories.numTimePoints);
% nID = 0;
for i = 1 : nTraj
    curTraj = allTrajectories(i);
    
    trajStats.trajLength = [trajStats.trajLength curTraj.lifetime];
    
    if curTraj.lifetime < minLengthTH
        continue;
    end        
        
    curTraj = fillNans(curTraj);
    
    nSubTraj = floor(curTraj.lifetime/minLengthTH);
    trajStats.nSubTraj = trajStats.nSubTraj + nSubTraj;
    startFrame = curTraj.segmentStartFrame; % subtrajectories
    
    for k = 1 : nSubTraj
        curSubTrajStart = startFrame + (k-1) * floor(curTraj.lifetime/nSubTraj);
        curSubTrajEnd = startFrame + k * floor(curTraj.lifetime/nSubTraj) - 1;
        if isfield(allTrajectoriesByStartTime{curSubTrajStart},'trajectories')
            nTrajStartTime = length(allTrajectoriesByStartTime{curSubTrajStart}.trajectories);
        else
            nTrajStartTime = 0;
        end
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).ind = i;
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).len = curSubTrajEnd - curSubTrajStart + 1;
        subTrajInds = (curSubTrajStart:curSubTrajEnd) - startFrame + 1;
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).x = curTraj.x(subTrajInds);
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).y = curTraj.y(subTrajInds);
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).startFrame = curSubTrajStart;
        allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).endFrame = curSubTrajEnd;
        
        %         nID = nID + 1;
        %         allTrajectoriesByStartTime{curSubTrajStart}.trajectories(nTrajStartTime+1).id = nID;
                
       
        %         startFrame = curTraj.segmentStartFrame;
        %         nTrajStartTime = length(allTrajectoriesByStartTime{startFrame});
        %         allTrajectoriesByStartTime{startFrame}.trajectories(nTrajStartTime+1).ind = i;
        %         allTrajectoriesByStartTime{startFrame}.trajectories(nTrajStartTime+1).len = curTraj.lifetime;
        %         allTrajectoriesByStartTime{startFrame}.trajectories(nTrajStartTime+1).x = curTraj.x;
        %         allTrajectoriesByStartTime{startFrame}.trajectories(nTrajStartTime+1).y = curTraj.y;
        %         allTrajectoriesByStartTime{startFrame}.trajectories(nTrajStartTime+1).endFrame = curTraj.endFrame;
        trajStats.nTrajFinal = trajStats.nTrajFinal + 1;
    end    
    trajStats.accTrajLength = trajStats.accTrajLength + curTraj.lifetime;
end

end

%%
function [curTraj] = fillNans(curTraj)
for i = 1 : length(curTraj.x)
    if isnan(curTraj.x(i))
        curTraj.x(i) = curTraj.x(i-1);
    end
    if isnan(curTraj.y(i))
        curTraj.y(i) = curTraj.y(i-1);
    end
end
end

%%
function [res] = multipleDetections(xs,ys,matchingThresholdInPixels)
DD = pdist2([xs',ys'],[xs',ys']);
res = true;
for i = 1 : length(xs)
    if min([DD(i,1:i-1) DD(i,i+1:end)]) < matchingThresholdInPixels
        res = false;
            return;
    end
    %     for j = i+1 : length(xs)
    %         if (abs(xs(j) - xs(i)) <= matchingThresholdInPixels) && (abs(ys(j) - ys(i)) <= matchingThresholdInPixels)
    %             res = false;
    %             return;
    %         end
    %     end
end
end