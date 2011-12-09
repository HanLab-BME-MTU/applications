function trajStats = classifyEdgeMotion(protSamples)


%assign edge position and protrusion vector standard deviations
edgePosStd = 1; %in pixels
protVecStd = sqrt(2*edgePosStd);

%get number of windows and number of frames
avgNormal = protSamples.avgNormal;
[numWindows,numFrames] = size(avgNormal);

%go over all windows
for iWin = 1 : numWindows
    
    %initialize motion classification vector
    classCurr = NaN(1,numFrames-1);
    
    %get current protrusion normal
    protNorm = avgNormal(iWin,1:end-1)';
    
    %integrate protrusion normals to get position over time
    posNorm = [0; cumsum(protNorm)];

    %determine sign of each protrusion normal
    protNormSign = sign(protNorm);
    
    %extract points of switching
    protNormSignDiff = diff(protNormSign);
    switchPoints = find(protNormSignDiff~=0) + 1;
    switchPoints = [1; switchPoints; numFrames];
    
    %calculate length of intervals between switching
    intervalLength = diff(switchPoints);
    
    %assign as significant any intervals that are >= 3 frame steps
    indxGood = find( intervalLength >= 3 );
    
    
    

    
end


%% ~~~ old stuff ~~~

%reserve memory
% traj = repmat(struct('distance',[],'time',[],'timePoints',[]),numWindows,1);
% 
% %convert protrusions into format for analysis
% for iWin = 1 : numWindows
%     
%     %extract protrusion vector normals
%     protVecNorm = protSamples.avgNormal(iWin,:)';
%     
%     %convert into position and store
%     posNorm = [0; cumsum(protVecNorm(1:end-1))];
%     posNorm = posNorm - min(posNorm) + 1;
%     traj(iWin).distance = [posNorm ones(numFrames,1)];
%     
%     %store time information
%     traj(iWin).time = [(1:numFrames)' zeros(numFrames,1)];
%     traj(iWin).timePoints = (1:numFrames)';
% 
% end
% 
% %run analysis
% trajStats = trajectoryAnalysis(traj);

 