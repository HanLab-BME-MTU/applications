function windowMotionType = classifyEdgeMotion(protSamples,doPlot)
%CLASSIFYEDGEMOTION classifies edge motion in each window based on protrusion vector samples
%
%SYNPOSIS edgeMotionType = classifyEdgeMotion(protSamples,doPlot)
%
%INPUT  protSample      : The protrusion samples as output by the windowing
%                         software.
%       doPlot          : Flag with value 1 to plot and 0 not to plot
%                         classification. Optional. Default: 0.
%OUTPUT windowMotionType: 2D array same size and protSamples.avgNormal
%                         indicating motion classification, where +1 means
%                         protrusion, -1 means retraction, 0 means pause,
%                         and NaN indicates windows without a prutrusion
%                         vector.
%
%Khuloud Jaqaman, December 2011

%% Input

if nargin < 1
    error('--classifyEdgeMotion: Please enter protrusion samples!');
end

if nargin < 2
    doPlot = 0;
end

%% Classification

%assign edge position and protrusion vector standard deviations
edgePosStd = 1; %in pixels
protVecStd = sqrt(2*edgePosStd);

%get number of windows and number of frames
avgNormal = protSamples.avgNormal;
[numWindows,numFrames] = size(avgNormal);

%initialize output
windowMotionType = NaN(size(avgNormal));

%go over all windows
for iWin = 1 : numWindows
    
    %initialize motion classification vector
    classCurr = NaN(1,numFrames-1);
    
    %get current protrusion normal
    protNorm = avgNormal(iWin,1:end-1)';
    
    %determine sign of each protrusion normal
    protNormSign = sign(protNorm);
    
    %extract points of switching
    protNormSignDiff = diff(protNormSign);
    switchPoints = find(protNormSignDiff~=0) + 1;
    switchPoints = [1; switchPoints; numFrames]; %#ok<AGROW>
    
    %calculate length of intervals between switching
    intervalLength = diff(switchPoints);
    
    %get number of intervals
    numIntervals = length(intervalLength);
    
    %calculate position change (i.e. cumulative protrusion/retraction) in
    %each interval
    intervalProt = NaN(numIntervals,1);
    for iInt = 1 : numIntervals
        intervalProt(iInt) = sum(protNorm(switchPoints(iInt):switchPoints(iInt+1)-1));
    end
    
    %initialize interval classification vector
    intervalClass = NaN(size(intervalLength));
    
    %assign as significant any intervals that are >= 3 frame steps
    indx = find(intervalLength >= 3);
    intervalClass(indx) = sign(intervalProt(indx)); % +1 = protrusion, -1 = retraction
    
    %go over intervals that are <= 2 frame steps and test for significance
    indx2 = find(intervalLength<=2);
    oneMinusPvalue = normcdf(abs(intervalProt(indx2)),0,protVecStd);
    indx = indx2(oneMinusPvalue>0.95);
    intervalClass(indx) = sign(intervalProt(indx));% +1 = protrusion, -1 = retraction
    indx = indx2(oneMinusPvalue<=0.95);
    intervalClass(indx) = 0; % 0 = pause
    
    %from intervalClass, assign classification to individual frame steps
    for iInt = 1 : numIntervals
        classCurr(switchPoints(iInt):switchPoints(iInt+1)-1) = intervalClass(iInt);
    end
    
    %save in output variable
    windowMotionType(iWin,1:end-1) = classCurr;
    
    if doPlot
        
        %integrate protrusion normals to get position over time for
        %plotting
        posNorm = nancumsumwithzero(protNorm);
        
        %plot protrusion in green, retraction in red and pause in cyan
        %intervals with missing points stay empty
        figure, hold on
        plot(posNorm,'kd'), myErrorbar(posNorm,edgePosStd*ones(numFrames,1))
        indx = find(classCurr==1);
        for iInt = indx
            plot([iInt iInt+1],posNorm(iInt:iInt+1),'g')
        end
        indx = find(classCurr==-1);
        for iInt = indx
            plot([iInt iInt+1],posNorm(iInt:iInt+1),'r')
        end
        indx = find(classCurr==0);
        for iInt = indx
            plot([iInt iInt+1],posNorm(iInt:iInt+1),'c')
        end
        
    end %(if doPlot)
    
end %(for iWin = 1 : numWindows)

%% ~~~ the end ~~~

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

 