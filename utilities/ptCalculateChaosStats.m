function [chaosStats, xAxis] = ptCalculateChaosStats (handles, radioButtons)
% ptCalculateChaosStats plots chaos theory statistics from MPMs. 
%
% SYNOPSIS       ptCalculateChaosStats (handles)
%
% INPUT          handles : a structure which contains the information from the GUI
%                radioButtons : some radiobutton values from the gui
%                
% OUTPUT         chaosStats : struct with following fields:
%                     : vector with clustering parameter values
%
% DEPENDENCIES   ptCalculateChaosStats  uses { ClusterQuantRipley.m }
%                                  
%                ptCalculateChaosStats is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version

% Get the latest data from the handles
MPM = handles.allMPM;
clusterProps = handles.allClusterProps;
validFrames = handles.allValidFrames;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
%increment = jobData(1).increment;
%numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Determine the movie with the most frames
% [shortestMPM, mpmLength] = ptMinMPMLength (MPM);
% minFrames = mpmLength / 2;
[shortestMovie, minFrames] = ptMinMovieLength (validFrames);

% Make sure we only process up to the shortest MPM else the averaging will
% not work correctly
if plotEndFrame > minFrames
    plotEndFrame = minFrames;
end

for jobCount = 1 : length(MPM) 

    % Get start and end frames and increment value
    startFrame = jobData(jobCount).firstimg;
    endFrame = jobData(jobCount).lastimg;
    increment = jobData(jobCount).increment;
    numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

    % Initialize properties counter depending on radiobutton value
    alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
    if ~alwaysCountFrom1
       MPMCount = ceil ((plotStartFrame - startFrame) / increment);
    else
       MPMCount = ceil ((plotStartFrame - 1) / increment);
    end

    % Initialize cropped MPM
    %ripMPM = zeros (size(MPM{jobCount},1),numberOfFrames*2);
    ripMPM = zeros(size(MPM{jobCount},1), 2*length(validFrames{jobCount}(1,:)));
    
    % Initialize the ripley clustering vectors
    ripleyClust = zeros (length(MPM), numberOfFrames);
    ripleyClustSlopePoint = zeros (length(MPM), numberOfFrames);
    
    % Initialize X-axis vector and iCount
    xAxis = zeros (1, numberOfFrames);
    %xAxis = zeros (1, minFrames);
    iCount = 0;

    % Go through every frame of the set to get x-axis and MPM values
    for frameCount = plotStartFrame : increment : plotEndFrame
      
       % Increase MPM counter
       MPMCount = MPMCount + 1;
      
       frameIndx = find(validFrames{jobCount}(1,:) == MPMCount);
       
       if isempty(frameIndx)
          continue;
       end

       % Store the frame number for display on the x-axis
       iCount = iCount + 1;
       if ~alwaysCountFrom1
           xAxis(iCount) = frameCount;
       else
           xAxis(iCount) = iCount;
       end

       % Since the ripley function doesn't know how to handle plot start
       % and end frames, we have to crop the MPM first
       ripMPM(:,2*iCount-1:2*iCount) = MPM{jobCount}(:,2*frameIndx-1 : 2*frameIndx);

    end   % for frameCount
    
    % Get row and colsizes
    rowSize = jobData(jobCount).rowsize;
    colSize = jobData(jobCount).colsize;
    frameSize = [colSize rowSize];
    
    % Calculate clustering parameter values
    fprintf (1, 'Performing Ripley clustering ...\n');
    [cpar,pvr,dpvr,cpar2] = ClusterQuantRipley (ripMPM, colSize, rowSize);
    
    % Store cpar value
    ripleyClust(jobCount,1:length(cpar)) = cpar;
    
    ripleyClustSlopePoint(jobCount,1:length(cpar2)) = cpar2;
    
end  % for jobCount = 1 : length(MPM) 

% Determine the last entry in xAxis
xIndex = find(xAxis);
lastEntry = xIndex(end);

% Calculate the average values over all MPMs
avgRipleyClust = sum(ripleyClust,1) / length(MPM);
avgripleyClustSlopePoint = sum(ripleyClustSlopePoint,1) / length(MPM);

% Prepare output data
chaosStats.ripleySlopeInclin = avgRipleyClust(1:lastEntry);
chaosStats.ripleySlopeStart = avgripleyClustSlopePoint(1:lastEntry);

% Also take valid xAxis entries
xAxis = xAxis(1:lastEntry);
