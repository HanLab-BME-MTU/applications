function [chaosStats, xAxis] = ptCalculateChaosStats (handles)
% ptCalculateChaosStats plots chaos theory statistics from MPMs. 
%
% SYNOPSIS       ptCalculateChaosStats (handles)
%
% INPUT          handles : a structure which contains the information from the GUI
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
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
increment = jobData(1).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Determine the movie with the most frames
[longestMPM, mpmLength] = ptMaxMPMLength (MPM);
maxFrames = mpmLength / 2;

% Initialize the ripley clustering vector
ripleyClust = zeros (length(MPM), numberOfFrames);

for jobCount = 1 : length(MPM) 

    % Get start and end frames and increment value
    startFrame = jobData(jobCount).firstimg;
    endFrame = jobData(jobCount).lastimg;
    increment = jobData(jobCount).increment;
    numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

    % Initialize properties counter depending on radiobutton value
    alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
    if ~alwaysCountFrom1
       MPMCount = ceil ((plotStartFrame - 1) / increment);
    else
       MPMCount = ceil ((plotStartFrame - 1) / increment);
    end

    % Initialize cropped MPM
    ripMPM = zeros (size(MPM{jobCount},1),numberOfFrames*2);
    
    % Initialize X-axis vector and iCount
    xAxis = zeros (1, numberOfFrames-1);
    iCount = 0;

    % Go through every frame of the set to get x-axis and MPM values
    for frameCount = plotStartFrame : increment : plotEndFrame

       % Increase MPM counter
       MPMCount = MPMCount + 1;
       
       % Store the frame number for display on the x-axis
       iCount = iCount + 1;
       xAxis (iCount) = frameCount;

       % Since the ripley function doesn't know how to handle plot start
       % and end frames, we have to crop the MPM first
       ripMPM(:,2*iCount-1:2*iCount) = MPM{jobCount}(:,2*MPMCount-1 : 2*MPMCount);
       
    end   % for frameCount
    
    % Get row and colsizes
    rowSize = jobData(jobCount).rowsize;
    colSize = jobData(jobCount).colsize;
    
    % Calculate clustering parameter values
    [cpar,pvr,dpvr] = ClusterQuantRipley (ripMPM, colSize, rowSize);
    
    % Store cpar value
    ripleyClust(jobCount,1:length(cpar)) = cpar;
    
end  % for jobCount = 1 : length(MPM) 

% Calculate the average over all MPMs
avgRipleyClust = sum(ripleyClust,1) / length(MPM);

% Prepare output data
chaosStats.ripleyClustering = avgRipleyClust;
