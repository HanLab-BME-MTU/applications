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
avgSingleDisp = handles.avgVelocityStats.avgSingleDisplacement;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
%increment = jobData(1).increment;
%numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Determine the movie with the most frames
[shortestMPM, mpmLength] = ptMinMPMLength (MPM);
minFrames = mpmLength / 2;

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
    ripMPM = zeros (size(MPM{jobCount},1),numberOfFrames*2);
    
    % Initialize the ripley clustering vector
    ripleyClust = zeros (length(MPM), numberOfFrames);
    simulationAvg = zeros (length(MPM), numberOfFrames);
    simulationMax = zeros (length(MPM), numberOfFrames);
    simulationMin = zeros (length(MPM), numberOfFrames);

    % Make sure the avg displacement vector is also of the right size; this is
    % dependent on the multframevelocity parameter as well
    multframevelocity = guiData.multframevelocity;
    avgDisplacement = zeros (1, numberOfFrames-1);
    avgDisplacement(multframevelocity:end) = avgSingleDisp;
    
    % Fill up the zeros at the start with the first real value
    if multframevelocity > 1
       avgDisplacement(1:multframevelocity-1) = avgDisplacement(multframevelocity);
    end
    
    % Initialize X-axis vector and iCount
    xAxis = zeros (1, numberOfFrames-1);
    iCount = 0;

    % Go through every frame of the set to get x-axis and MPM values
    for frameCount = plotStartFrame : increment : plotEndFrame

       % Increase MPM counter
       MPMCount = MPMCount + 1;
       
       % Store the frame number for display on the x-axis
       iCount = iCount + 1;
       if ~alwaysCountFrom1
           xAxis(iCount) = frameCount;
       else
           xAxis(iCount) = iCount;
       end

       % Since the ripley function doesn't know how to handle plot start
       % and end frames, we have to crop the MPM first
       ripMPM(:,2*iCount-1:2*iCount) = MPM{jobCount}(:,2*MPMCount-1 : 2*MPMCount);
       
    end   % for frameCount
    
    % Get row and colsizes
    rowSize = jobData(jobCount).rowsize;
    colSize = jobData(jobCount).colsize;
    frameSize = [colSize rowSize];
    
    % The following only has to be done if simulation runs are to be made
    if radioButtons.ripleysimplot
    
        % Calculate the average cell diameter from the clusterProps
        avgAreaFrame = zeros(1,size(clusterProps{jobCount},1));
        avgArea = zeros(1,length (clusterProps{jobCount}));
        for kCount = 1 : length (clusterProps{jobCount})
           for clustCount = 1 : size(find(clusterProps{jobCount}(:,2,kCount) > 0),1)
              avgAreaFrame(clustCount) = clusterProps{jobCount}(clustCount,3,kCount) / ...
                                         clusterProps{jobCount}(clustCount,2,kCount);
           end
           avgArea(kCount) = sum(avgAreaFrame) / length(avgAreaFrame);
        end
        avgAreaTotal = sum(avgArea) / length(avgArea);
        avgCellRadius = sqrt(avgAreaTotal/pi);

        % Do the calculations to obtain the simulation values
        confInterval = guiData.ripleyconfint; 
        [simav,simmax,simmin] = calculatedIndifference(ripMPM, avgDisplacement, 2*avgCellRadius, frameSize, confInterval);

        % Store simav, simmax and simmin values
        simulationAvg(jobCount,1:length(simav)) = simav;
        simulationMax(jobCount,1:length(simmax)) = simmax;
        simulationMin(jobCount,1:length(simmin)) = simmin;
    end
    
    % Calculate clustering parameter values
    fprintf (1, 'Performing Ripley clustering...\n');
    [cpar,pvr,dpvr] = ClusterQuantRipley (ripMPM, colSize, rowSize);
    
    % Store cpar value
    ripleyClust(jobCount,1:length(cpar)) = cpar;
    
end  % for jobCount = 1 : length(MPM) 

% The following only has to be done if simulation runs are to be made
if radioButtons.ripleysimplot
   % Calculate the average values over all MPMs 
   avgSimulationAvg = sum(simulationAvg,1) / length(MPM);
   avgSimulationMax = sum(simulationMax,1) / length(MPM);
   avgSimulationMin = sum(simulationMin,1) / length(MPM);
end

% Calculate the average values over all MPMs
avgRipleyClust = sum(ripleyClust,1) / length(MPM);

% Prepare output data
chaosStats.ripleyClustering = avgRipleyClust;

if radioButtons.ripleysimplot
   chaosStats.simulationAvg = avgSimulationAvg;
   chaosStats.simulationMax = avgSimulationMax;
   chaosStats.simulationMin = avgSimulationMin;
end
