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
% Johan de Rooij        jun 05          removing ripMPM mistake and
%                                       changing averaging strategie.
%                                       NB: increment = 1
%                                       user input is not used.
%                                       Has to be implemented.

% Get the latest data from the handles
MPM = handles.allMPM;
clusterProps = handles.allClusterProps;
validFrames = handles.allValidFrames;
jobData = handles.jobData;
guiData = handles.guiData;

% Check that all the images are available for all selected jobs, because we
% need row and colsizes. If not available use a row and colsize calculated
% from the coordinates in the MPM (max x and y values)
for jobCount = 1 : length(MPM)
    if ~isfield(jobData(jobCount),'rowsize') | ~isfield(jobData(jobCount),'colsize')
        if ~jobData(jobCount).imagesavailable
            % Split MPM in x (col) and y (row) values
            xMPM = MPM{jobCount}(:,1:2:end);
            yMPM = MPM{jobCount}(:,2:2:end);

            % Find the max x value
            jobData(jobCount).colsize = max(xMPM(:));
            jobData(jobCount).rowsize = max(yMPM(:));
        end
    end
end

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
% increment = jobData(1).increment;
% NB!! here we cheat, because down below, only increment = 1 will work!!
increment = 1;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

%% !! forget about the following untill after averaging!

% Determine the movie with the most frames (que? least frames?)
% [shortestMPM, mpmLength] = ptMinMPMLength (MPM);
% minFrames = mpmLength / 2;
%[shortestMovie, minFraceilmes] = ptMinMovieLength (validFrames);

% Make sure we only process up to the shortest MPM else the averaging will
% not work correctly
% if plotEndFrame > minFrames
%    plotEndFrame = minFrames;
% end

% Get drugtime point
drugTimepoint = guiData.drugtimepoint;

% Get max MPM length
[mpmNr, maxLength] = ptMaxMPMLength(MPM);

% Initialize the ripley clustering vectors
ripleyClust = zeros (length(MPM)+1, ceil((plotEndFrame-plotStartFrame)/increment)+1);
ripleyClustSlopePoint = zeros (length(MPM)+1, ceil((plotEndFrame-plotStartFrame)/increment)+1);
% put the frame numbers desired in the first row:
ripleyClust(1,:) = (plotStartFrame:increment:plotEndFrame);
ripleyClustSlopePoint(1,:) = (plotStartFrame:increment:plotEndFrame);
for jobCount = 1 : length(MPM) 

    % Get start and end frames and increment value
    startFrame = jobData(jobCount).firstimg;
    endFrame = jobData(jobCount).lastimg;
    % increment = jobData(jobCount).increment;
    % numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

    % Initialize properties counter depending on radiobutton value
    alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
    if ~alwaysCountFrom1
       MPMCount = ceil ((plotStartFrame - startFrame) / increment);
    else
       MPMCount = ceil ((plotStartFrame - 1) / increment);
    end
    
    % now you need to know whether the frame you wish to start from
    % actually was valid. because otherwise you cannot use it in ripley.
    MPMCount = MPMCount + 1;
    checkStartFrame = find(validFrames{jobCount}(1,:) == MPMCount);
    
    % if it was not valid, increase MPMCount and try the next frame.
    if isempty(checkStartFrame);
       while isempty (checkStartFrame);
           MPMCount = MPMCount + 1;
           checkStartFrame = find(validFrames{jobCount}(1,:) == MPMCount);
       end
       disp(['For Job  ',num2str(jobCount)]);
       disp(['NB!! Due to missing or non-segmentable frames,']); 
       disp(['starting ripley-clustering at frame ',num2str(MPMCount)]);
    end
    
    %same story for the end frame, but decrease counter?!
    Countert = plotEndFrame;
    checkEndFrame = find(validFrames{jobCount}(1,:) == Countert);
    if isempty(checkEndFrame);
        while isempty(checkEndFrame);
            Countert = Countert - 1;
            checkEndFrame = find(validFrames{jobCount}(1,:) == Countert);
        end
        disp(['For Job  ',num2str(jobCount)]);
        disp(['NB!! Due to missing or non-segmentable frames,']); 
        disp (['ending ripley-clustering at frame ',num2str(Countert)]);
    end
    
    % NB!!! here we need to add something that allows checking frames when
    % increment is not 1!! to get the right coordinates form the MPM and
    % store the right xAxis values as well!! for now, omly increment = 1
    % will work!! change also line 49 then!!
    
    % now get the Rip input matrix:
    ripMPM = MPM{jobCount}(:,(checkStartFrame*2)-1:checkEndFrame*2);
            
    
    % Initialize cropped MPM (not needed in johans code.)
    % ripMPM = zeros (size(MPM{jobCount},1),numberOfFrames*2);
    % ripMPM = zeros(size(MPM{jobCount},1), 2*length(validFrames{jobCount}(1,plotStartFrame:plotEndFrame)));
    
    % get xAxis values.
    xAxis = validFrames{jobCount}(1,checkStartFrame:checkEndFrame);
   
    
    % Since the ripley function doesn't know how to handle plot start
    % and end frames, we have to crop the MPM first
    % ripMPM = MPM{jobCount};
    
    % Get row and colsizes
    rowSize = jobData(jobCount).rowsize;
    colSize = jobData(jobCount).colsize;
    frameSize = [colSize rowSize];
    
    % Calculate clustering parameter values
    fprintf (1, 'Performing Ripley clustering job %d...\n', jobCount);
    [cpar1, cpar2,cpar3,pvr,dpvr] = ClusterQuantRipleyMC (ripMPM, colSize, rowSize, drugTimepoint,3);
    
    % Store cpar value in the right columns of the ripley vectors:
    for counterrows = 1 : length(xAxis);
        RipRow = find(ripleyClust(1,:) == xAxis(counterrows));
        ripleyClust(jobCount+1,RipRow) = cpar3(1,counterrows);
        ripleyClustSlopePoint(jobCount+1,RipRow) = cpar2(1,counterrows);
    end
    % now we can use this SlopeStartPoint to calculate a derivative per job
    % and also per frame..
    DerTemp = ripleyClustSlopePoint(jobCount+1,:);
    findZeros = find(DerTemp<0.001);
    DerTemp(findZeros) = NaN;
    for framecount = 1:(length(DerTemp)-1);
        DeltaRipStart(jobCount,framecount) = (abs(DerTemp(framecount+1) - DerTemp(framecount))/increment);
    end    
end  % for jobCount = 1 : length(MPM) 


% now, average all frames in the big ripley matrices, ignoring all zero
% entries in slopepoint, because they come from invalid frames.
% have to check this with Dinah!!, no slope = NaN right?? (not zero..)
% Dinah's comment: In the new Ripley function, for completely scattered
% distributions (cpar3 < 0) cpar2 is automatically set to nan

findZeros = find(ripleyClustSlopePoint<0.001);
ripleyClustSlopePoint(findZeros) = NaN;
ripleyClust(findZeros) = NaN;


nrOfJobs = length(MPM);

if nrOfJobs > 1;
    % take out the y values (row 1 is just framenr)
    ripClustSlopeTemp = ripleyClustSlopePoint(2:nrOfJobs+1,:);
    ripClustTemp = ripleyClust(2:nrOfJobs+1,:);

    % average, use nanmean to ignore NaNs. Is that OK in Clust? it is slightly
    % cheating in ClustDeltaRipStartSlopePoint.. Have to correct that!
    avgripleyClustSlopePoint = nanmean(ripClustSlopeTemp);
    avgripleyClust = nanmean(ripClustTemp);
else
    avgripleyClustSlopePoint = ripleyClustSlopePoint(2,:);
    avgripleyClust = ripleyClust(2,:);
end

% now we need to get rid of frames that have NaN and get the right xAxis
% values!

% first put frames and values back together:
ripStartTemp2 = zeros (2,numberOfFrames);
ripClustTemp2 = zeros (2,numberOfFrames);

ripStartTemp2(1,:) = ripleyClust(1,:);
ripClustTemp2(1,:) = ripleyClust(1,:);

ripStartTemp2(2,:) = avgripleyClustSlopePoint;
ripClustTemp2(2,:) = avgripleyClust;

% now take out the columns of these new matrices that have NaN for value.
% we could also choose to leave them in, but when there are not too many
% bad frames in a row, this will be better for the derivatie calculation!

findNaNStart = find(isnan(avgripleyClustSlopePoint));
findNaNClust = find(isnan(avgripleyClust));

for counterNaNs = 1 : length(findNaNStart);
    invalidrow = findNaNStart(counterNaNs);
    ripStartTemp2(:,invalidrow+1-counterNaNs) = [];
end

for counterNaNs = 1 : length(findNaNClust);
    invalidrow = findNaNClust(counterNaNs);
    ripClustTemp2(:,invalidrow+1-counterNaNs) = [];
end
    
% now separate xAxis and values again (for ptPlotChaosStats)
% and put them in a structure directly
chaosStats.ripleySlopeInclin = ripClustTemp2(2,:);
chaosStats.ripleySlopeStart = ripStartTemp2(2,:);

xAxis.Start = ripStartTemp2(1,:);
xAxis.Slope = ripClustTemp2(1,:);

% now for the derivative..a nrofjobs by nrofframes matrix..
if nrOfJobs > 1;
    AvgDRSTemp = nanmean(DeltaRipStart);
else
    AvgDRSTemp = DeltaRipStart;
end

% adding an xAxis to it:
AvgDRSTemp2 = zeros(2,length(AvgDRSTemp));
AvgDRSTemp2(1,:) = ripleyClustSlopePoint(1,1:length(AvgDRSTemp));
AvgDRSTemp2(2,:) = AvgDRSTemp;
% taking out the non-valid entries
findDerNaNs = find(isnan(AvgDRSTemp));
for DerNaNCount = 1 : length(findDerNaNs);
    badrow = findDerNaNs(DerNaNCount);
    AvgDRSTemp2(:,badrow+1-DerNaNCount) = [];
end
% running average, also: not correct, because it runs over non-consec
% frames. only for cleaner output purposes..
AvgDRS = filter(ones(1,5)/5,1,AvgDRSTemp2(2,:));

% Prepare output data
chaosStats.AvgDRS = AvgDRS;
xAxis.Der = AvgDRSTemp2(1,:);

%---------------------------------------------------------------------

function [xMax,yMax] = calcMaxMPMValues(MPM)
% This function calculates the maximum x and y values present in MPM (over
% all frames)

% Split MPM in x (col) and y (row) values
xMPM = MPM(:,1:2:end);
yMPM = MPM(:,2:2:end);

% Find the max x value
xMax = max(xMPM(:));
yMax = max(yMPM(:));
