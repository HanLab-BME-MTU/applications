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
% ripleyClust(1,:) = (plotStartFrame:increment:plotEndFrame);
% ripleyClustSlopePoint(1,:) = (plotStartFrame:increment:plotEndFrame);
MPMCount=0;
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
    %in case you want to use jacco's (JvR) ripley:
    % [cpar1, cpar2,cpar3,pvr,dpvr,RememberIndividualRipleys] = JvRClusterQuantRipleyMC (ripMPM, colSize, rowSize, drugTimepoint,1);
    fileseps=findstr(filesep, jobData(1,jobCount).jobpath); 
    savepathparent=jobData(1,jobCount).jobpath(1:fileseps(end)); 
    [RememberIndividualRipleys, jrVar, jrIntegral]=JvRclusterQuantRipleyMCtotalnew(ripMPM,colSize, rowSize, savepathparent);
    % in case you want to use dinahs:
    %[cpar3,pvr,dpvr,cpar2,RememberIndividualRipleys] = ClusterQuantRipley (ripMPM, colSize, rowSize, drugTimepoint,3);

    % if you want to see the individual ripley curves, than use the below
    % function. if you don't please comment!!
    % if you want to save them as well:
    % [RememberIndividualRipleys] = plotripleys (RememberIndividualRipleys);
    fileseps=findstr(filesep, jobData(1,jobCount).jobpath);
    savepathparent=jobData(1,jobCount).jobpath(1:fileseps(end));
    savedir=[savepathparent 'ripsperframejacco'];
    mkdir(savedir);
    try
        NameStruc=findstr('_', jobData(1,jobCount).imagename);
        NameNr=jobData(1,jobCount).imagename((NameStruc(1)+2):(NameStruc(end)-1));
    catch
        NameNr=jobCount
    end
    save ([savedir filesep 'ripsperframe' num2str(NameNr)],'RememberIndividualRipleys');

end  % for jobCount = 1 : length(MPM) 


% now we need to get rid of frames that have NaN and get the right xAxis
% values!

% first put frames and values back together:
ripStartTemp2 = zeros (2,numberOfFrames);
ripClustTemp2 = zeros (2,numberOfFrames);

ripStartTemp2(1,:) = xAxis(1,:);
ripClustTemp2(1,:) = xAxis(1,:);

ripStartTemp2(2,:) = jrVar;
ripClustTemp2(2,:) = jrIntegral;

    
% now separate xAxis and values again (for ptPlotChaosStats)
% and put them in a structure directly
chaosStats.ripleySlopeInclin = ripClustTemp2(2,:);
chaosStats.ripleySlopeStart = ripStartTemp2(2,:);

xAxis.Start = ripStartTemp2(1,:);
xAxis.Slope = ripClustTemp2(1,:);


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
