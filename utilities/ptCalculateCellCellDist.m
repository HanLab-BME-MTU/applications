function [cellCellDistStats, xAxis] = ptCalculateCellCellDist (handles)
% ptPlotHistValues plots information gathered in cellProps and
% clusterProps regarding cell-cell distances 
%
% SYNOPSIS       [cellCellDistStats, xAxis] = ptCalculateCellCellDist (handles)
%
% INPUT          handles : a structure which contains the information
%                            from the GUI (see below)
%                
% OUTPUT         cellCellDistStats : struct with the following fields:
%                    averageDist : vector with average inter-cell distance
%                xAxis : vector with x-axis values
%
% DEPENDENCIES   ptCalculateCellCellDist.m  uses { nothing }
%                                  
%                ptCalculateCellCellDist.m is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jul 04          First version of ptPlotHistValues
% Andre Kerstens        Aug 04          Added save function for figures
% Andre Kerstens        Sep 04          Complete rewrite of plot functions

% Get the latest data from the handles
MPM = handles.allMPM;
cellProps = handles.allCellProps;
clusterProps = handles.allClusterProps;
frameProps = handles.allFrameProps;
validFrames = handles.allValidFrames;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
binSize = guiData.binsize;
maxDistance = guiData.maxcellcelldist;

% Determine the movie with the most frames
%[longestMPM, mpmLength] = ptMaxMPMLength (MPM);
%maxFrames = mpmLength / 2;

% Determine the movie with the most frames
%[shortestMPM, mpmLength] = ptMinMPMLength (MPM);
%minFrames = mpmLength / 2;
[shortestMovie, minFrames] = ptMinMovieLength (validFrames);

% Make sure we only process up to the shortest MPM else the averaging will
% not work correctly
if plotEndFrame > minFrames
    plotEndFrame = minFrames;
end

% Get start and end frames and increment value
% startFrame = jobData(1).firstimg;
% endFrame = jobData(shortestMPM).lastimg;
% increment = jobData(1).increment;
% numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
startFrame = jobData(shortestMovie).firstimg;
endFrame = jobData(shortestMovie).lastimg;
increment = jobData(shortestMovie).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Get pixellength and frame interval
%frameInterval = round (jobData(1).timeperframe / 60);    % In minutes
pixelLength = jobData(1).mmpixel;

% Initialize properties counter depending on radiobutton value
alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
if ~alwaysCountFrom1
   propCount = ceil ((plotStartFrame - startFrame) / increment);
else
   propCount = ceil ((plotStartFrame - 1) / increment);
end

% Initialize average distance 
averageDist = zeros(1, numberOfFrames);

% Initialize X-axis vector and counter
xAxis = zeros (1, numberOfFrames);
iCount = 0;

% Calculate a number of statistics for every frame
for frameCount = plotStartFrame : increment : plotEndFrame
   
    % Update the properties counter
    propCount = propCount + 1;
    
    % Initialize average sum counter
    averageSum = 0;

    for jobCount = 1 : length (cellProps)
       
       % Find the index where this frame can be found
       frameIndx = find(validFrames{jobCount}(1,:) == propCount);

       if isempty(frameIndx)
           % Frame was bad and cannot be found in cellProps
           cells{jobCount} = [];
           distances{jobCount} = [];
       else
 
           % Remove the zero rows from cellProps 
           [notZeroEntryRows, notZeroEntryCols] = find (cellProps{jobCount}(:,:,frameIndx));
           notZeroEntryRows = unique (notZeroEntryRows);
           cells{jobCount} = cellProps{jobCount}(notZeroEntryRows,:,frameIndx);

           % Remove the zero rows from clusterProps 
           %[notZeroEntryRows, notZeroEntryCols] = find (clusterProps{jobCount}(:,:,frameIndx));
           %notZeroEntryRows = unique (notZeroEntryRows);
           %clusters{jobCount} = clusterProps{jobCount}(notZeroEntryRows,:,frameIndx);

           % Calculate all the distances between cells by doing a Delaunay
           % triangulation
           triangles = delaunay (cells{jobCount}(:,1), cells{jobCount}(:,2));

           % Close the loop by copying the first coordinate set to a new column
           triangles(:,4) = triangles(:,1);

           % Match pairs of points, describing one line (no longer three points, describing a
           % triangle) and throw out the lines (between two points) that occur twice 
           uniqTriang = cat (1, unique(triangles(:,1:2), 'rows'), ...
                                unique(triangles(:,2:3), 'rows'), ...
                                unique(triangles(:,3:4), 'rows'));

           % Prepare storage for the distances
           distances{jobCount} = zeros (length(uniqTriang), 1);

           % Calculate the distances
           for jCount = 1 : length(uniqTriang)   
              distances{jobCount}(jCount) = sqrt ((cells{jobCount}(uniqTriang(jCount,1), 1) - ...
                                                   cells{jobCount}(uniqTriang(jCount,2), 1))^2 + ...
                                                  (cells{jobCount}(uniqTriang(jCount,1), 2) - ...
                                                   cells{jobCount}(uniqTriang(jCount,2), 2))^2);  
           end
           
           % Increase counter used later to calculate average
           averageSum = averageSum + 1;
       end  % if isempty(frameIndx)
    end  % for jobCount = 1 : length (cellProps)

    if averageSum > 0   % Only go on if we have at least 1 good frame in the joblist
   
        % Update the x-axis vector and counter
        iCount = iCount + 1;
        if ~alwaysCountFrom1
           xAxis(iCount) = frameCount;
        else
           xAxis(iCount) = iCount;
        end
        
        % Cat all the matrices that we found together
        allDistances = cat (1, distances{:});
        allCells = cat (1, cells{:});

        % Calculate the average distance in micrometer
        averageDist(iCount) = (sum (allDistances) / length (allDistances)) * pixelLength;
        %averageDist(iCount) = (sum (allDistances) / averageSum) * pixelLength;

    end  % if averageSum > 0
end  % for frameCount = plotStartFrame : increment : plotEndFrame

% Prepare output values
cellCellDistStats.averageDist = averageDist(1:iCount);

% Make sure the x-axis has the correct length
xAxis = xAxis(1:iCount);