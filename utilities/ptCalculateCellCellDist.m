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
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
binSize = guiData.binsize;
maxDistance = guiData.maxcellcelldist;

% Determine the movie with the most frames
[longestMPM, mpmLength] = ptMaxMPMLength (MPM);
maxFrames = mpmLength / 2;

% Get start and end frames and increment value
startFrame = jobData(1).firstimg;
endFrame = jobData(longestMPM).lastimg;
increment = jobData(1).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Get pixellength and frame interval
frameInterval = round (jobData(1).timeperframe / 60);    % In minutes
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
    
   % Update the x-axis vector and counter
   iCount = iCount + 1;
   xAxis (iCount) = frameCount;

   for jobCount = 1 : length (cellProps)
 
       % Remove the zero rows from cellProps 
       [notZeroEntryRows, notZeroEntryCols] = find (cellProps{jobCount}(:,:,propCount));
       notZeroEntryRows = unique (notZeroEntryRows);
       cells{jobCount} = cellProps{jobCount}(notZeroEntryRows,:,propCount);

       % Remove the zero rows from clusterProps 
       [notZeroEntryRows, notZeroEntryCols] = find (clusterProps{jobCount}(:,:,propCount));
       notZeroEntryRows = unique (notZeroEntryRows);
       clusters{jobCount} = clusterProps{jobCount}(notZeroEntryRows,:,propCount);

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
   end  % jobCount = 1 : length (cellProps)
   
   % Cat all the matrices that we found together
   allDistances = cat (1, distances{:});
   allCells = cat (1, cells{:});
   
   % Calculate the average distance in micrometer
   averageDist(iCount) = (sum (allDistances) / length (allDistances)) * pixelLength;
   
   % Calculate a distance histogram
   %distanceHist (:,iCount) = hist (distances, binSize);
   
   
%    if frameCount == plotStartFrame
%       for kCount = 1 : size(cells,1)
%          distAllCells = sqrt ((allCells(kCount,1) - allCells(:,1)).^2 + ...
%                               (allCells(kCount,2) - allCells(:,2)).^2);
%    
%          % Kick out all the ones that are > maxDistance
%          neighbours = distAllCells (find (distAllCells < (3*str2num(maxDistance))));
%          
%          % Sort the neighbours vector
%          sortedNeighbours = sort(neighbours);
%          
%          % Kick out the zero entry (the one distanced with itself)
%          sortedNeighbours(find (sortedNeighbours == 0)) = [];
%          
%          % Take the 3 nearest neighbours, but leave out the zero entry
%          if length (sortedNeighbours) >= 3
%             nearNeighbours (kCount,:) = sortedNeighbours(1:3);
%          else
%             nearNeighbours(kCount,:) = [0 , 0 , 0];
%             if length (sortedNeighbours) > 0
%                nearNeighbours (kCount,1: length(sortedNeighbours)) = sortedNeighbours;
%             end
%          end
%       end
%    elseif frameCount > plotStartFrame
%       % Dummy
%    end
   
end 

% Prepare output values
cellCellDistStats.averageDist = averageDist;
