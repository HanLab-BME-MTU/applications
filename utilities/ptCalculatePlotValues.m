function [cellClusterStats, areaStats, perimeterStats, xAxis] = ptCalculatePlotValues (handles)
% ptCalculatePlotValues plots information gathered in cellProps and
% clusterProps during the cell tracking (PolyTrack) 
%
% SYNOPSIS       [cellClusterStats, areaStats, perimeterStats, xAxis] = ptCalculatePlotValues (handles)
%
% INPUT          handles : a structure which contains the information from the GUI 
%                
% OUTPUT         cellClusterStats : vector with cell and cluster stats per frame
%                areaStats : vector with area stats per frame
%                perimeterStats : vector with perimeter stats per frame
%                xAxis : vector with x-axis values
%
% DEPENDENCIES   ptCalculatePlotValues.m  uses { nothing }
%                                  
%                ptCalculatePlotValues.m is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file
% Andre Kerstens        Jul 04          Added pixelarea all cells percentage
% Andre Kerstens        Jul 04          Added plots for area/convex-hull-area
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

% Get row and colsizes for the convex hull calculation
rowSize = jobData(1).rowsize;
colSize = jobData(1).colsize;

% Calculate the total area of the frame (in um^2)
totalAreaFrame = (rowSize * colSize) * (pixelLength^2);

% Initialize matrices
cellAmount = zeros(1, numberOfFrames);
clusterAmount = zeros(1, numberOfFrames);
sumCellsInCluster = zeros(1, numberOfFrames);
cellsPerCluster = zeros(1, numberOfFrames);
singleCellAmount = zeros(1, numberOfFrames);
percentageSingleCells = zeros(1, numberOfFrames);
percentageClusteredCells = zeros(1, numberOfFrames);
sumClusterArea = zeros(1, numberOfFrames);
areaPerCluster = zeros(1, numberOfFrames);
sumSingleCellArea = zeros(1, numberOfFrames);
areaPerSingleCell = zeros(1, numberOfFrames);
totalAreaAllCells = zeros(1, numberOfFrames);
percentageAreaAllCells = zeros(1, numberOfFrames);
sumClusterPerimeter = zeros(1, numberOfFrames);
perimeterLength = zeros(1, numberOfFrames);
sumClusterArea = zeros(1, numberOfFrames);
perimeterDivArea = zeros(1, numberOfFrames);
convexHullData = zeros(3,numberOfFrames);

% Initialize properties counter
propCount = ceil ((plotStartFrame - startFrame) / increment);

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
   end
   
   % Cat all the matrices that we found together
   allCells = cat(1,cells{:});
   allClusters = cat(1,clusters{:});
   
   % Calculate the amount of all cells per frame
   cellAmount(iCount) = round(sum(allClusters(:, 2)) / length(cellProps));

   % Calculate the amount of clusters per frame (a cluster should contain at
   % least 2 nuclei (and therefore we use > 1)
   clusterAmount(iCount) = round(size(allClusters(find (allClusters(:,2) > 1)), 1) / length(cellProps));
   
   % Calculate the average amount of cells per cluster
   sumCellsInCluster(iCount) = sum (allClusters (find (allClusters(:,2) > 1), 2));
   if clusterAmount(iCount)
      cellsPerCluster(iCount) = round (sumCellsInCluster(iCount) / clusterAmount(iCount));
   else
      cellsPerCluster(iCount) = 0;
   end
   
   % Calculate the amount of single cells per frame. This is in principle a
   % cluster with only one nuclei
   singleCellAmount(iCount) = round(size(allClusters(find (allClusters (:,2) == 1)), 1) / length(cellProps));
   
   % Calculate the percentage of single cells
   percentageSingleCells(iCount) = (singleCellAmount(iCount) / ...
                                    (singleCellAmount(iCount) + sumCellsInCluster(iCount))) * 100.0;
                                    
   % Calculate the percentage of clustered cells
   percentageClusteredCells(iCount) = 100.0 - percentageSingleCells(iCount);
   
   % Calculate the average area per cluster
   sumClusterArea(iCount) = sum (allClusters (find (allClusters (:,2) > 1), 3));
   if clusterAmount(iCount) ~= 0
      areaPerCluster(iCount) = (sumClusterArea(iCount) / clusterAmount(iCount)) * (pixelLength^2);
   else
      areaPerCluster(iCount) = 0;
   end
   
   % Calculate the average area per single cell
   sumSingleCellArea(iCount) = sum (allClusters (find (allClusters (:,2) == 1), 3));
   if singleCellAmount(iCount) ~= 0
      areaPerSingleCell(iCount) = (sumSingleCellArea(iCount) / singleCellAmount(iCount)) * (pixelLength^2);
   else
      areaPerSingleCell(iCount) = 0;
   end
   
   % Calculate area taken up by single and clustered cells together (in um^2)
   totalAreaAllCells(iCount) = (sumSingleCellArea(iCount) + sumClusterArea(iCount)) * (pixelLength^2);
   percentageAreaAllCells(iCount) = (totalAreaAllCells(iCount) / totalAreaFrame) * 100.0;
   
   % Calculate average perimeter length
   sumClusterPerimeter(iCount) = sum (allClusters (find (allClusters (:,2) > 1), 4));
   if clusterAmount(iCount) ~= 0 
      perimeterLength(iCount) = sumClusterPerimeter(iCount) / clusterAmount(iCount);
   else
      perimeterLength(iCount) = 0;
   end
   
   % Calculate perimeter / area for clusters
   if sumClusterArea(iCount) ~= 0
      perimeterDivArea(iCount) = sumClusterPerimeter(iCount) / sumClusterArea(iCount);
   else
      perimeterDivArea(iCount) = 0;
   end
   
   % Initialize tmpConvexHullData
   clear tmpConvexHullData catConvexHullData;
   
   for jobCount = 1 : length(frameProps)
      % Calculate the area and convex-hull-area in um^2
      tmpConvexHullData{jobCount}(1) = frameProps{jobCount}(1, 1, iCount) * (pixelLength^2);
      tmpConvexHullData{jobCount}(2) = frameProps{jobCount}(1, 2, iCount) * (pixelLength^2);
   
      % Calculate the ratio area / convex-hull-area in percent
      tmpConvexHullData{jobCount}(3) = (tmpConvexHullData{jobCount}(1) / tmpConvexHullData{jobCount}(2)) * 100;
   end
   
   % Cat the results
   catConvexHullData = cat(1, tmpConvexHullData{:});
   
   % Average data
   convexHullData(1,iCount) = sum(catConvexHullData(:,1)) / length(frameProps);
   convexHullData(2,iCount) = sum(catConvexHullData(:,2)) / length(frameProps);
   convexHullData(3,iCount) = sum(catConvexHullData(:,3)) / length(frameProps);
end 

% Put the convex hull data in a matrix
%allConvexHullData = cat(1,convexHullData{:});

% Store all calculated values in the cellClusterStats struct
cellClusterStats.cellAmount = cellAmount;
cellClusterStats.clusterAmount = clusterAmount;
cellClusterStats.cellsPerCluster = cellsPerCluster;
cellClusterStats.singleCellAmount = singleCellAmount;
cellClusterStats.percentageSingleCells = percentageSingleCells;
cellClusterStats.percentageClusteredCells = percentageClusteredCells;

% Store all calculated values in the areaStats struct
areaStats.areaPerSingleCell = areaPerSingleCell;
areaStats.areaPerCluster = areaPerCluster;
areaStats.totalAreaAllCells = totalAreaAllCells;
areaStats.percentageAreaAllCells = percentageAreaAllCells;
areaStats.convexHullData = convexHullData;

% Store all calculated values in the perimeterStats struct
perimeterStats.perimeterLength = perimeterLength;
perimeterStats.perimeterDivArea = perimeterDivArea;