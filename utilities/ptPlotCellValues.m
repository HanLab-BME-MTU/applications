function ptCalculatePlotValues (ptPostpro)
% ptCalculatePlotValues plots information gathered in cellProps and
% clusterProps during the cell tracking (PolyTrack) 
%
% SYNOPSIS       ptCalculatePlotValues (ptPostpro)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
%                
% OUTPUT         None (plots are directly shown on the screen) 
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

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
rowSize = ptPostpro.rowsize;
colSize = ptPostpro.colsize;

% Calculate the total area of the frame
totalAreaFrame = rowSize * colSize;

% Get radio button values
cellClusterPlot = ptPostpro.cellclusterplot;
areaPlot = ptPostpro.areaplot;
perimeterPlot = ptPostpro.perimeterplot;

% Get the cell and cluster properties for all frames. These are in the format:
%
%    cellProps :
%      cellProps (:,1,:) = coord (:,1);
%	   cellProps (:,2,:) = coord (:,2);
%	   cellProps (:,3,:) = clusterNr (:);  (number of cluster - label)
%
%    clusterProps :
%      clusterProps (:,1,:) = uniqClusterNr (:);          (which cells are in this cluster)
%      clusterProps (:,2,:) = numberOfCells (:);          (how many cells in this cluster)
%      clusterProps (:,3,:) = clusterArea (:);            (the area of this cluster)
%      clusterProps (:,4,:) = clusterPerimeter (:);       (the length of the perimeter)
%      clusterProps (:,5,:) = clusterPerimeterElements (:); (how many elements does the perimeter exist of)
%
% These matrices are filled up with zero rows to make them all the same
% length (same principle as the M matrix)
%
cellProps = ptPostpro.cellProps;
clusterProps = ptPostpro.clusterProps;

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
    
   % Remove the zero rows from cellProps 
   [notZeroEntryRows, notZeroEntryCols] = find (cellProps (:,:,propCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   cells = cellProps (notZeroEntryRows,:,propCount);
   
   % Remove the zero rows from clusterProps 
   [notZeroEntryRows, notZeroEntryCols] = find (clusterProps (:,:,propCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   clusters = clusterProps (notZeroEntryRows,:,propCount);
   
   % Calculate the amount of all cells per frame
   cellAmount  (iCount) = sum (clusters (:, 2));
   
   % Calculate the amount of clusters per frame (a cluster should contain at
   % least 2 nuclei (and therefore we use > 1)
   clusterAmount (iCount) = size (clusters (find (clusters (:,2) > 1)), 1);
   
   % Calculate the average amount of cells per cluster
   sumCellsInCluster (iCount) = sum (clusters (find (clusters (:,2) > 1), 2));
   if clusterAmount (iCount)
      cellsPerCluster (iCount) = round ( sumCellsInCluster (iCount) / clusterAmount (iCount));
   else
      cellsPerCluster (iCount) = 0;
   end
   
   % Calculate the amount of single cells per frame. This is in principle a
   % cluster with only one nuclei
   singleCellAmount (iCount) = size (clusters (find (clusters (:,2) == 1)), 1);
   
   % Calculate the percentage of single cells
   percentageSingleCells (iCount) = (singleCellAmount (iCount) / ...
                                    (singleCellAmount (iCount) + sumCellsInCluster (iCount))) * 100.0;
                                    
   % Calculate the percentage of clustered cells
   percentageClusteredCells (iCount) = 100.0 - percentageSingleCells (iCount);
   
   % Calculate the average area per cluster
   sumClusterArea (iCount) = sum (clusters (find (clusters (:,2) > 1), 3));
   if clusterAmount (iCount) ~= 0
      areaPerCluster (iCount) = sumClusterArea (iCount) / clusterAmount (iCount);
   else
      areaPerCluster (iCount) = 0;
   end
   
   % Calculate the average area per single cell
   sumSingleCellArea (iCount) = sum (clusters (find (clusters (:,2) == 1), 3));
   if singleCellAmount (iCount) ~= 0
      areaPerSingleCell (iCount) = sumSingleCellArea (iCount) / singleCellAmount (iCount);
   else
      areaPerSingleCell (iCount) = 0;
   end
   
   % Calculate area taken up by single and clustered cells together
   totalAreaAllCells (iCount) = sumSingleCellArea (iCount) + sumClusterArea (iCount);
   percentageAreaAllCells (iCount) = (totalAreaAllCells (iCount) / totalAreaFrame) * 100.0;
   
   % Calculate average perimeter length
   sumClusterPerimeter (iCount) = sum (clusters (find (clusters (:,2) > 1), 4));
   if clusterAmount (iCount) ~= 0 
      perimeterLength (iCount) = sumClusterPerimeter (iCount) / clusterAmount (iCount);
   else
      perimeterLength (iCount) = 0;
   end
   
   % Calculate perimeter / area for clusters
   if sumClusterArea (iCount) ~= 0
      perimeterDivArea (iCount) = sumClusterPerimeter (iCount) / sumClusterArea (iCount);
   else
      perimeterDivArea (iCount) = 0;
   end
end 

if cellClusterPlot
   % Generate single cell and cluster plots if the users requested these
   ptPlotCellClusterStats (imageName, savePath, xAxis, cellAmount, clusterAmount, cellsPerCluster, ...
                           singleCellAmount, percentageSingleCells, percentageClusteredCells);
end   
if areaPlot
   % Generate area plots if the users requested these
   ptPlotAreaStats (imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster, percentageAreaAllCells);
end

if perimeterPlot
   % Generate perimater plots if the users requested these
   ptPlotPerimeterStats (imageName, savePath, xAxis, perimeterLength, perimeterDivArea);
end

if cellClusterPlot | areaPlot | perimeterPlot
   % For all the figures we want to keep the xAxis as well 
   save ('xAxis-CellStats.mat','xAxis');
end
