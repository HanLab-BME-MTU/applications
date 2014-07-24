function ptCalculateCellStatistics (ptPostpro)
% ptCalculateCellStatistics calculates statistics information gathered in cellProps and
% clusterProps during the cell tracking (PolyTrack) 
%
% SYNOPSIS       ptCalculateCellStatistics (ptPostpro)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
%                
% OUTPUT         None (values are saved on disk) 
%
% DEPENDENCIES   ptCalculateCellStatistics.m  uses { nothing }
%                                  
%                ptCalculateCellStatistics.m is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Created new m-file that calculates
%                                       the statistics before they are
%                                       plotted

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
SaveDir = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;

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

% Save MAT files for amount of cells, amount of clusters, cells per
% cluster, amount of single cells, amount of cells in clusters, perc.
% single cells and perc. clustered cells
cd (SaveDir);
save ('amountAllCells.mat','cellAmount');
save ('amountClusters.mat','clusterAmount');
save ('amountSingleCells.mat','singleCellAmount');
save ('cellsPerCluster.mat','cellsPerCluster');
save ('amountClusteredCells.mat','sumCellsInCluster');
save ('percentageClusteredCells.mat','percentageClusteredCells');
save ('percentageSingleCells.mat','percentageSingleCells');
save ('xAxis-CellStats.mat','xAxis');

% Save CSV files for amount of cells and perc. single cells
cd (SaveDir);
csvwrite ('amountAllCells.csv', [xAxis ; cellAmount]);
csvwrite ('percentageSingleCells.csv', [xAxis ; percentageSingleCells]);

cellAmount, clusterAmount, cellsPerCluster, singleCellAmount, percentageSingleCells, percentageClusteredCells, sumCellsInCluster
cellAndClusterAreaStats
ptPlotPerimeterStats


save ('xAxis-CellStats.mat','xAxis');
