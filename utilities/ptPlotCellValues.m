function ptPlotCellValues (ptPostpro)
% ptPlotCellValues plots information gathered in cellProps and
% clusterProps during the cell tracking (PolyTrack) 
%
% SYNOPSIS       ptPlotCellValues (ptPostpro)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotCellValues.m  uses { nothing }
%                                  
%                ptPlotCellValues.m is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file

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

% Generate single cell and cluster plots if the users requested these
if cellClusterPlot
    
   % Generate the figure and title     
   h_fig = figure; title (imageName);
   
   % Draw a subplot showing the amount of cells per frame
   ymax = max (cellAmount) + 1;
   subplot(2,2,1); plot (xAxis, cellAmount); 
   title ('Amount of Cells');
   xlabel ('Frames');
   ylabel ('# of cells');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
   
   % Draw a subplot showing the amount of clusters per frame
   ymax = max (clusterAmount) + 1;
   subplot(2,2,2); plot (xAxis, clusterAmount); 
   title ('Amount of Clusters');
   xlabel ('Frames');
   ylabel ('# of clusters');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
      
   % Draw a subplot showing the amount of cells per cluster
   ymax = max (cellsPerCluster) + 1;
   subplot (2,2,3); plot (xAxis, cellsPerCluster); 
   title ('Average Amount of Cells per Cluster');
   xlabel ('Frames');
   ylabel ('# cells per cluster');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
 
   % Draw a subplot showing the amount of single cells
   ymax = max (singleCellAmount) + 1;
   subplot (2,2,4); plot (xAxis, singleCellAmount); 
   title ('Amount of Single Cells');
   xlabel ('Frames');
   ylabel ('# of single cells');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
   
   
   % Generate a new figure for percentage single / clustered cels
   h_fig2 = figure; title (imageName);
   
   % Draw a subplot showing the percentage of single cells
   ymax = 100.0;   % 100% is the max we can get
   subplot (2,1,1); plot (xAxis, percentageSingleCells); 
   title ('Percentage Single Cells');
   xlabel ('Frames');
   ylabel ('% single cells');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
   
   % Draw a subplot showing the percentage of clustered cells
   ymax = 100.0;   % 100% is the max we can get
   subplot (2,1,2); plot (xAxis, percentageClusteredCells); 
   title ('Percentage Clustered Cells');
   xlabel ('Frames');
   ylabel ('% clustered cells');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
   
   % Save MAT files for amount of cells and perc. single cells
   cd (savePath);
   save ('amountAllCells.mat','cellAmount');
   save ('percentageSingleCells.mat','percentageSingleCells');
   
   % Save CSV files for amount of cells and perc. single cells
   cd (savePath);
   csvwrite ('amountAllCells.csv','cellAmount');
   csvwrite ('percentageSingleCells.csv','percentageSingleCells');
   
   % Save the figures in fig, eps and tif format
   hgsave (h_fig, [savePath filesep 'singleCellsAndClusterStats.fig']);
   print (h_fig, [savePath filesep 'singleCellsAndClusterStats.eps'], '-depsc2', '-tiff');
   print (h_fig, [savePath filesep 'singleCellsAndClusterStats.tif'], '-dtiff');         
end   

% Generate area plots if the users requested these
if areaPlot
    
   % Generate the figure and title      
   h_fig = figure; title (imageName);
   
   % Draw a subplot showing the avg area of a single cell    
   ymax = max (areaPerSingleCell) + 1;
   subplot (2,1,1); plot (xAxis, areaPerSingleCell);
   title ('Average Single Cell Area');
   xlabel ('Frames');
   ylabel ('Avg single cell area');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
        
   % Draw a subplot showing the avg area of a cluster
   ymax = max (areaPerCluster) + 1;
   subplot (2,1,2); plot (xAxis, areaPerCluster); 
   title ('Average Cluster Area');
   xlabel ('Frames');
   ylabel ('Avg cluster area');
   axis ([xAxis(1) xAxis(end) 0 ymax]);     
   
   % Save the figures in fig, eps and tif format     
   hgsave(h_fig,[savePath filesep 'cellAndClusterAreaStats.fig']);
   print(h_fig, [savePath filesep 'cellAndClusterAreaStats.eps'],'-depsc2','-tiff');
   print(h_fig, [savePath filesep 'cellAndClusterAreaStats.tif'],'-dtiff');
end

if perimeterPlot
        
   % Generate the figure and title
   h_fig = figure; title (imageName);
   
   % Draw a plot showing the avg perimeter length of clusters
   ymax = max (perimeterLength) + 1;
   subplot (2,1,1); plot (xAxis, perimeterLength); 
   title ('Average Perimeter Length of Clusters');
   xlabel ('Frames');
   ylabel ('Perimeter Length');
   axis ([xAxis(1) xAxis(end) 0 ymax]);
   
   % Draw a plot showing the perimeter of clusters divided by area
   ymax = max (perimeterDivArea);
   subplot (2,1,2); plot (xAxis, perimeterDivArea); 
   title ('Perimeter/Area of Clusters');
   xlabel ('Frames');
   ylabel ('Perimeter/Area');
   axis ([xAxis(1) xAxis(end) 0 ymax]);     
   
   % Save the figures in fig, eps and tif format        
   hgsave(h_fig,[savePath filesep 'areaPerimStats.fig']);
   print(h_fig, [savePath filesep 'areaPerimStats.eps'],'-depsc2','-tiff');
   print(h_fig, [savePath filesep 'areaPerimStats.tif'],'-dtiff');     
end
