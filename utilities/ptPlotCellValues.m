function ptPlotCellValues (ptPostpro, MPM)
% ptPlotCellValues plots information gathered in cellProps and
% clusterProps during the cell tracking (PolyTrack) 
%
% SYNOPSIS       ptPlotCellValues (ptPostpro)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
%                MPM       : matrix containing the cell tracks
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
% Andre Kerstens        May 04          Cleaned up source and renamed file

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.plotfirstimg;
endFrame = ptPostpro.plotlastimg;
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

% Calculate a number of statistics for every frame
for frameCount = startFrame : endFrame
   
   % Remove the zero rows from cellProps 
   [notZeroEntryRows, notZeroEntryCols] = find (cellProps (:,:,frameCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   cells = cellProps (notZeroEntryRows,:,frameCount);
   
   % Remove the zero rows from clusterProps 
   [notZeroEntryRows, notZeroEntryCols] = find (clusterProps (:,:,frameCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   clusters = clusterProps (notZeroEntryRows,:,frameCount);
   
   % Calculate the amount of clusters per frame (a cluster should contain at
   % least 2 nuclei (and therefore we use > 1)
   clusterAmount (frameCount) = size (clusters (find (clusters (:,2) > 1)), 1);
   
   % Calculate the average amount of cells per cluster
   sumCellsInCluster (frameCount) = sum (clusters (find (clusters (:,2) > 1), 2));
   if clusterAmount (frameCount)
      cellsPerCluster (frameCount) = round ( sumCellsInCluster (frameCount) / clusterAmount (frameCount));
   else
      cellsPerCluster (frameCount) = 0;
   end
   
   % Calculate the amount of single cells per frame. This is in principle a
   % cluster with only one nuclei
   singleCellAmount (frameCount) = size (clusters (find (clusters (:,2) == 1)), 1);
   
   % Calculate the percentage of single cells
   percentageSingleCells (frameCount) = (singleCellAmount (frameCount) / ...
                                        (singleCellAmount (frameCount) + sumCellsInCluster (frameCount))) * 100.0;
   
   % Calculate the average area per cluster
   sumClusterArea (frameCount) = sum (clusters (find (clusters (:,2) > 1), 3));
   if clusterAmount (frameCount) ~= 0
      areaPerCluster (frameCount) = sumClusterArea (frameCount) / clusterAmount (frameCount);
   else
      areaPerCluster (frameCount) = 0;
   end
   
   % Calculate the average area per single cell
   sumSingleCellArea (frameCount) = sum (clusters (find (clusters (:,2) == 1), 3));
   if singleCellAmount (frameCount) ~= 0
      areaPerSingleCell (frameCount) = sumSingleCellArea (frameCount) / singleCellAmount (frameCount);
   else
      areaPerSingleCell (frameCount) = 0;
   end
   
   % Calculate perimeter / area for clusters
   sumClusterPerimeter (frameCount) = sum (clusters (find (clusters (:,2) > 1), 4));
   if sumClusterPerimeter (frameCount) ~= 0
      perimeterDivArea (frameCount) = sumClusterPerimeter (frameCount) / sumClusterArea (frameCount);
   else
      perimeterDivArea (frameCount) = 0;
   end
end 

% Generate single cell and cluster plots if the users requested these
if cellClusterPlot
    
   % Generate the figure and title     
   h_fig = figure; title (imageName);
        
   % Draw a subplot showing the amount of clusters per frame
   ymax = max (clusterAmount) + 1;
   subplot(2,2,1); plot (clusterAmount); 
   title ('Amount of Clusters');
   xlabel ('Frames');
   ylabel ('# of clusters');
   axis ([0 endFrame 0 ymax]);
      
   % Draw a subplot showing the amount of cells per cluster
   ymax = max (cellsPerCluster) + 1;
   subplot (2,2,2); plot (cellsPerCluster); 
   title ('Average Amount of Cells per Cluster');
   xlabel ('Frames');
   ylabel ('# cells per cluster');
   axis ([0 endFrame 0 ymax]);
 
   % Draw a subplot showing the amount of single cells
   ymax = max (singleCellAmount) + 1;
   subplot (2,2,3); plot (singleCellAmount); 
   title ('Amount of Single Cells');
   xlabel ('Frames');
   ylabel ('# of single cells');
   axis ([0 endFrame 0 ymax]);
   
   % Draw a subplot showing the percentage of single cells
   ymax = 100.0;   % 100% is the max we can get
   subplot (2,2,4); plot (percentageSingleCells); 
   title ('Percentage Single Cells');
   xlabel ('Frames');
   ylabel ('% single cells');
   axis ([0 endFrame 0 ymax]);
   
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
   subplot (2,1,1); plot (areaPerSingleCell);
   title ('Average Single Cell Area');
   xlabel ('Frames');
   ylabel ('Avg single cell area');
   axis ([0 endFrame 0 ymax]);
        
   % Draw a subplot showing the avg area of a cluster
   ymax = max (areaPerCluster) + 1;
   subplot (2,1,2); plot (areaPerCluster); 
   title ('Average Cluster Area');
   xlabel ('Frames');
   ylabel ('Avg cluster area');
   axis ([0 endFrame 0 ymax]);
        
   % Save the figures in fig, eps and tif format     
   hgsave(h_fig,[savePath filesep 'cellAndClusterAreaStats.fig']);
   print(h_fig, [savePath filesep 'cellAndClusterAreaStats.eps'],'-depsc2','-tiff');
   print(h_fig, [savePath filesep 'cellAndClusterAreaStats.tif'],'-dtiff');
end

if perimeterPlot
        
   % Generate the figure and title
   h_fig = figure; title (imageName);
   
   % Draw a plot showing the perimeter of clusters divided by area
   ymax = max (perimeterDivArea);
   plot (perimeterDivArea); 
   title ('Perimeter/Area of Clusters');
   xlabel ('Frames');
   ylabel ('Perimeter/Area');
   axis ([0 endFrame 0 ymax]);
        
   % Save the figures in fig, eps and tif format        
   hgsave(h_fig,[savePath filesep 'areasPerim.fig']);
   print(h_fig, [savePath filesep 'areasPerim.eps'],'-depsc2','-tiff');
   print(h_fig, [savePath filesep 'areasPerim.tif'],'-dtiff');     
end
    

% h_fig=figure,title(handles.jobvalues.imagename);
% ymax=max(percentageSingle);
% subplot(1,1,1); plot(percentageSingle), title('Percentage of single cells (red - filtered, blue - non-filtered)');
% hold on;
% plot(percentageSingleFiltered, 'r');
% xlabel('Frames');
% ylabel('percentage of single cells');
% axis([0 endFrame 0 1]);
% hold off;
%         
%            
% hgsave(h_fig,[savePath filesep 'percentSingle.fig']);
% print(h_fig, [savePath filesep 'percentSingle.eps'],'-depsc2','-tiff');
% print(h_fig, [savePath filesep 'percentSingle.tif'],'-dtiff');

% cd(savePath);
% save ('newProps', 'newProps');
% save ('percentageSingle','percentageSingle');
