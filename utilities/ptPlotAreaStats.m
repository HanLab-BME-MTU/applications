function ptPlotAreaStats (imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster, percentageAreaAllCells)
% ptPlotAreaStats generates the plots for the cell and cluster area statistics
%
% SYNOPSIS       ptPlotAreaStats (imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster, percentageAreaAllCells)
%
% INPUT          imageName : name of the first image in the movie (used as title)
%                savePath : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                areaPerSingleCell : matrix containing avg area (in pixels) of single cells per frame
%                areaPerCluster : matrix containing avg area (in pixels) of clusters per frame
%                percentageAreaAllCells : matrix containing percentage of
%                                         the frame area taken up by cells
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotAreaStats  uses { nothing }
%                                  
%                ptPlotAreaStats is used by { ptCalculatePlotValues }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotAreaStats
% Andre Kerstens        Jul 04          Added percentage area all cells (against total frame area)

% Generate the figure and title      
h_fig = figure('Name', imageName);

% Draw a subplot showing the avg area of a single cell    
ymax = max (areaPerSingleCell) + 1;
subplot (2,1,1); plot (xAxis, areaPerSingleCell);
title ('Average Single Cell Area');
xlabel ('Frames');
ylabel ('Avg single cell area');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end
     
% Draw a subplot showing the avg area of a cluster
ymax = max (areaPerCluster) + 1;
subplot (2,1,2); plot (xAxis, areaPerCluster); 
title ('Average Cluster Area');
xlabel ('Frames');
ylabel ('Avg cluster area');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end

% Save the figures in fig, eps and tif format     
hgsave(h_fig,[savePath filesep 'cellAndClusterAreaStats.fig']);
print(h_fig, [savePath filesep 'cellAndClusterAreaStats.eps'],'-depsc2','-tiff');
print(h_fig, [savePath filesep 'cellAndClusterAreaStats.tif'],'-dtiff');

% Generate the figure and title      
h_fig = figure('Name', imageName);

% Draw a plot showing the percentage of the frame area taken up by all the cells
ymax = 100;  % 100%
plot (xAxis, percentageAreaAllCells); 
title ('Percentage of Frame Area taken up by All Cells');
xlabel ('Frames');
ylabel ('Percentage All Cells');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end

% Save the figures in fig, eps and tif format     
hgsave(h_fig,[savePath filesep 'percentageAreaAllCells.fig']);
print(h_fig, [savePath filesep 'percentageAreaAllCells.eps'],'-depsc2','-tiff');
print(h_fig, [savePath filesep 'percentageAreaAllCells.tif'],'-dtiff');
