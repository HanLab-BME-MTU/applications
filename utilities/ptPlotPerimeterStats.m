function ptPlotPerimeterStats (imageName, savePath, xAxis, perimeterLength, perimeterDivArea)
% ptPlotPerimeterStats generates the plots for the cell and cluster perimeter statistics
%
% SYNOPSIS       ptPlotPerimeterStats (imageName, savePath, xAxis, perimeterLength, perimeterDivArea)
%
% INPUT          imageName : name of the first image in the movie (used as title)
%                savePath : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                perimeterLength : matrix containing avg cluster perimeter length per frame
%                perimeterDivArea : matrix containing avg cluster perimeter length divided
%                                   by cluster area per frame
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotPerimeterStats  uses { nothing }
%                                  
%                ptPlotPerimeterStats is used by { ptCalculatePlotValues }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotPerimeterStats

% Generate the figure and title
h_fig = figure('Name', imageName);

% Draw a plot showing the avg perimeter length of clusters
ymax = max (perimeterLength) + 1;
subplot (2,1,1); plot (xAxis, perimeterLength); 
title ('Average Perimeter Length of Clusters');
xlabel ('Frames');
ylabel ('Perimeter Length');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end

% Draw a plot showing the perimeter of clusters divided by area
ymax = max (perimeterDivArea);
subplot (2,1,2); plot (xAxis, perimeterDivArea); 
title ('Perimeter/Area of Clusters');
xlabel ('Frames');
ylabel ('Perimeter/Area');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end

% Save the figures in fig, eps and tif format        
hgsave(h_fig,[savePath filesep 'areaPerimStats.fig']);
print(h_fig, [savePath filesep 'areaPerimStats.eps'],'-depsc2','-tiff');
print(h_fig, [savePath filesep 'areaPerimStats.tif'],'-dtiff');     
