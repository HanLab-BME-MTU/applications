function ptPlotPerimeterStats (radioButtons, imageName, SaveDir, xAxis, perimeterStats, windowSize)
% ptPlotPerimeterStats generates the plots for the cell and cluster perimeter statistics
%
% SYNOPSIS       ptPlotPerimeterStats (radioButtons, imageName, SaveDir, xAxis, perimeterStats, windowSize)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                perimeterStats : struct containing the following fields:
%                   perimeterLength : matrix containing avg cluster perimeter length per frame
%                   perimeterDivArea : matrix containing avg cluster perimeter length divided
%                                      by cluster area per frame
%                windowSize : size of the averaging window                
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotPerimeterStats  uses { nothing }
%                                  
%                ptPlotPerimeterStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotPerimeterStats
% Andre Kerstens        Sep 04          Complete rewrite of plot functions

% Get data from struct
perimeterLength = perimeterStats.perimeterLength;
perimeterDivArea = perimeterStats.perimeterDivArea;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raPerimeterLength = movingAverage (perimeterLength, windowSize, 'median');
    raPerimeterDivArea = movingAverage (perimeterDivArea, windowSize, 'median');
end

% Here is where all the plotting starts
if ~radioButtons.donotshowplots

    % Generate the figure and title
    h_fig = figure('Name', imageName);

    % Draw a plot showing the avg perimeter length of clusters
    ymax = max (perimeterLength) + 1;
    subplot (2,1,1); plot (xAxis, perimeterLength); 
        
    if radioButtons.runningaverage
       hold on; plot (xAxis, raPerimeterLength, 'r'); hold off;
    end
        
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
        
    if radioButtons.runningaverage
       hold on; plot (xAxis, raPerimeterDivArea, 'r'); hold off;
    end
        
    title ('Perimeter/Area of Clusters');
    xlabel ('Frames');
    ylabel ('Perimeter/Area');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format        
    hgsave(h_fig,[SaveDir filesep [imageName '_areaPerimStats.fig']]);
    % print(h_fig, [SaveDir filesep [imageName '_areaPerimStats.eps']],'-depsc2','-tiff');
    % print(h_fig, [SaveDir filesep [imageName '_areaPerimStats.tif']],'-dtiff');     
end  % if ~radioButtons.donotshowplots
