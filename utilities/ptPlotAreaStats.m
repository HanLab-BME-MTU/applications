function ptPlotAreaStats (radioButtons, plotName, SaveDir, xAxis, areaStats, windowSize)
% ptPlotAreaStats generates the plots for the cell and cluster area statistics
%
% SYNOPSIS       ptPlotAreaStats (radioButtons, plotName, SaveDir, xAxis, areaStats, windowSize)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                plotName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                areaStats : struct containing the following fields:
%                   areaPerSingleCell : matrix containing avg area (in um^2) of single cells per frame
%                   areaPerCluster : matrix containing avg area (in um^2) of clusters per frame
%                   totalAreaAllCells : matrix containing avg area (in um^2) of all cells per frame
%                   percentageAreaAllCells : matrix containing percentage of the frame area taken up by cells
%                   areaRatio : the average ratio of area/convex-hull-area for a frame
%                windowSize : size of the averaging window                
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotAreaStats  uses { nothing }
%                                  
%                ptPlotAreaStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotAreaStats
% Andre Kerstens        Jul 04          Added percentage area all cells (against total frame area)
% Andre Kerstens        Jul 04          Added average ratio of area/convex-hull-area
% Andre Kerstens        Sep 04          Complete rewrite of plot functions

% Get data from struct
areaPerSingleCell = areaStats.areaPerSingleCell;
areaPerCluster = areaStats.areaPerCluster;
totalAreaAllCells = areaStats.totalAreaAllCells;
percentageAreaAllCells = areaStats.percentageAreaAllCells;
areaRatio = areaStats.convexHullData;

% Initialize average areaRatio
raAreaRatio = zeros(size(areaRatio));

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAreaPerSingleCell = movingAverage (areaPerSingleCell, windowSize, 'median');
    raAreaPerCluster = movingAverage (areaPerCluster, windowSize, 'median');
    raTotalAreaAllCells = movingAverage (totalAreaAllCells, windowSize, 'median');
    raPercentageAreaAllCells = movingAverage (percentageAreaAllCells, windowSize, 'median');
    raAreaRatio(1,:) = movingAverage (areaRatio(1,:), windowSize, 'median');
    raAreaRatio(2,:) = movingAverage (areaRatio(2,:), windowSize, 'median');
    raAreaRatio(3,:) = movingAverage (areaRatio(3,:), windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.areaplot_2
    
    if ~radioButtons.donotshowplots

        % Generate the figure and title      
        h_fig = figure('Name', plotName);

        % Draw a subplot showing the avg area of a single cell    
        ymax = max (areaPerSingleCell) + 1;
        subplot (2,1,1); plot (xAxis, areaPerSingleCell);
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAreaPerSingleCell, 'r'); hold off;
        end
        
        title ('Average Single Cell Area');
        xlabel ('Frames');
        ylabel ('Avg single cell area (um^2)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a subplot showing the avg area of a cluster
        ymax = max (areaPerCluster) + 1;
        subplot (2,1,2); plot (xAxis, areaPerCluster); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAreaPerCluster, 'r'); hold off;
        end
        
        title ('Average Cluster Area');
        xlabel ('Frames');
        ylabel ('Avg cluster area (um^2)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
    
        % Save the figures in fig, eps and tif format     
        hgsave(h_fig,[SaveDir filesep [plotName '_averageCellAndClusterArea.fig']]);
        % print(h_fig, [SaveDir filesep [plotName '_averageCellAndClusterArea.eps']],'-depsc2','-tiff');
        % print(h_fig, [SaveDir filesep [plotName '_averageCellAndClusterArea.tif']],'-dtiff');
    end  % if ~radioButtons.donotshowplots
end

if radioButtons.areaplot_1

    if ~radioButtons.donotshowplots

        % We only shot plot this one if images are available for this job
        % If images were not available all values in totalAreaAllCells are 0
        if ~isempty(find(totalAreaAllCells))
        
            % Generate the figure and title      
            h_fig2 = figure('Name', plotName);

            % Draw a plot showing the total area taken up by all the cells
            ymax = max (totalAreaAllCells) + 1;  % 100%
            subplot (2,1,1); plot (xAxis, totalAreaAllCells); 

            if radioButtons.runningaverage
                hold on; plot (xAxis, raTotalAreaAllCells, 'r'); hold off;
            end

            title ('Total Area taken up by All Cells');
            xlabel ('Frames');
            ylabel ('Avg cell area (um^2)');
            if length (xAxis) > 1
               axis ([xAxis(1) xAxis(end) 0 ymax]);
            else
               axis ([xAxis(1) xAxis(1)+1 0 ymax]);
            end

            % Draw a plot showing the percentage of the frame area taken up by all the cells
            ymax = 100;  % 100%
            subplot (2,1,2); plot (xAxis, percentageAreaAllCells); 

            if radioButtons.runningaverage
                hold on; plot (xAxis, raPercentageAreaAllCells, 'r'); hold off;
            end

            title ('Percentage of Frame Area taken up by All Cells');
            xlabel ('Frames');
            ylabel ('Percentage All Cells');
            if length (xAxis) > 1
               axis ([xAxis(1) xAxis(end) 0 ymax]);
            else
               axis ([xAxis(1) xAxis(1)+1 0 ymax]);
            end

            % Save the figures in fig, eps and tif format     
            hgsave(h_fig2,[SaveDir filesep [plotName '_areaAllCells.fig']]);
            % print(h_fig2, [SaveDir filesep [plotName '_areaAllCells.eps']],'-depsc2','-tiff');
            % print(h_fig2, [SaveDir filesep [plotName '_areaAllCells.tif']],'-dtiff');
        end
    end  % if ~radioButtons.donotshowplotsend
end

if radioButtons.areaplot_3

    if ~radioButtons.donotshowplots

        % Generate the figure and title      
        h_fig3 = figure('Name', plotName);

        % Draw a plot showing the average cluster area
        ymax = 100;  % 100%
        subplot (2,1,1); plot (xAxis, areaRatio(3,:));
         
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAreaRatio(3,:), 'r'); hold off;
        end
        
        title ('Average Ratio Area / Convex-Hull-Area of Clusters');
        xlabel ('Frames');
        ylabel ('Ratio (%)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a plot showing the average cluster area and convex hull area
        ymax = max ([areaRatio(1,:) areaRatio(2,:)]) + 1;; 
        subplot (2,1,2); plot (xAxis, areaRatio(1,:), 'k');
        hold on, plot (xAxis, areaRatio(2,:), 'g'), hold off;
         
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAreaRatio(1,:), 'r'); hold off;
            hold on; plot (xAxis, raAreaRatio(2,:), 'r'); hold off;
        end
        
        title ('Average Area (red) and Convex Hull Area (green)');
        xlabel ('Frames');
        ylabel ('Area (um^2)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
    
        % Save the figures in fig, eps and tif format     
        hgsave(h_fig3,[SaveDir filesep [plotName '_averageAreaAndRatio.fig']]);
        % print(h_fig3, [SaveDir filesep [plotName '_averageAreaAndRatio.eps']],'-depsc2','-tiff');
        % print(h_fig3, [SaveDir filesep [plotName '_averageAreaAndRatio.tif']],'-dtiff');
    end  % if ~radioButtons.donotshowplots
    
    % Save CSV files
    csvwrite ([SaveDir filesep plotName '_averageAreaAndRatio.csv'], [xAxis ; areaRatio(3,:)]);
end
