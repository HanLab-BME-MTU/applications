function ptPlotAreaStats (ptPostpro, imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster, totalAreaAllCells, percentageAreaAllCells, areaRatio)
% ptPlotAreaStats generates the plots for the cell and cluster area statistics
%
% SYNOPSIS       ptPlotAreaStats (ptPostpro, imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster, totalAreaAllCells,
%                                 percentageAreaAllCells, areaRatio)
%
% INPUT          ptPostpro : postprocessing structure
%                imageName : name of the first image in the movie (used as title)
%                savePath : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                areaPerSingleCell : matrix containing avg area (in um^2) of single cells per frame
%                areaPerCluster : matrix containing avg area (in um^2) of clusters per frame
%                totalAreaAllCells : matrix containing avg area (in um^2) of all cells per frame
%                percentageAreaAllCells : matrix containing percentage of the frame area taken up by cells
%                areaRatio : the average ratio of area/convex-hull-area for a frame
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
% Andre Kerstens        Jul 04          Added average ratio of area/convex-hull-area

if ptPostpro.areaplot_2

    % Generate the figure and title      
    h_fig = figure('Name', imageName);

    % Draw a subplot showing the avg area of a single cell    
    ymax = max (areaPerSingleCell) + 1;
    subplot (2,1,1); plot (xAxis, areaPerSingleCell);
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
    title ('Average Cluster Area');
    xlabel ('Frames');
    ylabel ('Avg cluster area (um^2)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format     
    hgsave(h_fig,[savePath filesep 'averageCellAndClusterArea.fig']);
    print(h_fig, [savePath filesep 'averageCellAndClusterArea.eps'],'-depsc2','-tiff');
    print(h_fig, [savePath filesep 'averageCellAndClusterArea.tif'],'-dtiff');
end

if ptPostpro.areaplot_1

    % Generate the figure and title      
    h_fig2 = figure('Name', imageName);

    % Draw a plot showing the total area taken up by all the cells
    ymax = max (totalAreaAllCells) + 1;  % 100%
    subplot (2,1,1); plot (xAxis, totalAreaAllCells); 
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
    title ('Percentage of Frame Area taken up by All Cells');
    xlabel ('Frames');
    ylabel ('Percentage All Cells');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format     
    hgsave(h_fig2,[savePath filesep 'areaAllCells.fig']);
    print(h_fig2, [savePath filesep 'areaAllCells.eps'],'-depsc2','-tiff');
    print(h_fig2, [savePath filesep 'areaAllCells.tif'],'-dtiff');
end

if ptPostpro.areaplot_3

    % Generate the figure and title      
    h_fig3 = figure('Name', imageName);

    % Draw a plot showing the average cluster area
    ymax = 100;  % 100%
    subplot (2,1,1); plot (xAxis, areaRatio(:,3));
    title ('Average Ratio Area / Convex-Hull-Area of Clusters');
    xlabel ('Frames');
    ylabel ('Ratio (%)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Draw a plot showing the average cluster area and convex hull area
    ymax = max ([areaRatio(:,1); areaRatio(:,2)]) + 1;; 
    subplot (2,1,2); plot (xAxis, areaRatio(:,1), 'r');
    hold on, plot (xAxis, areaRatio(:,2), 'g'), hold off;
    title ('Average Area (red) and Convex Hull Area (green)');
    xlabel ('Frames');
    ylabel ('Area (um^2)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format     
    hgsave(h_fig3,[savePath filesep 'averageAreaAndRatio.fig']);
    print(h_fig3, [savePath filesep 'averageAreaAndRatio.eps'],'-depsc2','-tiff');
    print(h_fig3, [savePath filesep 'averageAreaAndRatio.tif'],'-dtiff');
    
    % Save CSV files
    csvwrite ('averageAreaAndRatio.csv', [xAxis ; averageAreaAndRatio]);
end