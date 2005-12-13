function ptPlotCellClusterStats (radioButtons, imageName, SaveDir, xAxis, cellClusterStats, windowSize)
% ptPlotCellClusterStats generates the plots for the cell and cluster statistics
%
% SYNOPSIS       ptPlotCellClusterStats (ptPostpro, imageName, SaveDir, xAxis, cellClusterStats, windowSize)
%
% INPUT          ptPostpro : postpro struct
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                cellClusterStats : struct containing the fields:
%
%                   cellAmount : matrix containing amount of all cells per frame
%                   clusterAmount : matrix containing amount of clusters per frame
%                   cellsPerCluster : matrix containing avg amount of cells per cluster per frame
%                   singleCellAmount : matrix containing amount of single cells per frame
%                   percentageSingleCells: matrix containing the percentage of single cells per frame
%                   percentageClusteredCells: matrix containing the percentage of clustered cells per frame
%                windowSize : size of the averaging window                
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotCellClusterStats  uses { nothing }
%                                  
%                ptPlotCellClusterStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotCellClusterStats
% Andre Kerstens        Sep 04          Complete rewrite of plot functions

% Get all data from struct
cellAmount = cellClusterStats.cellAmount;
clusterAmount = cellClusterStats.clusterAmount;
cellsPerCluster = cellClusterStats.cellsPerCluster;
singleCellAmount = cellClusterStats.singleCellAmount;
percentageSingleCells = cellClusterStats.percentageSingleCells;
percentageClusteredCells = cellClusterStats.percentageClusteredCells;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raCellAmount = movingAverage (cellAmount, windowSize, 'median');
    raClusterAmount = movingAverage (clusterAmount, windowSize, 'median');
    raCellsPerCluster = movingAverage (cellsPerCluster, windowSize, 'median');
    raSingleCellAmount = movingAverage (singleCellAmount, windowSize, 'median');
    raPercentageSingleCells = movingAverage (percentageSingleCells, windowSize, 'median');
    raPercentageClusteredCells = movingAverage (percentageClusteredCells, windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.cellclusterplot_1

    if ~radioButtons.donotshowplots
        
        % Generate the figure and title     
        h_fig = figure('Name', imageName);

        % Draw a subplot showing the amount of cells per frame
        ymax = max (cellAmount) + 1;
        subplot(2,2,1); plot (xAxis, cellAmount); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raCellAmount, 'r'); hold off;
        end
        
        title ('Amount of Cells');
        xlabel ('Frames');
        ylabel ('# of cells');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a subplot showing the amount of clusters per frame
        ymax = max (clusterAmount) + 1;
        subplot(2,2,2); plot (xAxis, clusterAmount); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raClusterAmount, 'r'); hold off;
        end
        
        title ('Amount of Clusters');
        xlabel ('Frames');
        ylabel ('# of clusters');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a subplot showing the amount of cells per cluster
        ymax = max (cellsPerCluster) + 1;
        subplot (2,2,3); plot (xAxis, cellsPerCluster); 
         
        if radioButtons.runningaverage
            hold on; plot (xAxis, raCellsPerCluster, 'r'); hold off;
        end
        
        title ('Average Amount of Cells per Cluster');
        xlabel ('Frames');
        ylabel ('# cells per cluster');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

           % Draw a subplot showing the amount of single cells
        ymax = max (singleCellAmount) + 1;
        subplot (2,2,4); plot (xAxis, singleCellAmount); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raSingleCellAmount, 'r'); hold off;
        end
        
        title ('Amount of Single Cells');
        xlabel ('Frames');
        ylabel ('# of single cells');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
        
        % Save the figures in fig, eps and tif format
        hgsave (h_fig, [SaveDir filesep [imageName '_amountAllCells.fig']]);
        % % print (h_fig, [SaveDir filesep [imageName '_amountAllCells.eps']], '-depsc2', '-tiff');
        % print (h_fig, [SaveDir filesep [imageName '_amountAllCells.tif']], '-dtiff', '-zbuffer'); 
    end  % if ~radioButtons.donotshowplots
    
    % Save MAT files for amount of cells and perc. single cells
    save ([SaveDir filesep imageName '_amountAllCells.mat'],'cellAmount');
    
    % Save CSV files for amount of cells and amount single cells
    csvwrite ([SaveDir filesep imageName '_amountAllCells.csv'], [xAxis ; cellAmount]);
    csvwrite ([SaveDir filesep imageName '_amountSingleCells.csv'], [xAxis ; singleCellAmount]);
end

if radioButtons.cellclusterplot_2

    if ~radioButtons.donotshowplots
            
        % Generate a new figure for percentage single / clustered cels
        h_fig2 = figure('Name', imageName);

        % Draw a subplot showing the percentage of single cells
        ymax = 100.0;   % 100% is the max we can get
        subplot (2,1,1); plot (xAxis, percentageSingleCells); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raPercentageSingleCells, 'r'); hold off;
        end
        
        title ('Percentage Single Cells');
        xlabel ('Frames');
        ylabel ('% single cells');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a subplot showing the percentage of clustered cells
        ymax = 100.0;   % 100% is the max we can get
        subplot (2,1,2); plot (xAxis, percentageClusteredCells); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raPercentageClusteredCells, 'r'); hold off;
        end
        
        title ('Percentage Clustered Cells');
        xlabel ('Frames');
        ylabel ('% clustered cells');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
        
        % Save the figures in fig, eps and tif format
        hgsave (h_fig2, [SaveDir filesep [imageName '_percentageSingleClusteredCells.fig']]);
        % print (h_fig2, [SaveDir filesep [imageName '_percentageSingleClusteredCells.eps']], '-depsc2', '-tiff');
        % print (h_fig2, [SaveDir filesep [imageName '_percentageSingleClusteredCells.tif']], '-dtiff'); 
    end  % if ~radioButtons.donotshowplots
    
    % Save MAT files for amount of cells and perc. single cells
    save ([SaveDir filesep imageName '_percentageSingleCells.mat'],'percentageSingleCells');

    % Save CSV files for amount of cells and perc. single cells
    csvwrite ([SaveDir filesep imageName '_percentageSingleCells.csv'], [xAxis ; percentageSingleCells]);
end
