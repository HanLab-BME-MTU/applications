function ptPlotCellClusterStats (ptPostpro, imageName, savePath, xAxis, cellAmount, clusterAmount, cellsPerCluster, singleCellAmount, percentageSingleCells, percentageClusteredCells)
% ptPlotCellClusterStats generates the plots for the cell and cluster statistics
%
% SYNOPSIS       ptPlotCellClusterStats (ptPostpro, imageName, savePath, xAxis, cellAmount, 
%                                        clusterAmount, cellsPerCluster, singleCellAmount, 
%                                        percentageSingleCells, percentageClusteredCells)
%
% INPUT          ptPostpro : postpro struct
%                imageName : name of the first image in the movie (used as title)
%                savePath : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                cellAmount : matrix containing amount of all cells per frame
%                clusterAmount : matrix containing amount of clusters per frame
%                cellsPerCluster : matrix containing avg amount of cells per cluster per frame
%                singleCellAmount : matrix containing amount of single cells per frame
%                percentageSingleCells: matrix containing the percentage of single cells per frame
%                percentageClusteredCells: matrix containing the percentage of clustered cells per frame
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotCellClusterStats  uses { nothing }
%                                  
%                ptPlotCellClusterStats is used by { ptPlotCellValues }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 04          Initial release of ptPlotCellClusterStats

if ptPostpro.cellclusterplot_1

    % Generate the figure and title     
    h_fig = figure('Name', imageName);

    % Draw a subplot showing the amount of cells per frame
    ymax = max (cellAmount) + 1;
    subplot(2,2,1); plot (xAxis, cellAmount); 
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
    title ('Amount of Single Cells');
    xlabel ('Frames');
    ylabel ('# of single cells');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end
    
    % Save MAT files for amount of cells and perc. single cells
    cd (savePath);
    save ('amountAllCells.mat','cellAmount');
    
    % Save CSV files for amount of cells and perc. single cells
    cd (savePath);
    csvwrite ('amountAllCells.csv', [xAxis ; cellAmount]);
    
    % Save the figures in fig, eps and tif format
    hgsave (h_fig, [savePath filesep 'amountAllCells.fig']);
    print (h_fig, [savePath filesep 'amountAllCells.eps'], '-depsc2', '-tiff');
    print (h_fig, [savePath filesep 'amountAllCells.tif'], '-dtiff'); 
end

if ptPostpro.cellclusterplot_2
    
    % Generate a new figure for percentage single / clustered cels
    h_fig2 = figure('Name', imageName);

    % Draw a subplot showing the percentage of single cells
    ymax = 100.0;   % 100% is the max we can get
    subplot (2,1,1); plot (xAxis, percentageSingleCells); 
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
    title ('Percentage Clustered Cells');
    xlabel ('Frames');
    ylabel ('% clustered cells');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save MAT files for amount of cells and perc. single cells
    cd (savePath);
    save ('percentageSingleCells.mat','percentageSingleCells');

    % Save CSV files for amount of cells and perc. single cells
    cd (savePath);
    csvwrite ('percentageSingleCells.csv', [xAxis ; percentageSingleCells]);

    % Save the figures in fig, eps and tif format
    hgsave (h_fig2, [savePath filesep 'percentageSingleClusteredCells.fig']);
    print (h_fig2, [savePath filesep 'percentageSingleClusteredCells.eps'], '-depsc2', '-tiff');
    print (h_fig2, [savePath filesep 'percentageSingleClusteredCells.tif'], '-dtiff'); 
end
