function ptPlotCellCellDist (radioButtons, imageName, SaveDir, xAxis, cellCellDistStats, windowSize)
%  ptPlotCellCellDist generates the plots for the inter cell distances
%
% SYNOPSIS       ptPlotCellCellDist (radioButtons, imageName, SaveDir, xAxis, cellCellDistStats, windowSize)
%
% INPUT          ptPostpro : postprocessing structure
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                cellCellDistStats : struct containing:
%                   averageDist : vector with average inter-cell distances
%                windowSize : size of the averaging window                
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotCellCellDist  uses { nothing }
%                                  
%                ptPlotCellCellDist is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial release of ptPlotCellCellDist

% Extract values from the data struct
averageDist = cellCellDistStats.averageDist;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAverageDist = movingAverage (averageDist, windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.cellcelldistplot_1
    
    if ~radioButtons.donotshowplots

        % Generate the figure and title     
        h_fig = figure('Name', imageName);

        % Draw a plot showing the avg distance between nuclei per frame
        ymax = max (averageDist) + 1;
        plot (xAxis, averageDist); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAverageDist, 'r'); hold off;
        end
        
        title ('Average Distance between Cells (triangulated)');
        xlabel ('Frames');
        ylabel ('Avg distance (um)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Save the figures in fig, eps and tif format     
        hgsave (h_fig,[SaveDir filesep [imageName '_AvgDistanceBetweenCells.fig']]);
        print (h_fig, [SaveDir filesep [imageName '_AvgDistanceBetweenCells.eps']],'-depsc2','-tiff');
        print (h_fig, [SaveDir filesep [imageName '_AvgDistanceBetweenCells.tif']],'-dtiff');
    end  % if ~radioButtons.donotshowplots    
end
