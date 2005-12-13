function ptPlotSingleSpeedStats (radioButtons, imageName, SaveDir, xAxis, velocitySingleStats, windowSize)
% ptPlotSingleSpeedStats generates the plots for the velocity statistics
%
% SYNOPSIS       ptPlotSingleSpeedStats (radioButtons, imageName, SaveDir, xAxis, velocityStats)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                velocityStats : struct with the following fields:
%                   velAllCellsHigherThanAvgSingleCells : vector with
%                        amount of all cells with velocity higher than avg single cells
%                windowSize : size of the averaging window  
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotsSingleSpeedStats  uses { nothing }
%                                  
%                ptPlotSingleSpeedStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

% Get data from struct
velAllCellsHigherThanAvgSingleCells = velocitySingleStats.velAllCellsHigherThanAvgSingleCells;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raVelAllCellsHigherThanAvgSingleCells = movingAverage (velAllCellsHigherThanAvgSingleCells, windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.speedplot_1

    if ~radioButtons.donotshowplots

        % Generate the plot that shows all cells velocity against single cell
        % velocity (in percent)
        h_fig3 = figure('Name', imageName);

        % Draw the plot
        ymax = 100;
        plot (xAxis, velAllCellsHigherThanAvgSingleCells); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raVelAllCellsHigherThanAvgSingleCells, 'r'); hold off;
        end
        
        title ('% All Cells with Velocity higher than Average Single Cell Velocity');
        xlabel ('Frames');
        ylabel ('Percentage (%)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Save the figures in fig, eps and tif format        
        hgsave (h_fig3,[SaveDir filesep [imageName '_velAllCellsHigherThanAvgSingleCells.fig']]);
        % print (h_fig3, [SaveDir filesep [imageName '_velAllCellsHigherThanAvgSingleCells.eps']],'-depsc2','-tiff');
        % print (h_fig3, [SaveDir filesep [imageName '_velAllCellsHigherThanAvgSingleCells.tif']],'-dtiff');
    end  % if ~radioButtons.donotshowplots
    
    % Save MAT files for percentage velocity higher than single cell speed
    save ([SaveDir filesep imageName '_velAllCellsHigherThanAvgSingleCells.mat'],'velAllCellsHigherThanAvgSingleCells');

    % Save CSV files for percentage velocity higher than single cell speed
    csvwrite ([SaveDir filesep imageName '_velAllCellsHigherThanAvgSingleCells.csv'], [xAxis ; velAllCellsHigherThanAvgSingleCells]);
end
