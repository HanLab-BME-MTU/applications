function ptPlotSpeedStats (radioButtons, imageName, SaveDir, xAxis, avgVelocityStats, windowSize)
% ptPlotSpeedStats generates the plots for the velocity statistics
%
% SYNOPSIS       ptPlotSpeedStats (radioButtons, imageName, SaveDir, xAxis, velocityStats)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                velocityStats : struct with the following fields:
%
%                   avgVelocity : vector with avg velocity all cells
%                   avgVelocitySquared : vector with squared avg velocity all cells
%                   avgSingleVelocity : vector with avg single cell vel.
%                   avgClusteredVelocity : vector with avg clustered cell vel.
%                windowSize : size of the averaging window                
%
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotsSpeedStats  uses { nothing }
%                                  
%                ptPlotSpeedStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

% Get data from struct
avgVelocity = avgVelocityStats.avgVelocity;
avgVelocitySquared = avgVelocityStats.avgVelocitySquared;
avgSingleVelocity = avgVelocityStats.avgSingleVelocity;
avgClusteredVelocity = avgVelocityStats.avgClusteredVelocity;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAvgVelocity = movingAverage (avgVelocity, windowSize, 'median');
    raAvgVelocitySquared = movingAverage (avgVelocitySquared, windowSize, 'median');
    raAvgSingleVelocity = movingAverage (avgSingleVelocity, windowSize, 'median');
    raAvgClusteredVelocity = movingAverage (avgClusteredVelocity, windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.speedplot_2

    if ~radioButtons.donotshowplots
        
        % Generate the avg velocity plot (all cells)
        h_fig = figure('Name', imageName);

        % Draw a plot showing average velocity of all cells
        ymax = max (avgVelocity) + 1;
        subplot (2,1,1); plot (xAxis, avgVelocity); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgVelocity, 1, [2 1 1]);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveAllVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesAllVelocity.csv'], est);
           end
        end
        
        title ('Avg Velocity All Cells');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        ymax = max (avgVelocitySquared) + 1;
        subplot (2,1,2); plot (xAxis, avgVelocitySquared); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgVelocitySquared, 'r'); hold off;
        end
                
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgVelocitySquared, 1, [2 1 2]);
           hold off;
        end        
        
        title ('Avg Squared Velocity All Cells');
        xlabel ('Frames');
        ylabel ('Velocity^2 (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Save the figures in fig, eps and tif format        
        hgsave (h_fig,[SaveDir filesep [imageName '_avgVelocityAllCells.fig']]);
        print (h_fig, [SaveDir filesep [imageName '_avgVelocityAllCells.eps']],'-depsc2','-tiff');
        print (h_fig, [SaveDir filesep [imageName '_avgVelocityAllCells.tif']],'-dtiff');      

        % Generate the figure and title
        h_fig2 = figure('Name', imageName);

        % Draw a subplot showing the avg velocity of a single cell    
        ymax = max (avgSingleVelocity) + 1;
        subplot (2,1,1); plot (xAxis, avgSingleVelocity);
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgSingleVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgSingleVelocity, 1, [2 1 1]);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveSingleVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesSingleVelocity.csv'], est);
           end
        end  
        
        title ('Avg Single Cell Velocity');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end   

        % Draw a subplot showing the avg velocity of a cluster
        ymax = max (avgClusteredVelocity) + 1;
        subplot (2,1,2); plot (xAxis, avgClusteredVelocity); 

        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgClusteredVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgClusteredVelocity, 1, [2 1 2]);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveClusteredVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesClusteredVelocity.csv'], est);
           end
        end  
        
        title ('Avg Clustered Cell Velocity');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
                
        % Save the figures in fig, eps and tif format     
        hgsave (h_fig2,[SaveDir filesep [imageName '_avgSingleAndClusterVelocity.fig']]);
        print (h_fig2, [SaveDir filesep [imageName '_avgSingleAndClusterVelocity.eps']],'-depsc2','-tiff');
        print (h_fig2, [SaveDir filesep [imageName '_avgSingleAndClusterVelocity.tif']],'-dtiff');
    end  % if ~radioButtons.donotshowplots
        
    % Save MAT files for avg all, single and clustered cell velocity
    cd (SaveDir);
    save ([imageName '_avgCellVelocity.mat'],'avgVelocity');
    save ([imageName '_avgSingleCellVelocity.mat'],'avgSingleVelocity');
    save ([imageName '_avgClusteredCellVelocity.mat'],'avgClusteredVelocity');

    % Save CSV files for avg all, single and clustered cell velocity
    csvwrite ([imageName '_avgCellVelocity.csv'], [xAxis ; avgVelocity]);
    csvwrite ([imageName '_avgSingleCellVelocity.csv'], [xAxis ; avgSingleVelocity]);
    csvwrite ([imageName '_avgClusteredCellVelocity.csv'], [xAxis ; avgClusteredVelocity]);
end
