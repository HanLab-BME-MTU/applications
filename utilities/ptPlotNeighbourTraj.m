function ptPlotNeighbourTraj (radioButtons, imageName, SaveDir, xAxis, neighTrajStats, windowSize)
% ptPlotNeighbourTraj plots neighbour traj. information gathered in MPM. 
%
% SYNOPSIS       ptPlotNeighbourTraj (imageName, SaveDir, xAxisNeigh, neighTrajStats, windowSize)
%
% INPUT          radioButtons : values of radiobuttons on the gui
%                imageName : Name that will be used as the plot title
%                SaveDir : directory where plots will be stored
%                xAxisNeigh : vector with x-axis values
%                neighTrajStats : vector with neighbourhood trajectories
%                windowSize : size of the averaging window 
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotNeighbourTraj  uses {nothing}
%                                  
%                ptPlotNeighbourTraj is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial version

% Fetch the input data
avgTrajFrame = neighTrajStats.avgTrajFrame;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAvgTrajFrame = movingAverage (avgTrajFrame, windowSize, 'median');
end

% Here's where the plotting starts
if ~radioButtons.donotshowplots

    % Generate the neighbour traj. plot 
    h_fig = figure('Name', imageName);

    % Draw a plot showing neighbour traj. 
    ymax = max (avgTrajFrame) + (0.1*max (avgTrajFrame));
    plot (xAxis, avgTrajFrame, 'c'); 
    
    if radioButtons.plotestimate
       hold on;
       [yPlot, est] = ptPlotEstimate (xAxis, avgTrajFrame, 1);
       hold off;
    end
        
    if radioButtons.runningaverage
        hold on; plot (xAxis, raAvgTrajFrame, 'b'); hold off;
    end
        
    title ('Avg Neighbour Trajectory Velocity');
    xlabel ('Frames');
    ylabel ('Avg Traj. Vel. (um/min)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format        
    hgsave (h_fig,[SaveDir filesep [imageName '_avgNeighbourTrajLength.fig']]);
    print (h_fig, [SaveDir filesep [imageName '_avgNeighbourTrajLength.eps']],'-depsc2','-tiff');
    print (h_fig, [SaveDir filesep [imageName '_avgNeighbourTrajLength.tif']],'-dtiff');      
end

% Save CSV files
csvwrite ([SaveDir filesep imageName '_avgNeighbourTrajLength.csv'], [xAxis ; avgTrajFrame]);

if radioButtons.plotestimate
   if ~radioButtons.donotshowplots 
       csvwrite ([SaveDir filesep imageName '_fittedCurveNeighbourTraj.csv'], [xAxis ; yPlot]);
       csvwrite ([SaveDir filesep imageName '_curveEstimatesNeighbourTraj.csv'], est);
   end
end
