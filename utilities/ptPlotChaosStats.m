function ptPlotChaosStats (radioButtons, imageName, SaveDir, xAxis, chaosStats, windowSize)
% ptPlotChaosStats plots clustering info from MPMs. 
%
% SYNOPSIS       ptPlotChaosStats (imageName, SaveDir, xAxisNeigh, chaosStats, windowSize)
%
% INPUT          radioButtons : values of radiobuttons on the gui
%                imageName : Name that will be used as the plot title
%                SaveDir : directory where plots will be stored
%                xAxisNeigh : vector with x-axis values
%                chaosStats : vector with clustering info
%                windowSize : size of the averaging window 
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotChaosStats uses { nothing }
%                                  
%                ptPlotChaosStats is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version

% Fetch the input data
ripleySlopeInclin = chaosStats.ripleySlopeInclin;
ripleySlopeStart = chaosStats.ripleySlopeStart;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raRipleySlopeInclin = movingAverage (ripleySlopeInclin, windowSize, 'median');
end

% Calculate average values if we have to
if radioButtons.runningaverage    
    raRipleySlopeStart = movingAverage (ripleySlopeStart, windowSize, 'median');
end

% Here's where the plotting starts
if ~radioButtons.donotshowplots

    % Generate the chaos plot
    h_fig = figure('Name', imageName);

    % Draw a plot showing clustering parameter
    ymax = max(ripleySlopeInclin) + (0.1*max (ripleySlopeInclin));
   
    % Plot the running average
    if radioButtons.runningaverage
        plot (xAxis, ripleySlopeInclin, 'c', xAxis, raRipleySlopeInclin, 'b'); 
        legend('Inclination value','Avg inclination value',1);
    else
        plot (xAxis, ripleySlopeInclin, 'c'); 
        legend('Inclination value',1);
    end

    % If needed show the fitted trapezoid on the plot
    if radioButtons.plotestimate
       hold on;
       [yPlot, est] = ptPlotEstimate (xAxis, ripleySlopeInclin, -1);
       hold off;
    end
            
    title ('Avg Inclination Value');
    xlabel ('Frames');
    ylabel ('Avg Inclin. Value');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig,[SaveDir filesep [imageName '_avgRipleySlopeInclin.fig']]);
    print (h_fig, [SaveDir filesep [imageName '_avgRipleySlopeInclin.eps']],'-depsc2','-tiff');
    print (h_fig, [SaveDir filesep [imageName '_avgRipleySlopeInclin.tif']],'-dtiff');      

    
    % Generate the chaos plot
    h_fig2 = figure('Name', imageName);

    % Draw a plot showing clustering parameter
    ymax = max(ripleySlopeStart) + (0.1*max (ripleySlopeStart));
   
    % Plot the running average
    if radioButtons.runningaverage
        plot (xAxis, ripleySlopeStart, 'c', xAxis, raRipleySlopeStart, 'b'); 
        legend('Slope start value','Avg slope start value',2);
    else
        plot (xAxis, ripleySlopeStart, 'c'); 
        legend('Slope start value',2);
    end
            
    title ('Avg Slope Start Value');
    xlabel ('Frames');
    ylabel ('Avg Slope Start Value');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig2,[SaveDir filesep [imageName '_avgRipleySlopeStart.fig']]);
    print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.eps']],'-depsc2','-tiff');
    print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.tif']],'-dtiff');      
end

% Save CSV files
csvwrite ([SaveDir filesep imageName '_avgRipleySlopeInclination.csv'], [xAxis ; ripleySlopeInclin]);
csvwrite ([SaveDir filesep imageName '_avgRipleySlopeStart.csv'], [xAxis ; ripleySlopeStart]);

