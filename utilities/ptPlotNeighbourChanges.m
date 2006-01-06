function ptPlotNeighbourChanges (radioButtons, imageName, SaveDir, xAxis, neighChangeStats, windowSize, drugTimepoint)
% ptPlotNeighbourChanges plots neighbour change info from MPM. 
%
% SYNOPSIS       ptPlotNeighbourChanges (imageName, SaveDir, xAxisNeigh, neighChangeStats, windowSize)
%
% INPUT          radioButtons : values of radiobuttons on the gui
%                imageName : Name that will be used as the plot title
%                SaveDir : directory where plots will be stored
%                xAxisNeigh : vector with x-axis values
%                neighChangeStats : vector with neighbourhood interactions
%                windowSize : size of the averaging window 
%                drugTimepoint : frame number where the drug/EGF is applied 
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotNeighbourChanges  uses {nothing}
%                                  
%                ptPlotNeighbourhoodChanges is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial version

% Fetch the input data
avgNbChange = neighChangeStats.avgNbChange;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAvgNbChange = movingAverage (avgNbChange, windowSize, 'median');
end

% Here's where the plotting starts
if ~radioButtons.donotshowplots

    % Generate the neighbour change plot (all cells)
    h_fig2 = figure('Name', imageName);

    % Draw a plot showing neighborhood change
    ymax = max (avgNbChange) + (0.1*max (avgNbChange));
    plot (xAxis, avgNbChange, 'b'); 
        
    % If needed show the fitted trapezoid on the plot
    if radioButtons.plotestimate
       hold on;
       [yPlot, est] = ptPlotEstimate (xAxis, avgNbChange, 1, [], drugTimepoint);
       hold off;
    end
    
    if radioButtons.runningaverage
        hold on; plot (xAxis, raAvgNbChange, 'c'); hold off;
    end
        
    title ('Avg Neighbour Interaction Change');
    xlabel ('Frames');
    ylabel ('Avg Neighbour Change');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig2,[SaveDir filesep [imageName '_avgNeighbourInteractionChange.fig']]);
    print (h_fig2, [SaveDir filesep [imageName '_avgNeighbourInteractionChange.eps']],'-depsc2','-tiff');
    print (h_fig2, [SaveDir filesep [imageName '_avgNeighbourInteractionChange.tif']],'-dtiff'); 
end

% Save CSV files
csvwrite ([SaveDir filesep imageName '_avgNeighbourInteractionChange.csv'], [xAxis ; avgNbChange]);

if radioButtons.plotestimate
    if ~radioButtons.donotshowplots 
        if ~isempty(yPlot)
             csvwrite ([SaveDir filesep imageName '_fittedCurveNeighbourChange.csv'], [xAxis ; yPlot]);
        end
        if ~isempty(est)
             csvwrite ([SaveDir filesep imageName '_curveEstimatesNeighbourChange.csv'], est);
        end
    end
end

