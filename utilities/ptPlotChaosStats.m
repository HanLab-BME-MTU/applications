function ptPlotChaosStats (radioButtons, imageName, savePath, xAxis, chaosStats, windowSize)
% ptPlotChaosStats plots clustering info from MPMs. 
%
% SYNOPSIS       ptPlotChaosStats (imageName, savePath, xAxisNeigh, chaosStats, windowSize)
%
% INPUT          radioButtons : values of radiobuttons on the gui
%                imageName : Name that will be used as the plot title
%                savePath : directory where plots will be stored
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
ripleyClustering = chaosStats.ripleyClustering;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raRipleyClustering = movingAverage (ripleyClustering, windowSize, 'median');
end

% Here's where the plotting starts
if ~radioButtons.donotshowplots

    % Generate the neighbour change plot (all cells)
    h_fig = figure('Name', imageName);

    % Draw a plot showing average velocity of all cells
    ymax = max(ripleyClustering) + (0.1*max (ripleyClustering));
    plot (xAxis, ripleyClustering); 
        
    if radioButtons.runningaverage
        hold on; plot (xAxis, raRipleyClustering, 'r'); hold off;
    end
        
    title ('Avg Clustering Paremeter Value');
    xlabel ('Frames');
    ylabel ('Avg Clust. Value');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig,[savePath filesep [imageName '_avgRipleyClustering.fig']]);
    print (h_fig, [savePath filesep [imageName '_avgRipleyClustering.eps']],'-depsc2','-tiff');
    print (h_fig, [savePath filesep [imageName '_avgRipleyClustering.tif']],'-dtiff');      
end

% Save CSV files
csvwrite ([savePath filesep imageName '_avgRipleyClustering.csv'], [xAxis ; ripleyClustering]);

