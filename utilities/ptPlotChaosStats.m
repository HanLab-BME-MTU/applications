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
ripleyClustering = chaosStats.ripleyClustering;

if radioButtons.ripleysimplot
   simulationAvg = chaosStats.simulationAvg;
   simulationMax = chaosStats.simulationMax;
   simulationMin = chaosStats.simulationMin;
end

% Calculate average values if we have to
if radioButtons.runningaverage    
    raRipleyClustering = movingAverage (ripleyClustering, windowSize, 'median');
end

% Here's where the plotting starts
if ~radioButtons.donotshowplots

    % Generate the neighbour change plot (all cells)
    h_fig = figure('Name', imageName);

    % Draw a plot showing average velocity of all cells
    ymax1 = max(ripleyClustering) + (0.1*max (ripleyClustering));
    if radioButtons.ripleysimplot
       ymax2 = max(simulationMax) + (0.1*max (simulationMax));
    else
       ymax2 = 0;
    end
    ymax = max([ymax1 ymax2]);

    % Plot the running average
    if radioButtons.runningaverage
        plot (xAxis, ripleyClustering, 'c', xAxis, raRipleyClustering, 'r'); 
        legend('Clust value','Avg clust value',3);
    else
        plot (xAxis, ripleyClustering, 'c'); 
        legend('Clust value',3);
    end
    
    % Plot the simulation results as well if the user wants this
    if radioButtons.ripleysimplot
       hold on; 
       plot (xAxis, simulationAvg, 'g');
       plot (xAxis, simulationMax, '--');
       plot (xAxis, simulationMin, '.');  
       if radioButtons.runningaverage
          legend('Clust value','Avg clust value','Simulated avg','Simulated max','Simulated min',3);
       else
          legend('Clust value','Simulated avg','Simulated max','Simulated min',3);
       end
       hold off;
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
    hgsave (h_fig,[SaveDir filesep [imageName '_avgRipleyClustering.fig']]);
    print (h_fig, [SaveDir filesep [imageName '_avgRipleyClustering.eps']],'-depsc2','-tiff');
    print (h_fig, [SaveDir filesep [imageName '_avgRipleyClustering.tif']],'-dtiff');      
end

% Save CSV files
csvwrite ([SaveDir filesep imageName '_avgRipleyClustering.csv'], [xAxis ; ripleyClustering]);

if radioButtons.ripleysimplot
    csvwrite ([SaveDir filesep imageName '_avgRipleySimulationValues.csv'], ...
              [xAxis ; simulationAvg ; simulationMax ; simulationMin]);
end

