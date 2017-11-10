function ptPlotChaosStats (radioButtons, imageName, SaveDir, xAxis, chaosStats, windowSize, drugTimepoint)
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
%                drugTimepoint : frame number where the drug/EGF is applied 
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
% Johan de Rooij        jun 05          to match new ptCalculateChaosStats

% Fetch the input data
ripleySlopeInclin = chaosStats.ripleySlopeInclin;
ripleySlopeStart = chaosStats.ripleySlopeStart;
AvgDRS = chaosStats.AvgDRS;
xAxisStart = xAxis.Start;
xAxisSlope = xAxis.Slope;
xAxisDer = xAxis.Der;

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
    ymin = min(ripleySlopeInclin) - (0.1*max (ripleySlopeInclin));
   
    % Plot the running average
    if radioButtons.runningaverage
        plot (xAxisSlope, ripleySlopeInclin, 'b', xAxisSlope, raRipleySlopeInclin, 'c'); 
        legend('Inclination value','Avg inclination value',1);
    else
        plot (xAxisSlope, ripleySlopeInclin, 'b'); 
        legend('Inclination value',1);
    end

    % If needed show the fitted trapezoid on the plot
    if radioButtons.plotestimate
       hold on;
       [yPlot1, est1] = ptPlotEstimate (xAxisSlope, ripleySlopeInclin, -1, [], drugTimepoint);
       hold off;
    end
            
    title ('Avg Inclination Value');
    xlabel ('Frames');
    ylabel ('Avg Inclin. Value');
    if length (xAxisSlope) > 1
       axis ([xAxisSlope(1) xAxisSlope(end) 0 ymax]);
    else
       axis ([xAxisSlope(1) xAxisSlope(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig,[SaveDir filesep [imageName '_avgRipleySlopeInclin.fig']]);
    % print (h_fig, [SaveDir filesep [imageName '_avgRipleySlopeInclin.eps']],'-depsc2','-tiff');
    % print (h_fig, [SaveDir filesep [imageName '_avgRipleySlopeInclin.tif']],'-dtiff');      

    
    % Generate the chaos plot
    h_fig2 = figure('Name', imageName);

    % Draw a plot showing clustering parameter
    ymax = max(ripleySlopeStart) + (0.1*max (ripleySlopeStart));
   
    % Plot the running average
    if radioButtons.runningaverage
        plot (xAxisStart, ripleySlopeStart, 'b', xAxisStart, raRipleySlopeStart, 'c'); 
        legend('Slope start value','Avg slope start value',2);
    else
        plot (xAxisStart, ripleySlopeStart, 'b'); 
        %legend('Slope start value',2);
    end

    % If needed show the fitted line on the plot
    if radioButtons.plotestimate
       hold on;
       [xPlot2, yPlot2, sigma2, est2] = ptEstimateLine (xAxisStart, ripleySlopeStart, [], drugTimepoint);
       hold off;
    end
                        
    title ('Avg Slope Start Value');
    xlabel ('Frames');
    ylabel ('Avg Slope Start Value');
    if length (xAxisStart) > 1
       axis ([xAxisStart(1) xAxisStart(end) 0 ymax]);
    else
       axis ([xAxisStart(1) xAxisStart(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig2,[SaveDir filesep [imageName '_avgRipleySlopeStart.fig']]);
    % print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.eps']],'-depsc2','-tiff');
    % print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.tif']],'-dtiff'); 
    
    % generate the average plot
    h_fig3 = figure('Name', imageName);

    % Draw a plot showing clustering parameter
    ymax = max(AvgDRS) + (0.1*max (AvgDRS));
    plot (xAxisDer, AvgDRS, 'b');
              
    % If needed show the fitted line on the plot
    % not sure this works right away, so cut it out for now..
    % if radioButtons.plotestimate
    %    hold on;
    %    [xPlot2, yPlot2, sigma2, est2] = ptEstimateLine (xAxis, ripleySlopeStart, [], drugTimepoint);
    %    hold off;
    % end
                        
    title ('AvgDerOfRipStart');
    xlabel ('Frames');
    ylabel ('AvgDerOfRipStart');
    if length (xAxisDer) > 1
       axis ([xAxisDer(1) xAxisDer(end) 0 ymax]);
    else
       axis ([xAxisDer(1) xAxisDer(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format  
    hgsave (h_fig3,[SaveDir filesep [imageName '_AvgDerOfRipStart.fig']]);
    % print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.eps']],'-depsc2','-tiff');
    % print (h_fig2, [SaveDir filesep [imageName '_avgRipleySlopeStart.tif']],'-dtiff'); 
    
    
end

if ~radioButtons.donotshowplots
    % Save CSV files
    csvwrite ([SaveDir filesep imageName '_avgRipleySlopeInclination.csv'], [xAxisSlope ; ripleySlopeInclin]);
    csvwrite ([SaveDir filesep imageName '_avgRipleySlopeStart.csv'], [xAxisStart ; ripleySlopeStart]);
    csvwrite ([SaveDir filesep imageName '_AvgDerOfRipStart.csv'], [xAxisDer ; AvgDRS]);
end

if radioButtons.plotestimate
    if ~radioButtons.donotshowplots 
        csvwrite ([SaveDir filesep imageName '_fittedCurveRipleySlopeInclination.csv'], [xAxisSlope ; yPlot1]);
        csvwrite ([SaveDir filesep imageName '_curveEstimatesRipleySlopeInclination.csv'], est1);
        
        csvwrite ([SaveDir filesep imageName '_fittedCurveRipleySlopeStart.csv'], [xPlot2 ; yPlot2]);
        csvwrite ([SaveDir filesep imageName '_curveEstimatesRipleySlopeStart.csv'], est2);
    end
end
end
