function ptPlotSpeedVarStats (radioButtons, imageName, savePath, xAxis, velocityVarStats, windowSize)
% ptPlotSpeedVarStats generates the plots for the velocity statistics
%
% SYNOPSIS       ptPlotSpeedVarStats (radioButtons, imageName, savePath, xAxis, velocityVarStats)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name of the first image in the movie (used as title)
%                savePath : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                velocityVarStats : struct with the following fields:
%                   varVelocity : vector with velocity variance all cells
%                   varSingleCellVelocity : vector with velocity variance single cells
%                   varClusteredCellVelocity : vector with velocity variance clustered cells
%                windowSize : size of the averaging window  
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotsSpeedVarStats  uses { nothing }
%                                  
%                ptPlotSpeedVarStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

% Get data from struct
varVelocity = velocityVarStats.varVelocity;
varSingleCellVelocity = velocityVarStats.varSingleCellVelocity;
varClusteredCellVelocity = velocityVarStats.varClusteredCellVelocity;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raVarVelocity = movingAverage (varVelocity, windowSize, 'median');
    raSingleCellVelocity = movingAverage (varSingleCellVelocity, windowSize, 'median');
    raClusteredCellVelocity = movingAverage (varClusteredCellVelocity, windowSize, 'median');
end

% Here's where the plotting starts
if radioButtons.speedplot_3

    if ~radioButtons.donotshowplots

        % Generate the figure and title for the variance plot
        h_fig4 = figure('Name', imageName);

        % Draw a subplot showing the velocity variance of all cells    
        ymax = max (varVelocity)+1;
        subplot (3,1,1); plot (xAxis, varVelocity);
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raVarVelocity, 'r'); hold off;
        end
        
        title ('Variance Velocity All Cells');
        xlabel ('Frames');
        ylabel ('Variance');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end   

        % Draw a subplot showing the velocity variance of single cells
        ymax = max (varSingleCellVelocity)+1;
        subplot (3,1,2); plot (xAxis, varSingleCellVelocity); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raSingleCellVelocity, 'r'); hold off;
        end
        
        title ('Variance Single Cell Velocity');
        xlabel ('Frames');
        ylabel ('Variance');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Draw a subplot showing the velocity variance of clustered cells
        ymax = max (varClusteredCellVelocity)+1;
        subplot (3,1,3); plot (xAxis, varClusteredCellVelocity); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raClusteredCellVelocity, 'r'); hold off;
        end
        
        title ('Variance Clustered Cell Velocity');
        xlabel ('Frames');
        ylabel ('Variance');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Save the figures in fig, eps and tif format     
        hgsave (h_fig4,[savePath filesep [imageName '_varSingleAndClusterVelocity.fig']]);
        print (h_fig4, [savePath filesep [imageName '_varSingleAndClusterVelocity.eps']],'-depsc2','-tiff');
        print (h_fig4, [savePath filesep [imageName '_varSingleAndClusterVelocity.tif']],'-dtiff');
    end   % if ~radioButtons.donotshowplots
    
    % Save MAT files for avg all, single and clustered cell variance
    save ([savePath filesep imageName '_varCellVelocity.mat'],'varVelocity');
    save ([savePath filesep imageName '_varSingleCellVelocity.mat'],'varSingleCellVelocity');
    save ([savePath filesep imageName '_varClusteredCellVelocity.mat'],'varClusteredCellVelocity');

    % Save CSV files for avg all, single and clustered cell variance
    csvwrite ([savePath filesep imageName '_varCellVelocity.csv'], [xAxis ; varVelocity]);
    csvwrite ([savePath filesep imageName '_varSingleCellVelocity.csv'], [xAxis ; varSingleCellVelocity]);
    csvwrite ([savePath filesep imageName '_varClusteredCellVelocity.csv'], [xAxis ; varClusteredCellVelocity]);
end
