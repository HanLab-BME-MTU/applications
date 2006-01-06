function ptPlotVelocityHist (radioButtons, imageName, SaveDir, xAxis, velocityHistStats)
% ptPlotSpeedVarStats generates the velocity histogram plots
%
% SYNOPSIS       ptPlotHist (radioButtons, imageName, SaveDir, xAxis, velocityHistStats)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name for the title of the plots
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                velocityHistStats : struct with the following fields:
%
%                   velocityHist : vector with velocity histogram values
%                   histVelAllCells : cell containing velocity histograms per frame
%                   histVelSingleCells : cell containing single cell velocity histograms per frame
%                   histVelClusteredCells : cell containing clustered cell velocity histograms per frame
%                   maxVelocity : vector with max velocity values
%                   maxSingleCellVelocity : vector with max velocity values for single cells
%                   maxClusteredCellVelocity : vector with max velocity values for clustered cells
%                   binSize : value for bin size
%                
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotVelocityHist  uses { nothing }
%                                  
%                ptPlotVelocityHist is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

% Get data from the input struct
velocityHist = velocityHistStats.velocityHist;
binSize = velocityHistStats.binSize;

if radioButtons.speedplot_4

    if ~radioButtons.donotshowplots

        % Generate the bin vector by using the GUI field binsize
        bin = [1:binSize];

        % Generate a 3D displacement histogram using 3-d bars which is chopped up in bins
        h_fig = figure;
        bar3 (bin, velocityHist, 0.5, 'detached');
        title ([imageName '3D Velocity Histogram (binned)']);

        % Save this figure to disk as fig, eps and tiff
        hgsave (h_fig, [SaveDir filesep [imageName '_3dBinnedVelocityHistogram.fig']]);
        print (h_fig, [SaveDir filesep [imageName '_3dBinnedVelocityHistogram.eps']], '-depsc2', '-tiff');
        print (h_fig, [SaveDir filesep [imageName '_3dBinnedVelocityHistogram.tif']], '-dtiff');
    end  % if ~radioButtons.donotshowplots
end

% Save histogram for velocity all cells
if radioButtons.allcellshist
    for iCount = 1 : length(velocityHistStats.histVelAllCells)
       filename = [SaveDir filesep 'histVelAllCells_' num2str(iCount) '.mat'];
       histVect = velocityHistStats.histVelAllCells{iCount};
       save (filename, 'histVect');
    end
    maxVelocity = velocityHistStats.maxVelocity;
    save ([SaveDir filesep 'maxVelAllCells.mat'], 'maxVelocity');
end
            
% Save the velocity hist for single cells
if radioButtons.singlecellshist
    for iCount = 1 : length(velocityHistStats.histVelSingleCells)
       filename = [SaveDir filesep 'histVelSingleCells_' num2str(iCount) '.mat'];
       histVect = velocityHistStats.histVelSingleCells{iCount};
       save (filename, 'histVect');
    end
    maxSingleCellVelocity = velocityHistStats.maxSingleCellVelocity;
    save ([SaveDir filesep 'maxVelSingleCells.mat'], 'maxSingleCellVelocity');
end
     
% Save the velocity hist for clustered cells
if radioButtons.clusteredcellshist
    for iCount = 1 : length(velocityHistStats.histVelClusteredCells)
       filename = [SaveDir filesep 'histVelClusteredCells_' num2str(iCount) '.mat'];
       histVect = velocityHistStats.histVelClusteredCells{iCount};
       save (filename, 'histVect');
    end
    maxClusteredCellVelocity = velocityHistStats.maxClusteredCellVelocity;
    save ([SaveDir filesep 'maxVelClusteredCells.mat'], 'maxClusteredCellVelocity');
end

% For all the figures we want to keep the xAxis as well 
save ([SaveDir filesep 'velocityHist.mat'], 'velocityHist');
