function ptPlotSpeedValues (ptPostpro, MPM)
% ptPlotSpeedValues plots speed information gathered in MPM. 
%
% SYNOPSIS       ptPlotSpeedValues (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI
%                MPM       : matrix containing the cell tracks
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotSpeedValues  uses {nothing}
%                                  
%                ptPlotSpeedValues is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file
% Andre Kerstens        Jun 04          Changed ymax to 100% for percentage all cells
% Andre Kerstens        Aug 04          Fixed bug in velocity squared and made loop more 
%                                       effective (less calculations)

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagenamenotiff;
increment = ptPostpro.increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
cellProps = ptPostpro.cellProps;
clusterProps = ptPostpro.clusterProps;
binSize = ptPostpro.binsize;
frameInterval = round (ptPostpro.timeperframe / 60);    % In minutes
pixelLength = ptPostpro.mmpixel;

% In case the velocity over multiple frames is needed the following
% variable should be used
multipleFrameVelocity = ptPostpro.multFrameVelocity;

% Initialize the displacement and x-axis matrix
%displacementAll = zeros (size (MPM,1), numberOfFrames);
avgDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgSingleDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgClusteredDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);

% Initialize MPM counter
MPMCount = ceil ((plotStartFrame - startFrame) / increment) + 1;

% Initialize X-axis vector
xAxis = zeros (1, numberOfFrames-multipleFrameVelocity);
iCount = 0;

% Initialize storage
varVelocity = [];
varSingleCellVelocity = [];
varClusteredCellVelocity = [];
avgDisplacement = [];
avgSingleDisplacement = [];
avgClusteredDisplacement = [];
avgVelocity = [];
avgVelocitySquared = [];
avgSingleVelocity = [];
avgClusteredVelocity = [];
velAllCellsHigherThanAvgSingleCells = [];
displacementHist = [];
velocityHist = [];
maxVelocity = [];

% Go through every frame of the set. Start at the second frame
% because only from there we can start calculating a displacement
for frameCount = plotStartFrame+increment : increment : plotEndFrame
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % If the velocity is calculated over multiple frames we should have
   % enough frames to do this
   if MPMCount > multipleFrameVelocity
   
      % Store the frame number for display on the x-axis
      iCount = iCount + 1;
      xAxis (iCount) = frameCount;
   
      % Get the coordinates from the previous (dependent on the value of
      % multipleFrameVelocity) and the current frame
      prevCoordinates = MPM (:, 2*(MPMCount-multipleFrameVelocity)-1 : 2*(MPMCount-multipleFrameVelocity));
      curCoordinates = MPM (:, 2*MPMCount-1 : 2*MPMCount);
      coordinates = [prevCoordinates curCoordinates];
   
      % If coordinates in the current frame have started a track or the
      % ones in the previous frame ended a track, we cannot
      % calculate a displacement so we should throw these out
      coordinates (find ((coordinates (:,1) == 0 & coordinates (:,2) == 0) | ...
                         (coordinates (:,3) == 0 & coordinates (:,4) == 0)),:) = [];
      
      % From the cluster and cell properties we can find out which cells
      % are single and which belong to a cluster. First initialize two
      % matrices for this
      singleCellIndex = zeros (size (coordinates,1), 1);
      clusteredCellIndex = zeros (size (coordinates,1), 1);
      
      % Initialize single and clustered cell coordinates counters
      singleCellCount = 0;
      clusteredCellCount = 0;
      
      % Find the single cell and clustered coordinates and store these
      for jCount = 1 : size (coordinates,1) 

         % Find the index of the cell in CellProps
         cellIndex = find (cellProps (:,1,MPMCount) == coordinates (jCount,3) & ...
                           cellProps (:,2,MPMCount) == coordinates (jCount,4));
         
		 % In case the cell wasn't found in cellProps (maybe it was on the background?)
         % we can skip over this part
         if ~isempty (cellIndex)
             
            % Find this cluster in clusterProps
            clusterIndex = find (clusterProps (:,1,MPMCount) == cellProps (cellIndex,3,MPMCount));
         
            % Test whether this cell is in a cluster or not
            if clusterProps (clusterIndex,2,MPMCount) == 1  % Single cell found
            
               % Increase the counter for the single cell matrix
               singleCellCount = singleCellCount + 1; 
            
               % Store the coordinates
               singleCellIndex (singleCellCount,:) = jCount;
            else    % The cell was in a cluster
            
               % Increase the counter for the clustered cell matrix
               clusteredCellCount = clusteredCellCount + 1; 
            
               % Store the coordinates
               clusteredCellIndex (clusteredCellCount,:) = jCount;
            end
         end  % if ~isempty
      end  % for jCount = 1 : size (coordinates,1)
      
      % Remove all the zero entries in both matrices
      singleCellIndex (find (singleCellIndex(:,1) == 0), :) = [];
      clusteredCellIndex (find (clusteredCellIndex(:,1) == 0), :) = [];
      
      % Initialize displacement matrices
      displacement = zeros (size (coordinates,1),1);
      velocity = zeros (size (coordinates,1),1);
      singleCellsDisplacement = zeros (size (singleCellIndex,1),1);
      singleCellsVelocity = zeros (size (singleCellIndex,1),1);
      clusteredCellsDisplacement = zeros (size (clusteredCellIndex,1),1);
      clusteredCellsVelocity = zeros (size (clusteredCellIndex,1),1);
      
      
      % Calculate the displacement for all cells
      displacement = sqrt ((coordinates(:,1) - coordinates(:,3)).^2 + ...
                           (coordinates(:,2) - coordinates(:,4)).^2);
                      
      % Calculate true velocity in um/min for all cells
      velocity = (displacement * pixelLength) / (multipleFrameVelocity * frameInterval);
            
      % Calculate the displacement for single cells
      singleCellsDisplacement = displacement(singleCellIndex,:);
      
      % Calculate true velocity in um/min for all cells
      singleCellVelocity = velocity(singleCellIndex,:);
      
      % Calculate the displacement for clustered cells
      clusteredCellsDisplacement = displacement(clusteredCellIndex);
      
      % Calculate true velocity in um/min for all cells
      clusteredCellVelocity = velocity(clusteredCellIndex);
      
      % Calculate the variance of the velocity for all cells
      if ~isempty (velocity)
         varVelocity (iCount) = var (velocity);
      else
         varVelocity (iCount) = 0;
      end
      
      % Calculate the variance of the velocity for single cells
      if ~isempty (singleCellVelocity)
         varSingleCellVelocity (iCount) = var (singleCellVelocity);
      else
         varSingleCellVelocity (iCount) = 0;
      end
      
      % Calculate the variance of the velocity for clustered cells
      if ~isempty (clusteredCellVelocity)
         varClusteredCellVelocity (iCount) = var (clusteredCellVelocity);
      else
         varClusteredCellVelocity (iCount) = 0;
      end
                     
      % From this the average displacement for all cells can be calculated
      if length (displacement) > 0
         avgDisplacement (iCount) = sum (displacement) / length (displacement);
      else
         avgDisplacement (iCount) = 0;
      end

      % From this the average displacement for single cells can be calculated
      if length (singleCellsDisplacement) > 0
         avgSingleDisplacement (iCount) = sum (singleCellsDisplacement) / length (singleCellsDisplacement);
      else
         avgSingleDisplacement (iCount) = 0;
      end

      % From this the average displacement for clustered cells can be calculated
      if length (clusteredCellsDisplacement) > 0
         avgClusteredDisplacement (iCount) = sum (clusteredCellsDisplacement) / length (clusteredCellsDisplacement);
      else
         avgClusteredDisplacement (iCount) = 0;
      end
      
      % From this the average velocity for all cells can be calculated
      if length (displacement) > 0
         avgVelocity (iCount) = sum (velocity) / length (velocity);
         avgVelocitySquared(iCount) = avgVelocity(iCount)^2;
      else
         avgVelocity (iCount) = 0;
         avgVelocitySquared(iCount) = 0;
      end

      % From this the average velocity for single cells can be calculated
      if length (singleCellsDisplacement) > 0
         avgSingleVelocity (iCount) = sum (singleCellVelocity) / length (singleCellVelocity);
      else
         avgSingleVelocity (iCount) = 0;
      end

      % From this the average velocity for clustered cells can be calculated
      if length (clusteredCellsDisplacement) > 0
         avgClusteredVelocity (iCount) = sum (clusteredCellVelocity) / length (clusteredCellVelocity);
      else
         avgClusteredVelocity (iCount) = 0;
      end
      
      % Calculate how many cells (of all cells) have a velocity higher or
      % equal than the average single cell velocity (in percent)
      if length (velocity) ~= 0
         velAllCellsHigherThanAvgSingleCells (iCount) = (length (find (velocity >= avgSingleVelocity (iCount))) / ...
                                                        length (velocity)) * 100.0;
      else
         velAllCellsHigherThanAvgSingleCells (iCount) = 0;
      end
      
      % Generate a displacement histogram
      %displacementHist (:,iCount) = hist (displacement, binSize);
      
      % Generate a velocity histogram
      velocityHist (:,iCount) = hist (velocity, binSize);
      
      % Store max velocity of the frame
      maxVelocity(iCount) = max (velocity);
      
   end   % if MPMCount > multipleFrameVelocity 
      
end   % for frameCount

% Calculate the maximum velocity value
maxVel = max (maxVelocity);

% Generate the bin vector by using the GUI field binsize
bin = [];
for binCount = 1 : binSize
   binCol = round (binCount * maxVel / binSize);
   bin = [bin binCol];
end


% Here is where all the plotting starts

if ptPostpro.speedplot_2

    % Generate the avg velocity plot (all cells)
    h_fig = figure('Name', imageName);

    % Draw a plot showing average velocity of all cells
    ymax = max (avgVelocity) + 1;
    subplot (2,1,1); plot (xAxis, avgVelocity); 
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
    title ('Avg Squared Velocity All Cells');
    xlabel ('Frames');
    ylabel ('Velocity^2 (um/min)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format        
    hgsave (h_fig,[savePath filesep [imageName '_avgVelocityAllCells.fig']]);
    print (h_fig, [savePath filesep [imageName '_avgVelocityAllCells.eps']],'-depsc2','-tiff');
    print (h_fig, [savePath filesep [imageName '_avgVelocityAllCells.tif']],'-dtiff');      

    % Generate the figure and title
    h_fig2 = figure('Name', imageName);

    % Draw a subplot showing the avg velocity of a single cell    
    ymax = max (avgSingleVelocity) + 1;
    subplot (2,1,1); plot (xAxis, avgSingleVelocity);
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
    title ('Avg Clustered Cell Velocity');
    xlabel ('Frames');
    ylabel ('Velocity (um/min)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save MAT files for avg all, single and clustered cell velocity
    cd (savePath);
    save ([imageName '_avgCellVelocity.mat'],'avgVelocity');
    save ([imageName '_avgSingleCellVelocity.mat'],'avgSingleVelocity');
    save ([imageName '_avgClusteredCellVelocity.mat'],'avgClusteredVelocity');

    % Save CSV files for avg all, single and clustered cell velocity
    csvwrite ([imageName '_avgCellVelocity.csv'], [xAxis ; avgVelocity]);
    csvwrite ([imageName '_avgSingleCellVelocity.csv'], [xAxis ; avgSingleVelocity]);
    csvwrite ([imageName '_avgClusteredCellVelocity.csv'], [xAxis ; avgClusteredVelocity]);

    % Save the figures in fig, eps and tif format     
    hgsave (h_fig2,[savePath filesep [imageName '_avgSingleAndClusterVelocity.fig']]);
    print (h_fig2, [savePath filesep [imageName '_avgSingleAndClusterVelocity.eps']],'-depsc2','-tiff');
    print (h_fig2, [savePath filesep [imageName '_avgSingleAndClusterVelocity.tif']],'-dtiff');
end


if ptPostpro.speedplot_1

    % Generate the plot that shows all cells velocity against single cell
    % velocity (in percent)
    h_fig3 = figure('Name', imageName);

    % Draw the plot
    ymax = 100;
    plot (xAxis, velAllCellsHigherThanAvgSingleCells); 
    title ('% All Cells with Velocity higher than Average Single Cell Velocity');
    xlabel ('Frames');
    ylabel ('Percentage (%)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save MAT files for percentage velocity higher than single cell speed
    cd (savePath);
    save ([imageName '_velAllCellsHigherThanAvgSingleCells.mat'],'velAllCellsHigherThanAvgSingleCells');

    % Save CSV files for percentage velocity higher than single cell speed
    csvwrite ([imageName '_velAllCellsHigherThanAvgSingleCells.csv'], [xAxis ; velAllCellsHigherThanAvgSingleCells]);

    % Save the figures in fig, eps and tif format        
    hgsave (h_fig3,[savePath filesep [imageName '_velAllCellsHigherThanAvgSingleCells.fig']]);
    print (h_fig3, [savePath filesep [imageName '_velAllCellsHigherThanAvgSingleCells.eps']],'-depsc2','-tiff');
    print (h_fig3, [savePath filesep [imageName '_velAllCellsHigherThanAvgSingleCells.tif']],'-dtiff');
end


if ptPostpro.speedplot_3

    % Generate the figure and title for the variance plot
    h_fig4 = figure('Name', imageName);

    % Draw a subplot showing the velocity variance of all cells    
    ymax = max (varVelocity)+1;
    subplot (3,1,1); plot (xAxis, varVelocity);
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
    title ('Variance Clustered Cell Velocity');
    xlabel ('Frames');
    ylabel ('Variance');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save MAT files for avg all, single and clustered cell variance
    cd (savePath);
    save ([imageName '_varCellVelocity.mat'],'varVelocity');
    save ([imageName '_varSingleCellVelocity.mat'],'varSingleCellVelocity');
    save ([imageName '_varClusteredCellVelocity.mat'],'varClusteredCellVelocity');

    % Save CSV files for avg all, single and clustered cell variance
    csvwrite ([imageName '_varCellVelocity.csv'], [xAxis ; varVelocity]);
    csvwrite ([imageName '_varSingleCellVelocity.csv'], [xAxis ; varSingleCellVelocity]);
    csvwrite ([imageName '_varClusteredCellVelocity.csv'], [xAxis ; varClusteredCellVelocity]);

    % Save the figures in fig, eps and tif format     
    hgsave (h_fig4,[savePath filesep [imageName '_varSingleAndClusterVelocity.fig']]);
    print (h_fig4, [savePath filesep [imageName '_varSingleAndClusterVelocity.eps']],'-depsc2','-tiff');
    print (h_fig4, [savePath filesep [imageName '_varSingleAndClusterVelocity.tif']],'-dtiff');
end


% %Generate a 3D velocity histogram using the surface function
% h_fig = figure('Name', imageName);
% surf (displacementHist);
% title ('3D Velocity Histogram');
% 
% % Save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep [imageName '_3dDisplacementHistogram.fig']]);
% print (h_fig, [savePath filesep [imageName '_3dDisplacementHistogram.eps']], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep [imageName '_3dDisplacementHistogram.tif']], '-dtiff');


% Generate a 3D velocity histogram using 3-d bars which is chopped up in bins with size binSize

% Generate a 3D displacement histogram using 3-d bars which is chopped up in bins
h_fig = figure;
bar3 (bin, velocityHist, 0.5, 'detached');
title ('3D Velocity Histogram (binned)');

% Save this figure to disk as fig, eps and tiff
hgsave (h_fig, [savePath filesep [imageName '_3dBinnedVelocityHistogram.fig']]);
print (h_fig, [savePath filesep [imageName '_3dBinnedVelocityHistogram.eps']], '-depsc2', '-tiff');
print (h_fig, [savePath filesep [imageName '_3dBinnedVelocityHistogram.tif']], '-dtiff');


% For all the figures we want to keep the xAxis as well 
cd (savePath);
save ([imageName '_xAxis-Velocity.mat'],'xAxis');
save ([imageName '_velocityHist.mat'], 'velocityHist');

