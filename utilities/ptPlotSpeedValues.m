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

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
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
      singleCellCoords = zeros (size (coordinates,1), 4);
      clusteredCellCoords = zeros (size (coordinates,1), 4);
      
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
               singleCellCoords (singleCellCount,:) = coordinates (jCount,:);
            else    % The cell was in a cluster
            
               % Increase the counter for the clustered cell matrix
               clusteredCellCount = clusteredCellCount + 1; 
            
               % Store the coordinates
               clusteredCellCoords (clusteredCellCount,:) = coordinates (jCount,:);
            end
         end  % if ~isempty
      end  % for jCount = 1 : size (coordinates,1)
      
      % Remove all the zero entries in both matrices
      singleCellCoords (find (singleCellCoords(:,1) == 0), :) = [];
      clusteredCellCoords (find (clusteredCellCoords(:,1) == 0), :) = [];
      
      % Initialize displacement matrices
      displacement = zeros (size (coordinates,1),1);
      singleCellsDisplacement = zeros (size (singleCellCoords,1),1);
      clusteredCellsDisplacement = zeros (size (clusteredCellCoords,1),1);
      
      % Calculate the displacement for all cells
      displacement = sqrt ((coordinates(:,1) - coordinates(:,3)).^2 + ...
                           (coordinates(:,2) - coordinates(:,4)).^2);
                      
      % Calculate true velocity in um/min for all cells
      velocity = (displacement * pixelLength) / (multipleFrameVelocity * frameInterval);
      
      % Calculate the displacement for single cells
      singleCellsDisplacement = sqrt ((singleCellCoords(:,1) - singleCellCoords(:,3)).^2 + ...
                                      (singleCellCoords(:,2) - singleCellCoords(:,4)).^2);
      
      % Calculate true velocity in um/min for all cells
      singleCellVelocity = (singleCellsDisplacement * pixelLength) / (multipleFrameVelocity * frameInterval);                
      
      % Calculate the displacement for clustered cells
      clusteredCellsDisplacement = sqrt ((clusteredCellCoords(:,1) - clusteredCellCoords(:,3)).^2 + ...
                                         (clusteredCellCoords(:,2) - clusteredCellCoords(:,4)).^2);   
      
      % Calculate true velocity in um/min for all cells
      clusteredCellVelocity = (clusteredCellsDisplacement * pixelLength) / (multipleFrameVelocity * frameInterval);                               
      
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
         avgSingleVelocity (iCount) = sum (singleCellsDisplacement) / length (singleCellsDisplacement);
      else
         avgSingleVelocity (iCount) = 0;
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
      else
         avgVelocity (iCount) = 0;
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
      displacementHist (:,iCount) = hist (displacement, binSize);
      
      % Generate a velocity histogram
      velocityHist (:,iCount) = hist (velocity, binSize);
      
   end   % if MPMCount > multipleFrameVelocity 
      
end   % for frameCount


% % Generate the avg displacement plot (all cells)
% h_fig = figure; title (imageName);
%    
% % Draw a plot showing average displacement of all cells
% ymax = max (avgDisplacement)+1;
% plot (xAxis, avgDisplacement); 
% title ('Avg Displacement All Cells');
% xlabel ('Frames');
% ylabel ('Displacement');
% if length (xAxis) > 1
%    axis ([xAxis(1) xAxis(end) 0 ymax]);
% else
%    axis ([xAxis(1) xAxis(1)+1 0 ymax]);
% end
% 
% % Save the figures in fig, eps and tif format        
% hgsave (h_fig,[savePath filesep 'avgDisplacementAllCells.fig']);
% print (h_fig, [savePath filesep 'avgDisplacementAllCells.eps'],'-depsc2','-tiff');
% print (h_fig, [savePath filesep 'avgDisplacementAllCells.tif'],'-dtiff');      



% % Generate the figure and title
% h_fig = figure; title (imageName);
%    
% % Draw a subplot showing the avg displacement of a single cell    
% ymax = max (avgSingleDisplacement)+1;
% subplot (2,1,1); plot (xAxis, avgSingleDisplacement);
% title ('Avg Single Cell Displacement');
% xlabel ('Frames');
% ylabel ('Displacement');
% if length (xAxis) > 1
%    axis ([xAxis(1) xAxis(end) 0 ymax]);
% else
%    axis ([xAxis(1) xAxis(1)+1 0 ymax]);
% end   
% 
% % Draw a subplot showing the avg displacement of a cluster
% ymax = max (avgClusteredDisplacement)+1;
% subplot (2,1,2); plot (xAxis, avgClusteredDisplacement); 
% title ('Avg Clustered Cell Displacement');
% xlabel ('Frames');
% ylabel ('Displacement');
% if length (xAxis) > 1
%    axis ([xAxis(1) xAxis(end) 0 ymax]);
% else
%    axis ([xAxis(1) xAxis(1)+1 0 ymax]);
% end
% 
% % Save MAT files for avg all, single and clustered cell displacement
% cd (savePath);
% save ('avgCellDisplacement.mat','avgDisplacement');
% save ('avgSingleCellDisplacement.mat','avgSingleDisplacement');
% save ('avgClusteredCellDisplacement.mat','avgClusteredDisplacement');
%    
% % Save CSV files for avg all, single and clustered cell displacement
% csvwrite ('avgCellDisplacement.csv', [xAxis ; avgDisplacement]);
% csvwrite ('avgSingleCellDisplacement.csv', [xAxis ; avgSingleDisplacement]);
% csvwrite ('avgClusteredCellDisplacement.csv', [xAxis ; avgClusteredDisplacement]);
%    
% % Save the figures in fig, eps and tif format     
% hgsave (h_fig,[savePath filesep 'avgSingleAndClusterDisplacement.fig']);
% print (h_fig, [savePath filesep 'avgSingleAndClusterDisplacement.eps'],'-depsc2','-tiff');
% print (h_fig, [savePath filesep 'avgSingleAndClusterDisplacement.tif'],'-dtiff');


% Generate the avg velocity plot (all cells)
h_fig = figure('Name', imageName);
   
% Draw a plot showing average velocity of all cells
ymax = max (avgVelocity) + 1;
plot (xAxis, avgVelocity); 
title ('Avg Velocity All Cells');
xlabel ('Frames');
ylabel ('Velocity (um/min)');
if length (xAxis) > 1
   axis ([xAxis(1) xAxis(end) 0 ymax]);
else
   axis ([xAxis(1) xAxis(1)+1 0 ymax]);
end

% Save the figures in fig, eps and tif format        
hgsave (h_fig,[savePath filesep 'avgVelocityAllCells.fig']);
print (h_fig, [savePath filesep 'avgVelocityAllCells.eps'],'-depsc2','-tiff');
print (h_fig, [savePath filesep 'avgVelocityAllCells.tif'],'-dtiff');      


% Generate the figure and title
h_fig = figure('Name', imageName);
   
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
save ('avgCellVelocity.mat','avgVelocity');
save ('avgSingleCellVelocity.mat','avgSingleVelocity');
save ('avgClusteredCellVelocity.mat','avgClusteredVelocity');
   
% Save CSV files for avg all, single and clustered cell velocity
csvwrite ('avgCellVelocity.csv', [xAxis ; avgVelocity]);
csvwrite ('avgSingleCellVelocity.csv', [xAxis ; avgSingleVelocity]);
csvwrite ('avgClusteredCellVelocity.csv', [xAxis ; avgClusteredVelocity]);
   
% Save the figures in fig, eps and tif format     
hgsave (h_fig,[savePath filesep 'avgSingleAndClusterVelocity.fig']);
print (h_fig, [savePath filesep 'avgSingleAndClusterVelocity.eps'],'-depsc2','-tiff');
print (h_fig, [savePath filesep 'avgSingleAndClusterVelocity.tif'],'-dtiff');


% Generate the plot that shows all cells velocity against single cell
% velocity (in percent)
h_fig = figure('Name', imageName);
   
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
save ('velAllCellsHigherThanAvgSingleCells.mat','velAllCellsHigherThanAvgSingleCells');
   
% Save CSV files for percentage velocity higher than single cell speed
csvwrite ('velAllCellsHigherThanAvgSingleCells.csv', [xAxis ; velAllCellsHigherThanAvgSingleCells]);

% Save the figures in fig, eps and tif format        
hgsave (h_fig,[savePath filesep 'velAllCellsHigherThanAvgSingleCells.fig']);
print (h_fig, [savePath filesep 'velAllCellsHigherThanAvgSingleCells.eps'],'-depsc2','-tiff');
print (h_fig, [savePath filesep 'velAllCellsHigherThanAvgSingleCells.tif'],'-dtiff');


% Generate the figure and title for the variance plot
h_fig = figure('Name', imageName);
   
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
save ('varCellVelocity.mat','varVelocity');
save ('varSingleCellVelocity.mat','varSingleCellVelocity');
save ('varClusteredCellVelocity.mat','varClusteredCellVelocity');
   
% Save CSV files for avg all, single and clustered cell variance
csvwrite ('varCellVelocity.csv', [xAxis ; varVelocity]);
csvwrite ('varSingleCellVelocity.csv', [xAxis ; varSingleCellVelocity]);
csvwrite ('varClusteredCellVelocity.csv', [xAxis ; varClusteredCellVelocity]);
   
% Save the figures in fig, eps and tif format     
hgsave (h_fig,[savePath filesep 'varSingleAndClusterVelocity.fig']);
print (h_fig, [savePath filesep 'varSingleAndClusterVelocity.eps'],'-depsc2','-tiff');
print (h_fig, [savePath filesep 'varSingleAndClusterVelocity.tif'],'-dtiff');


% %Generate a 3D velocity histogram using the surface function
% h_fig = figure('Name', imageName);
% surf (displacementHist);
% title ('3D Velocity Histogram');
% 
% % Save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep '3dDisplacementHistogram.fig']);
% print (h_fig, [savePath filesep '3dDisplacementHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep '3dDisplacementHistogram.tif'], '-dtiff');


% Generate a 3D velocity histogram using 3-d bars which is chopped up in bins with size binSize

% % Calculate the maximum displacement value
% maxDisp = max (max (displacementHist));
% 
% % Generate the bin vector by using the GUI field binsize
% bin = [];
% for binCount = 1 : binSize
%    binCol = round (binCount * maxDisp / binSize);
%    bin = [bin binCol];
% end
% 
% % Generate a 3D displacement histogram using 3-d bars which is chopped up in bins
% h_fig = figure;
% %bar3 (bin, displacementHist, 0.5, 'detached');
% bar3 (bin, displacementHist, 0.5, 'detached');
% title ('3D Displacement Histogram (binned)');
% 
% % Save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep '3dBinnedDisplacementHistogram.fig']);
% print (h_fig, [savePath filesep '3dBinnedDisplacementHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep '3dBinnedDisplacementHistogram.tif'], '-dtiff');


% For all the figures we want to keep the xAxis as well 
save ('xAxis-CellDisp.mat','xAxis');




% % Prepare displacement matrix for clusters and single cells
% clusterDisplacement = zeros(size(handles.clustercells));
% singleDisplacement = zeros(size(handles.singlecells));
% 
% % Prepare vector that holds single and cluster cell numbers
% thresholdedCells = zeros(numberOfFrames,4);
% curCoordinatesSingleCell = zeros(numberOfFrames,2);
% curCoordinatesClusters = zeros(numberOfFrames,2);
% curCoordinates = zeros(numberOfFrames,4);
% 
% % Setup the binning used to define granularity of 3D velocity histogram
% % This number is based on the max track distance
% binSize = [0 : 1 : (ptPostpro.maxsearch + 10)];
% 
% % Now we go through every frame of the set
% for frameCount = startFrame : (endFrame - 1)
% 
%    % See if the user selected a specific range of cells and act accordingly
%    if ~isempty(handles.selectedcells) 
%       whatCells = handles.selectedcells;
%    else  
%       % If nothing has been selected use all the cells in this frame
%       whatCells = [1:length(handles.MPM(:,(2*frameCount-1)))]' ;
%    end
% 	    
%    % Take four columns of MPM [x1,y1,x2,y2]. Basically every row is two
%    % sets of coordinates, corresponding to two subsequent frames.
%    % For every frame take these 2 x and 2 y coordinates out of MPM and assign
%    % these to curCoordinates. From these the displacement can be calculated 
%    % from frame to frame later on
%    %curCoordinates = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount + 2));
% 
%    % The following is an experiment for the poster to the take the displacement
%    % over 3 frames instead of 1 (x1 to x4)
%    % AK: This should be made a variable in the gui
%    % Since to calculate the average velocity we need the frame to frame
%    % displacement, we keep that as well
%    curCoordinates1 = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount));
%    curCoordinates2 = handles.MPM (whatCells, (2 * frameCount + 5):(2 * frameCount + 6));
%    curCoordinates = [curCoordinates1 curCoordinates2];
%    curCoordinatesVel = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount + 2));
% 
%    % Find all the zeros in the curCoordinates matrix: these can then be filtered out
%    % Note: We do not erase the rows with zeros in them. In this way,
%    % the row index of displacement corresponds to the row index of MPM
%    [zeroRows,cols] = find(curCoordinates == 0);
%    [zeroRowsVel,colsVel] = find(curCoordinatesVel == 0);
%    
%    % Since we find row indexes for every colum (4 times the same) this has 
%    % to be filtered a bit
%    zeroRows = unique (zeroRows);
%    zeroRowsVel = unique (zeroRowsVel);
% 
%    % Now make sure all the rows with zeros have their colums zero'd out. We don't 
%    % delete these since we want to keep curCoordinates the same size as MPM
%    curCoordinates(zeroRows,:) = 0;
%    curCoordinatesVel(zeroRowsVel,:) = 0;
% 
%    % Calculate how many cells we have in this frame
%    numberOfCells (frameCount,1) = (length(whatCells) - length(zeroRows));
% 
%    % Calculate the displacement the cells travelled from this frame to the next
%    displacement (:,frameCount) = sqrt((curCoordinates(:,1)-curCoordinates(:,3)).^2 + (curCoordinates(:,2) - curCoordinates(:,4)).^2);
%    displacementVel (:,frameCount) = sqrt((curCoordinatesVel(:,1)-curCoordinatesVel(:,3)).^2 + (curCoordinatesVel(:,2) - curCoordinatesVel(:,4)).^2);
%    
%    % Find the indexes of all rows in clustercells which are non-zero and
%    % calculate the displacement of all the cells in those clusters
%    realClusterRows = find (handles.clustercells(:,frameCount));
%    clusterDisplacement (1:length (realClusterRows), frameCount) = displacement (handles.clustercells (realClusterRows, frameCount), frameCount);
% 
%    % Find the indexes of all rows in singlecells which are non-zero and
%    % calculate the displacement of all those cells who are NOT in a cluster
%    realSingleRows = find (handles.singlecells (:, frameCount));
%    singleDisplacement (1:length (realSingleRows), frameCount) = displacement (handles.singlecells (realSingleRows, frameCount), frameCount);
%    
%    % Calculate the average displacement for all selected cells (since we sum up, the zeros do not bother us) 
%    if (length(whatCells) - length(zeroRows)) ~= 0
%       averageDisplacementVel (frameCount,1) = sum (displacementVel(:,frameCount)) / (length(whatCells) - length(zeroRowsVel));
%       averageDisplacement (frameCount,1) = sum (displacement(:,frameCount)) / (length(whatCells) - length(zeroRows));
%    else
%       averageDisplacementVel (frameCount,1) = 0;
%        averageDisplacement (frameCount,1) = 0;
%    end
%    
%    % Calculate the average displacement for cells in a cluster
%    if (size (find (handles.clustercells (:,frameCount)), 1)) ~= 0
%       averageClusterDisplacement (frameCount, 1) = sum (clusterDisplacement(:, frameCount)) / ...
%                                                (size (find (handles.clustercells (:,frameCount)), 1));
%    else
%       averageClusterDisplacement (frameCount, 1) = 0;
%    end
% 
%    % Calculate the average displacement for cells NOT in a cluster
%    if (size (find (handles.singlecells (:,frameCount)), 1)) ~= 0
%       averageSingleDisplacement (frameCount, 1) = sum (singleDisplacement(:, frameCount)) / ...
%                                               (size (find (handles.singlecells (:,frameCount)), 1));
%    else
%       averageSingleDisplacement (frameCount, 1) = 0;
%    end
%    
%    % Calculate the variance of the displacements. Here we do have to worry about zeros
%    % and that's why we remove them first using the find function (finds all non-zero
%    % elements in displacement
%    realDisplacement = displacement (find (displacement (:,frameCount)), frameCount);
%    velocityVariance (frameCount) = var (realDisplacement);
%    
%    % Now we'll calculate single and cluster cell numbers via another method
%    thresholdValue = averageDisplacement (frameCount,1);
%    thresholdedCells (frameCount,1) = sum (realDisplacement > thresholdValue);  % Store nr of single cells
%    thresholdedCells (frameCount,2) = sum (realDisplacement <= thresholdValue); % Store nr of cluster cells
%    thresholdedCells (frameCount,3) = thresholdedCells (frameCount,1) / length (realDisplacement); % perc. single cells
%    thresholdedCells (frameCount,4) = thresholdedCells (frameCount,2) / length (realDisplacement); % perc. cluster cells
%    
%    % Generate the 2D histograms and normalize the result. This is preparation for the 3D 
%    % one which will be dealt with later
%    if (length(whatCells) - length(zeroRows)) ~= 0
%       velocityHist (:, frameCount) = hist (realDisplacement, binSize)';
%       velocityHist (:, frameCount) = velocityHist (:,frameCount) ./ (length(whatCells) - length(zeroRows));
%    else
%       velocityHist (:, frameCount) = 0;
%    end
% end

% Now we start the actual plotting of the velocity related stuff
% Every plot is shown on the screen and saved to disk in three formats: fig, tiff and eps

% Generate a 3D velocity histogram using the surface function
% h_fig = figure;
% surf (velocityHist);
% title ('3D velocity histogram');
% 
% % And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep '3dVelocityHistogram.fig']);
% print (h_fig, [savePath filesep '3dVelocityHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep '3dVelocityHistogram.tif'], '-dtiff');
% 
% % Generate a 3D velocity histogram using 3-d bars which is chopped up in bins with size binSize
% h_fig = figure;
% bar3 (binSize, velocityHist);
% title ('3D velocity histogram using bins');
% 
% % And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep '3dBinnedVelocityHistogram.fig']);
% print (h_fig, [savePath filesep '3dBinnedVelocityHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep '3dBinnedVelocityHistogram.tif'], '-dtiff');
% 
% % Generate a plot with the number of cells found per frame
% h_fig = figure, plot (numberOfCells), title ('Number of cells');
% xlabel ('Frames');
% ylabel ('# Cells');
% ymax = max (numberOfCells);
% axis ([startFrame endFrame 0 ymax]);
% 
% % And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep 'numberOfCells.fig']);
% print (h_fig, [savePath filesep 'numberOfCells.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep 'numberOfCells.tif'], '-dtiff');
% 
% Generate a plot with both average velocity and variance of single cells and 
% cells within clusters in one figure
% h_fig = figure;
% subplot (2,1,1); plot (averageClusterDisplacement), title ('Average velocity of cells within clusters');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max (averageClusterDisplacement);
% axis ([startFrame endFrame 0 ymax]);
% subplot (2,1,2); plot (averageSingleDisplacement), title('Average velocity of single cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max(averageSingleDisplacement);
% axis ([startFrame endFrame 0 ymax]);

% And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep 'velocitySingleCluster.fig']);
% print (h_fig, [savePath filesep 'velocitySingleCluster.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep 'velocitySingleCluster.tif'], '-dtiff');

% Generate a plot with both average velocity and variance of cells in one figure
% h_fig = figure;
% subplot (2,1,1); plot (averageDisplacement), title('Average velocity of cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max(averageDisplacement);
% axis ([startFrame endFrame 0 ymax]);
% subplot (2,1,2); plot (velocityVariance), title('Variance of avarage velocity of cells');
% xlabel ('Frames');
% ylabel ('variance');
% ymax = max(velocityVariance);
% axis ([startFrame endFrame 0 ymax]);

% And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [savePath filesep 'velocityVariance.fig']);
% print (h_fig, [savePath filesep 'velocityVariance.eps'], '-depsc2', '-tiff');
% print (h_fig, [savePath filesep 'velocityVariance.tif'], '-dtiff');

% Andre test
% velocityAll = realDisplacement';
% velocityAll = displacement (find (displacement(:,1)),:);
% transposedDisplacement = displacement';
% save ('transposedDisplacement','transposedDisplacement');
% h_fig = figure;
% plot (transposedDisplacement,'rx'), title('Velocity of cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');

% Lets build up a histogram as well
%dispHistogram = hist (realDisplacement (find (transposedDisplacement)), binSize);
% for jFrame = startFrame : (endFrame - 1)
%    dispHistogram (:,jFrame) = hist (displacement (find (displacement(:,jFrame))), binSize);
% end   
% h_fig = figure;
% plot (dispHistogram, 'b'), title('Histogram of velocity of cell');
% ylabel ('Frames');
% xlabel ('# Cells');
% save ('dispHistogram','dispHistogram');

%xmax = max(displacement);
%axis ([0 xmax startFrame endFrame]);
% h_fig = figure;
% frame=1
% for frame = 1 : 204
%     for cell = 1 : 666
%         plot (displacement(cell,frame),frame, 'x');
%         hold on
%     end
% end
% hold off
% xlabel ('displacement per frame (in pixel)');
% ylabel ('Frames');
% title('Velocity of cells');

% Save all of the calculated values to disk
% cd (savePath);
% save ('velocityVariance', 'velocityVariance');
% save ('displacement', 'displacement');
% save ('displacementVel', 'displacementVel');
% save ('singleDisplacement', 'singleDisplacement');
% save ('clusterDisplacement', 'clusterDisplacement');
% save ('averageDisplacement', 'averageDisplacement');
% save ('averageDisplacementVel', 'averageDisplacementVel');
% save ('averageSingleDisplacement', 'averageSingleDisplacement');
% save ('averageClusterDisplacement', 'averageClusterDisplacement');
% save ('thresholdedCells','thresholdedCells');
