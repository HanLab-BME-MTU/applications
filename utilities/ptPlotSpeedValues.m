function ptPlotSpeedValues (ptPostpro, MPM)
% ptPlotSpeedValues plots speed information gathered in MPM. 
%
% SYNOPSIS       ptPlotSpeedValues (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
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

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.plotfirstimg;
endFrame = ptPostpro.plotlastimg;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
numberOfFrames = endFrame - startFrame + 1;
 
% Initialize the displacement matrix
displacement = zeros(size(handles.MPM,1), numberOfFrames);

% Prepare displacement matrix for clusters and single cells
clusterDisplacement = zeros(size(handles.clustercells));
singleDisplacement = zeros(size(handles.singlecells));

% Prepare vector that holds single and cluster cell numbers
thresholdedCells = zeros(numberOfFrames,4);
curCoordinates1 = zeros(numberOfFrames,2);
curCoordinates2 = zeros(numberOfFrames,2);
curCoordinates = zeros(numberOfFrames,4);

% Setup the binning used to define granularity of 3D velocity histogram
% This number is based on the max track distance
binSize = [0 : 1 : (ptPostpro.maxsearch + 10)];

% Now we go through every frame of the set
for frameCount = startFrame : (endFrame - 1)

   % See if the user selected a specific range of cells and act accordingly
   if ~isempty(handles.selectedcells) 
      whatCells = handles.selectedcells;
   else  
      % If nothing has been selected use all the cells in this frame
      whatCells = [1:length(handles.MPM(:,(2*frameCount-1)))]' ;
   end
	    
   % Take four columns of MPM [x1,y1,x2,y2]. Basically every row is two
   % sets of coordinates, corresponding to two subsequent frames.
   % For every frame take these 2 x and 2 y coordinates out of MPM and assign
   % these to curCoordinates. From these the displacement can be calculated 
   % from frame to frame later on
   %curCoordinates = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount + 2));

   % The following is an experiment for the poster to the take the displacement
   % over 3 frames instead of 1 (x1 to x4)
   % AK: This should be made a variable in the gui
   % Since to calculate the average velocity we need the frame to frame
   % displacement, we keep that as well
   curCoordinates1 = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount));
   curCoordinates2 = handles.MPM (whatCells, (2 * frameCount + 5):(2 * frameCount + 6));
   curCoordinates = [curCoordinates1 curCoordinates2];
   curCoordinatesVel = handles.MPM (whatCells, (2 * frameCount - 1):(2 * frameCount + 2));

   % Find all the zeros in the curCoordinates matrix: these can then be filtered out
   % Note: We do not erase the rows with zeros in them. In this way,
   % the row index of displacement corresponds to the row index of MPM
   [zeroRows,cols] = find(curCoordinates == 0);
   [zeroRowsVel,colsVel] = find(curCoordinatesVel == 0);
   
   % Since we find row indexes for every colum (4 times the same) this has 
   % to be filtered a bit
   zeroRows = unique (zeroRows);
   zeroRowsVel = unique (zeroRowsVel);

   % Now make sure all the rows with zeros have their colums zero'd out. We don't 
   % delete these since we want to keep curCoordinates the same size as MPM
   curCoordinates(zeroRows,:) = 0;
   curCoordinatesVel(zeroRowsVel,:) = 0;

   % Calculate how many cells we have in this frame
   numberOfCells (frameCount,1) = (length(whatCells) - length(zeroRows));

   % Calculate the displacement the cells travelled from this frame to the next
   displacement (:,frameCount) = sqrt((curCoordinates(:,1)-curCoordinates(:,3)).^2 + (curCoordinates(:,2) - curCoordinates(:,4)).^2);
   displacementVel (:,frameCount) = sqrt((curCoordinatesVel(:,1)-curCoordinatesVel(:,3)).^2 + (curCoordinatesVel(:,2) - curCoordinatesVel(:,4)).^2);
   
   % Find the indexes of all rows in clustercells which are non-zero and
   % calculate the displacement of all the cells in those clusters
   realClusterRows = find (handles.clustercells(:,frameCount));
   clusterDisplacement (1:length (realClusterRows), frameCount) = displacement (handles.clustercells (realClusterRows, frameCount), frameCount);

   % Find the indexes of all rows in singlecells which are non-zero and
   % calculate the displacement of all those cells who are NOT in a cluster
   realSingleRows = find (handles.singlecells (:, frameCount));
   singleDisplacement (1:length (realSingleRows), frameCount) = displacement (handles.singlecells (realSingleRows, frameCount), frameCount);
   
   % Calculate the average displacement for all selected cells (since we sum up, the zeros do not bother us) 
   if (length(whatCells) - length(zeroRows)) ~= 0
      averageDisplacementVel (frameCount,1) = sum (displacementVel(:,frameCount)) / (length(whatCells) - length(zeroRowsVel));
      averageDisplacement (frameCount,1) = sum (displacement(:,frameCount)) / (length(whatCells) - length(zeroRows));
   else
      averageDisplacementVel (frameCount,1) = 0;
       averageDisplacement (frameCount,1) = 0;
   end
   
   % Calculate the average displacement for cells in a cluster
   if (size (find (handles.clustercells (:,frameCount)), 1)) ~= 0
      averageClusterDisplacement (frameCount, 1) = sum (clusterDisplacement(:, frameCount)) / ...
                                               (size (find (handles.clustercells (:,frameCount)), 1));
   else
      averageClusterDisplacement (frameCount, 1) = 0;
   end

   % Calculate the average displacement for cells NOT in a cluster
   if (size (find (handles.singlecells (:,frameCount)), 1)) ~= 0
      averageSingleDisplacement (frameCount, 1) = sum (singleDisplacement(:, frameCount)) / ...
                                              (size (find (handles.singlecells (:,frameCount)), 1));
   else
      averageSingleDisplacement (frameCount, 1) = 0;
   end
   
   % Calculate the variance of the displacements. Here we do have to worry about zeros
   % and that's why we remove them first using the find function (finds all non-zero
   % elements in displacement
   realDisplacement = displacement (find (displacement (:,frameCount)), frameCount);
   velocityVariance (frameCount) = var (realDisplacement);
   
   % Now we'll calculate single and cluster cell numbers via another method
   thresholdValue = averageDisplacement (frameCount,1);
   thresholdedCells (frameCount,1) = sum (realDisplacement > thresholdValue);  % Store nr of single cells
   thresholdedCells (frameCount,2) = sum (realDisplacement <= thresholdValue); % Store nr of cluster cells
   thresholdedCells (frameCount,3) = thresholdedCells (frameCount,1) / length (realDisplacement); % perc. single cells
   thresholdedCells (frameCount,4) = thresholdedCells (frameCount,2) / length (realDisplacement); % perc. cluster cells
   
   % Generate the 2D histograms and normalize the result. This is preparation for the 3D 
   % one which will be dealt with later
   if (length(whatCells) - length(zeroRows)) ~= 0
      velocityHist (:, frameCount) = hist (realDisplacement, binSize)';
      velocityHist (:, frameCount) = velocityHist (:,frameCount) ./ (length(whatCells) - length(zeroRows));
   else
      velocityHist (:, frameCount) = 0;
   end
end

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
cd (savePath);
save ('velocityVariance', 'velocityVariance');
save ('displacement', 'displacement');
save ('displacementVel', 'displacementVel');
save ('singleDisplacement', 'singleDisplacement');
save ('clusterDisplacement', 'clusterDisplacement');
save ('averageDisplacement', 'averageDisplacement');
save ('averageDisplacementVel', 'averageDisplacementVel');
save ('averageSingleDisplacement', 'averageSingleDisplacement');
save ('averageClusterDisplacement', 'averageClusterDisplacement');
save ('thresholdedCells','thresholdedCells');

% update the handles struct
guidata (hObject, handles);
