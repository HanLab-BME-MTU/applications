function speed(hObject)
% speed plots information gathered in MPM. Certain details
%              get derived from MPM first
%
% SYNOPSIS       speed(hObject)
%
% INPUT          hObject : handle to an object within PolyTrack_PP
%                inputs are fetched directly from the GUI (PolyTrack_PP)
%                MPM
%                singlecells
%                clustercells   
%                lastimage...
%
% OUTPUT         saves plots and calc velocities to disk          
%
% DEPENDENCIES   speed  uses {nothing}
%                                  
%                speed is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04
handles=guidata(hObject);

% Starting frame relativ to first frame analysed with polytrack
startFrame = round((handles.postpro.analfirstimg - handles.jobvalues.firstimage) / handles.jobvalues.increment) + 1;

% Make sure we don't get negative frame numbers
if startFrame < 1
    startFrame = 1;
end

% Last frame relative to first frame analysed with polytrack
stopFrame = floor((handles.postpro.anallastimg - handles.jobvalues.firstimage) / handles.jobvalues.increment + 0.00001) + 1;

% Make sure we stop at the last image and not further than that one
if stopFrame > handles.jobvalues.lastimage
    stopFrame = handles.jobvalues.lastimage;
end

% Get the path where to save the results
saveAllPath = handles.postpro.saveallpath;

numberOfFrames = stopFrame - startFrame + 1;
 
% Initialize the selected-cells matrix
whatCells = [];

% See if the user specified a cell range that should be used and allocate
% the space for these
if ~isempty(handles.selectedcells)
   whatCells = zeros(size(handles.selectedcells,1));
   displacement = zeros(size(handles.selectedcells,1));
   %dispHistogram = zeros(size(handles.selectedcells,1));
else
   whatCells = zeros(size(handles.MPM,1),2);
   displacement = zeros(size(handles.MPM,1), numberOfFrames);
   %dispHistogram = zeros(size(handles.MPM,1), numberOfFrames);
end

% Prepare displacement matrix for clusters and single cells
clusterDisplacement = zeros(size(handles.clustercells));
singleDisplacement = zeros(size(handles.singlecells));

% Prepare vector that holds single and cluster cell numbers
thresholdedCells = zeros(numberOfFrames,4);

% Setup the binning used to define granularity of 3D velocity histogram
% This number is based on the max track distance
binSize = [0 : 1 : (handles.jobvalues.maxsearch + 10)];

% Now we go through every frame of the set
for iFrame = startFrame : (stopFrame - 1)

   % See if the user selected a specific range of cells and act accordingly
   if ~isempty(handles.selectedcells) 
      whatCells = handles.selectedcells;
   else  
      % If nothing has been selected use all the cells in this frame
      whatCells = [1:length(handles.MPM(:,(2*iFrame-1)))]' ;
   end
	    
   % Take four columns of MPM [x1,y1,x2,y2]. Basically every row is two
   % sets of coordinates, corresponding to two subsequent frames.
   % For every frame take these 2 x and 2 y coordinates out of MPM and assign
   % these to curCoordinates. From these the displacement can be calculated 
   % from frame to frame later on
   curCoordinates = handles.MPM (whatCells, (2*iFrame-1):(2*iFrame+2));

   % Find all the zeros in the curCoordinates matrix: these can then be filtered out
   % Note: We do not erase the rows with zeros in them. In this way,
   % the row index of displacement corresponds to the row index of MPM
   [zeroRows,cols] = find(curCoordinates == 0);

   % Since we find row indexes for every colum (4 times the same) this has 
   % to be filtered a bit
   zeroRows = unique (zeroRows);

   % Now make sure all the rows with zeros have their colums zero'd out. We don't 
   % delete these since we want to keep curCoordinates the same size as MPM
   curCoordinates(zeroRows,:) = 0;

   % Calculate how many cells we have in this frame
   numberOfCells (iFrame,1) = (length(whatCells) - length(zeroRows));

   % Calculate the displacement the cells travelled from this frame to the next
   displacement (:,iFrame) = sqrt((curCoordinates(:,1)-curCoordinates(:,3)).^2 + (curCoordinates(:,2) - curCoordinates(:,4)).^2);

   % Find the indexes of all rows in clustercells which are non-zero and
   % calculate the displacement of all the cells in those clusters
   realClusterRows = find (handles.clustercells(:,iFrame));
   clusterDisplacement (1:length (realClusterRows), iFrame) = displacement (handles.clustercells (realClusterRows, iFrame), iFrame);

   % Find the indexes of all rows in singlecells which are non-zero and
   % calculate the displacement of all those cells who are NOT in a cluster
   realSingleRows = find (handles.singlecells (:, iFrame));
   singleDisplacement (1:length (realSingleRows), iFrame) = displacement (handles.singlecells (realSingleRows, iFrame), iFrame);
   
   % Calculate the average displacement for all selected cells (since we sum up, the zeros do not bother us) 
   if (length(whatCells) - length(zeroRows)) ~= 0
      averageDisplacement (iFrame,1) = sum (displacement(:,iFrame)) / (length(whatCells) - length(zeroRows));
   else
      averageDisplacement (iFrame,1) = 0;
   end
   
   % Calculate the average displacement for cells in a cluster
   if (size (find (handles.clustercells (:,iFrame)), 1)) ~= 0
      averageClusterDisplacement (iFrame, 1) = sum (clusterDisplacement(:, iFrame)) / ...
                                               (size (find (handles.clustercells (:,iFrame)), 1));
   else
      averageClusterDisplacement (iFrame, 1) = 0;
   end

   % Calculate the average displacement for cells NOT in a cluster
   if (size (find (handles.singlecells (:,iFrame)), 1)) ~= 0
      averageSingleDisplacement (iFrame, 1) = sum (singleDisplacement(:, iFrame)) / ...
                                              (size (find (handles.singlecells (:,iFrame)), 1));
   else
      averageSingleDisplacement (iFrame, 1) = 0;
   end
   
   % Calculate the variance of the displacements. Here we do have to worry about zeros
   % and that's why we remove them first using the find function (finds all non-zero
   % elements in displacement
   realDisplacement = displacement (find (displacement (:,iFrame)), iFrame);
   velocityVariance (iFrame) = var (realDisplacement);
   
   % Now we'll calculate single and cluster cell numbers via another method
   thresholdValue = averageDisplacement (iFrame,1);
   thresholdedCells (iFrame,1) = sum (realDisplacement > thresholdValue);  % Store nr of single cells
   thresholdedCells (iFrame,2) = sum (realDisplacement <= thresholdValue); % Store nr of cluster cells
   thresholdedCells (iFrame,3) = thresholdedCells (iFrame,1) / length (realDisplacement); % perc. single cells
   thresholdedCells (iFrame,4) = thresholdedCells (iFrame,2) / length (realDisplacement); % perc. cluster cells
   
   % Generate the 2D histograms and normalize the result. This is preparation for the 3D 
   % one which will be dealt with later
   if (length(whatCells) - length(zeroRows)) ~= 0
      velocityHist (:, iFrame) = hist (realDisplacement, binSize)';
      velocityHist (:, iFrame) = velocityHist (:,iFrame) ./ (length(whatCells) - length(zeroRows));
   else
      velocityHist (:, iFrame) = 0;
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
% hgsave (h_fig, [saveAllPath filesep '3dVelocityHistogram.fig']);
% print (h_fig, [saveAllPath filesep '3dVelocityHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep '3dVelocityHistogram.tif'], '-dtiff');
% 
% % Generate a 3D velocity histogram using 3-d bars which is chopped up in bins with size binSize
% h_fig = figure;
% bar3 (binSize, velocityHist);
% title ('3D velocity histogram using bins');
% 
% % And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [saveAllPath filesep '3dBinnedVelocityHistogram.fig']);
% print (h_fig, [saveAllPath filesep '3dBinnedVelocityHistogram.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep '3dBinnedVelocityHistogram.tif'], '-dtiff');
% 
% % Generate a plot with the number of cells found per frame
% h_fig = figure, plot (numberOfCells), title ('Number of cells');
% xlabel ('Frames');
% ylabel ('# Cells');
% ymax = max (numberOfCells);
% axis ([startFrame stopFrame 0 ymax]);
% 
% % And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [saveAllPath filesep 'numberOfCells.fig']);
% print (h_fig, [saveAllPath filesep 'numberOfCells.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep 'numberOfCells.tif'], '-dtiff');
% 
% Generate a plot with both average velocity and variance of single cells and 
% cells within clusters in one figure
% h_fig = figure;
% subplot (2,1,1); plot (averageClusterDisplacement), title ('Average velocity of cells within clusters');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max (averageClusterDisplacement);
% axis ([startFrame stopFrame 0 ymax]);
% subplot (2,1,2); plot (averageSingleDisplacement), title('Average velocity of single cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max(averageSingleDisplacement);
% axis ([startFrame stopFrame 0 ymax]);

% And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [saveAllPath filesep 'velocitySingleCluster.fig']);
% print (h_fig, [saveAllPath filesep 'velocitySingleCluster.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep 'velocitySingleCluster.tif'], '-dtiff');

% Generate a plot with both average velocity and variance of cells in one figure
% h_fig = figure;
% subplot (2,1,1); plot (averageDisplacement), title('Average velocity of cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max(averageDisplacement);
% axis ([startFrame stopFrame 0 ymax]);
% subplot (2,1,2); plot (velocityVariance), title('Variance of avarage velocity of cells');
% xlabel ('Frames');
% ylabel ('variance');
% ymax = max(velocityVariance);
% axis ([startFrame stopFrame 0 ymax]);

% And save this figure to disk as fig, eps and tiff
% hgsave (h_fig, [saveAllPath filesep 'velocityVariance.fig']);
% print (h_fig, [saveAllPath filesep 'velocityVariance.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep 'velocityVariance.tif'], '-dtiff');

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
% for jFrame = startFrame : (stopFrame - 1)
%    dispHistogram (:,jFrame) = hist (displacement (find (displacement(:,jFrame))), binSize);
% end   
% h_fig = figure;
% plot (dispHistogram, 'b'), title('Histogram of velocity of cell');
% ylabel ('Frames');
% xlabel ('# Cells');
% save ('dispHistogram','dispHistogram');

%xmax = max(displacement);
%axis ([0 xmax startFrame stopFrame]);
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
cd (saveAllPath);
save ('velocityVariance', 'velocityVariance');
save ('displacement', 'displacement');
save ('singleDisplacement', 'singleDisplacement');
save ('clusterDisplacement', 'clusterDisplacement');
save ('averageDisplacement', 'averageDisplacement');
save ('averageSingleDisplacement', 'averageSingleDisplacement');
save ('averageClusterDisplacement', 'averageClusterDisplacement');
save ('thresholdedCells','thresholdedCells');

% update the handles struct
guidata (hObject, handles);
