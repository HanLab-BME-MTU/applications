function [M, clusterProps, cellProps, imageCount] = ptTrackCells (ptJob, jobNumber)
% ptTrackCells finds and links coordinates in a serie of phase contrast images 
%
% SYNOPSIS       [M, clusterProps, cellProps, imageCount] = ptTrackCells (ptJob, jobNumber)
%
% INPUT          ptJob : a structure which contains the data for the job to process. See below 
%                        for the exact structure details
%                jobNumber : which job is currently being dealt with. This is used to print
%                            some status information on the matlab command line
%
% OUTPUT         M : described in ptTrackLinker
%                cellProps : 
%                  cellProps(:,1) = coord(:,1);
%	 	           cellProps(:,2) = coord(:,2);
%		           cellProps(:,3) = clustermember(:);
%                clusterProps :
%                  clusterProps (:,1) = uniqClusterNr (:);          (which cells are in this cluster)
%                  clusterProps (:,2) = numberOfCells (:);          (how many cells in this cluster)
%                  clusterProps (:,3) = clusterArea (:);            (the area of this cluster)
%                  clusterProps (:,4) = clusterPerimeter (:);       (the length of the perimeter)
%                  clusterProps (:,5) = clusterPerimeterElements (:); (how many elements does the perimeter exist of)
%                imageCount : the last image that was processed (also needed in case of crashes)
%
% DEPENDENCIES   ptTrackCells uses { imClusterSeg
%				     ptTrackLinker
%				     ptCheckMinimalCellDistance
%				     ptFindNuclei
%				     ptCalculateCellArea
%				     ptFindHalos
%				     ptFindCellsByTemplate }
%                                  
%                ptTrackCells is used by { PolyTrack }
%
% REMARK         the ptJobs structure looks as follows:
% 			imagedirectory : where are the images located
% 			imagename : what are the images called (sort of a template)
% 			imagenameslist : list of images within imagedirectory with imagename
% 			firstimage : which images shall we start with (refers to imagenameslist)
% 			lastimage : which images shall be the last one (refers to imagenameslist)
% 			intensitymax : highest value image can have (calc from bitdepth)
% 			minsize : minimal size of nuclei 
% 			maxsize : maximal size of nuclei
% 			minsdist : minimal distance between two cells
% 			increment : image to image increment (imagenameslist)
% 			noiseparameter : used to calculate threshold within ptFindCellsByTemplate for ignoring certain pixels
% 			savedirectory : where shall calculated data be saved to
% 			sizetemplate : what size should a template have
% 			boxsize : what size should the searcharea(template matching) have
% 			timestepslide : hao many timesteps should retrospectively be analysed during programm execution
% 			minedge : minimal distance to edge to still use templatesearch
% 			mincorrqualtempl : minimal quality of correlation(template) 
% 			mintrackcorrqual : minimal quality of correlation (track comparison)
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source, make function independent of GUI handle

% Tell the user that we've started
fprintf (1, 'Starting analysis of job number %d:\n', jobNumber);

% Get all the info from the ptJob object
imageDirectory       = ptJob.imagedirectory;
imageName            = ptJob.imagename;
increment            = ptJob.increment;
startFrame           = ptJob.firstimage;
endFrame             = ptJob.lastimage;
imageNameList        = ptJob.imagenameslist;
maxSearch            = ptJob.maxsearch;
saveDirectory        = ptJob.savedirectory;
percentBackground    = ptJob.noiseparameter;
sizeTemplate         = ptJob.sizetemplate;
boxSizeImage         = ptJob.boxsize;
minSizeNuc           = ptJob.minsize;
maxSizeNuc           = ptJob.maxsize;
intensityMax         = ptJob.intensityMax;
minEdge              = ptJob.minedge;
minimalQualityCorr   = ptJob.mincorrqualtempl;
minTrackCorr         = ptJob.mintrackcorrqual;
minDistCellToCell    = ptJob.minsdist;
timeStepSlide        = ptJob.timestepslide;
distanceToCellArea   = round (minDistCellToCell / 2);
erodeDiskSize        = 15;
emptyM               = zeros (1,4);      % All zeros M matrix entry
emptyCell            = zeros (1,3);
emptyCluster         = zeros (1,5);

% Check that the first and last image numbers are actually the right way around
if ~(endFrame > startFrame)
   fprintf (1, 'ptTrackCells: Last image # is smaller than first image # in job %d\n', jobNumber);
   return;
else		% lastImaNum > firstImaNum
    
   % Set up the loop counter
   loopCount = 0;

   % Initialize the storage for lost cells
   lostCells = [];

   % loop from the first frame to the last frame using the given increment
   for imageCount = startFrame : increment : endFrame
       
      % Print the image number that is being processed
      fprintf (1, '     Processing image # %d on %s\n', imageCount, datestr(now));
      
      % Increase the loop counter
      loopCount = loopCount + 1;
      
      % Read the image to be processed from disk using ptGetProcessedImage
      % Go to the directory where the images are
      if ~exist (imageDirectory, 'dir')
         fprintf (1, 'Image directory %s does not exist. Exiting...\n', imageDirectory);
         return;
      end
      cd (imageDirectory);

      % Get the filename of the image with number jImageNum
      imageFilename = char(imageNameList(imageCount));

      % Read the current image and normalize the intensity values to [0..1]
      tempImage = imreadnd2 (imageFilename, 0, intensityMax);

      % Process the image: subtract background, increase contrast between halos and nuclei
      [newImage, backgroundLevel] = ptGetProcessedImage (tempImage, intensityMax, 21);

      % Calculate the size of the image
      [imgHeight, imgWidth] = size (newImage);
      
      % Segment the image; if necessary multiple times
      nucleiArea = 1;
      escapeCount = 0; maxCount = 3;
      while ((nucleiArea > 0.5) || (haloArea > 0.5)) && (escapeCount < maxCount)

         % Calculate the variation for mu0 in case we segment more than once
         % The second time start with the initial mu0
         if escapeCount == 0
            imageMinimum = min (min (newImage));
            nucleusLevel = imageMinimum + abs(0.1 * imageMinimum);
            imageMaximum = max (max (newImage));
            haloLevel = imageMaximum - abs(0.1 * imageMaximum);
            mu0 = [nucleusLevel ; backgroundLevel ; haloLevel];
         else
            fprintf (1, 'ptTrackCells: Bad segmentation: adjusting mu0 parameter: %d %d %d\n', mu0(1), mu0(2), mu0(3));
            mu0 = mu0 + [0.03 ; -0.02 ; -0.03];
         end

         % Segment the image again
         [segmentedImage, dummy, mu0Calc] = imClusterSeg (newImage, 0, 'method', 'kmeans', 'k_cluster', 3, 'mu0', mu0);

         % The number of pixels in the nuclei should be relatively small compared to all the pixels in the image.
         % If large, something went wrong and we'll have to do the clustering again. Build in some sort
         % of mechanism that we don't hang in here forever.
         bwNuclei = zeros (size (segmentedImage));
         bwNuclei (find (segmentedImage == 1)) = 1;
         nucleiArea = bwarea (bwNuclei) / (size (segmentedImage, 1) * size (segmentedImage, 2));
         % Do the same for the halos since the total area shouldn't be to big as well
         bwHalo = zeros (size (segmentedImage));
         bwHalo (find (segmentedImage == 3)) = 1;
         haloArea = bwarea (bwHalo) / (size (segmentedImage, 1) * size (segmentedImage, 2));

         % Increase the counter so that we eventually come out of the loop
         escapeCount = escapeCount + 1;
      end   % while 

      % Find all the cell nuclei coordinates
      [nucCoord, imgNuclei] = ptFindNuclei (segmentedImage, minSizeNuc, maxSizeNuc);

      % If the user specified his/her own coordinates on the gui, use these instead
      if imageCount == startFrame
         % Get the previously initiated coordinates if they exist
         if ~isempty (ptJob.coordinatespicone);
            nucCoord = ptJob.coordinatespicone;

            % Just in case there are any zeros in there, let's remove them because they will only cause trouble
            nucCoord = nucCoord (unique (find (nucCoord(:,1) ~= 0)),:);
         end
      end

      % Find all the halo coordinates
      [haloCoord, wouldBeNucCoord, imgHalo] = ptFindHalos (segmentedImage, minSizeNuc, maxSizeNuc);

      % Make sure the minimum cell to cell distance is valid
      newCoord = ptCheckMinimalCellDistance (nucCoord, wouldBeNucCoord, minDistCellToCell);

      % Check whether we can find some more cells based on average cell size and cluster areas
      % that have no coordinates in them: if these are big enough we label them as cells
      [avgCoord, clusterImage, labeledCellImage] = ptFindCoordFromClusters (segmentedImage, newCoord, minSizeNuc);

      % If we found any new cells add them to the already found ones
      if ~isempty (avgCoord)
         newCoord = cat (1, newCoord, avgCoord);
      end

      % From frame 2 we should start matching coordinates
      if imageCount > startFrame
         % Match the current coordinates with the ones found in the previous frame
         matchedCells = ptTracker (previousCoord, newCoord, maxSearch, maxSearch);

         % Use template matching to find cells that where not found in the BMTNN match
	     unmatchedCells = find (matchedCells (:,3) == 0 & matchedCells (:,4) == 0);
	     unmatchedCellsCoord = matchedCells (unmatchedCells,1:2);
	 
         if ~isempty (unmatchedCells)
            couldNotMatch = [];
            for jCount =  1 : size (unmatchedCells, 1)
               unmatchedCellCoord = unmatchedCellsCoord(jCount,:);
               % Check whether the lost cell was near the edge of the image (specified by minEdge)
               % in this case we assume it wandered out of sight and we don't try to template match it
               if abs (unmatchedCellCoord(1,1)) > minEdge & abs (unmatchedCellCoord(1,1) - imgHeight) > minEdge & ...
                  abs (unmatchedCellCoord(1,2)) > minEdge & abs (unmatchedCellCoord(1,2) - imgWidth) > minEdge
                  
                  % Do the template matching
                  [templateCellCoord, correlation] = ptFindCellsByTemplate (unmatchedCellCoord, previousImage, ...
                                                             newImage, backgroundLevel,percentBackground, ...
                                                             sizeTemplate, boxSizeImage);
                                                                     
                  % Make sure that we only accept cells with a minimum correlation
	              if correlation > minimalQualityCorr
                   
                     % Check that the newly found cell is far enough away from
                     % other cells. If it is not disregard it.
                     if (min(sqrt ((newCoord(:,1) - templateCellCoord(1,1)).^2 + ...
                                   (newCoord(:,2) - templateCellCoord(1,2)).^2))) > minDistCellToCell

                        % Update the newCoord array (since we found a new cell after all)
                        newCoord (end+1,1) = templateCellCoord (1,1);
                        newCoord (end,2) = templateCellCoord (1,2);
                     end
                  else
                     % We didn't find a template match which means we lost the
                     % cell: add it to the lost cell list. We might find it again later
                     % so that we can close the gap in the track
                     % The current M entry is different from the current frame nr: recalculate
                     currentMEntry = imageCount;
                     if isempty (lostCells)
                        lostCells = [unmatchedCellCoord, currentMEntry];
                     else
                        lostCells = cat (1, lostCells, [unmatchedCellCoord, currentMEntry]);
                     end
                  end   % if correlation
               end   % if abs (
            end   % for jCount
         end   % if ~isempty (unmatchedCells)
         
         % Ofcourse we can do this the other way around as well: find old cells by template
         unmatchedNewCells = find (matchedCells (:,1) == 0 & matchedCells (:,2) == 0);
	     unmatchedNewCellsCoord = matchedCells (unmatchedNewCells,3:4);
	 
         if ~isempty (unmatchedNewCells)
            couldNotMatch = []; 
            for jCount =  1 : size (unmatchedNewCells, 1)
               unmatchedNewCellCoord = unmatchedNewCellsCoord(jCount,:);
               [templateCellCoord, correlation] = ptFindCellsByTemplate (unmatchedNewCellCoord, newImage, ...
                                                                         previousImage, backgroundLevel, ...
                                                                         percentBackground, sizeTemplate, ...
                                                                         boxSizeImage);
                                                                     
               % Make sure that we only accept cells with a minimum correlation
	       if correlation > minimalQualityCorr
                  if (min(sqrt ((newCoord(:,1) - templateCellCoord(1,1)).^2 + ...
                                (newCoord(:,2) - templateCellCoord(1,2)).^2))) > minDistCellToCell   

                     % Add these coordinates to the previously found ones
                     previousCoord (end+1,1) = templateCellCoord (1,1);
                     previousCoord (end,2) = templateCellCoord (1,2);
                  end
               end   % if correlation
            end   % for jCount
	     end   % if ~isempty (unmatchedCells)
         
         % Make sure the cells in newCoord and previousCoord are far enough
         % away from each other
         newCoord = ptCheckMinimalCellDistance (newCoord, [], minDistCellToCell);
         previousCoord = ptCheckMinimalCellDistance (previousCoord, [], minDistCellToCell);
         
         % Now that we have new newCoord and new PreviousCoord coordinates, we can do a renewed match
         matchedCells = ptTracker (previousCoord, newCoord, maxSearch, maxSearch);

         % Add a number of zero rows to the matchedCells matrix to make space for the coordinates
         % that we later will calculate to close tracks
         zeroRows = zeros (size (lostCells,1), 4);
         matchedCells = [matchedCells ; zeroRows];

         % Add the newly found matches to the M matrix and make sure all the rows and colums 
         % are correct for track linking
         M (1 : size (emptyM, 1), 1 : size (emptyM, 2), loopCount - 1) = emptyM;
         M (1 : size (matchedCells, 1), 1 : size (matchedCells, 2), loopCount - 1) = matchedCells;

         % There's one more thing to do: see if we can match lost cells in previous frames to any of
         % the new ones found in this frame. Allow this for a max history of timeStepSlide frames.
         newCells = matchedCells (find (matchedCells (:,1) == 0 & matchedCells (:,2) == 0 & ...
                                        matchedCells (:,3) ~= 0 & matchedCells (:,4) ~= 0),3:4);
         if ~isempty (lostCells) & ~isempty (newCells)        
            lostCellsToMatch = lostCells (find (lostCells (:,3) >= max ((loopCount - timeStepSlide), 1)), :);
            if ~isempty (lostCellsToMatch)
               matchedLostCells = ptTracker (lostCellsToMatch (:,1:2), newCells, maxSearch, maxSearch);
            
               % Kick out all entries that have not been matched (they get a chance later)
               matchedLostCells (find ((matchedLostCells (:,1) == 0 & matchedLostCells (:,2) == 0) | ...
                                       (matchedLostCells (:,3) == 0 & matchedLostCells (:,4) == 0)),:) = [];
            
               % If we did find any matches it means we can close a gap in the track of a cell
               % in the M matrix
               if ~isempty (matchedLostCells)
                  clusterDir = [saveDirectory filesep 'info'];
                  % The current M entry is different from the current frame nr: recalculate
                  currentMEntry = ceil ((imageCount - startFrame) / increment);
                  [M, lostCells] = ptCloseGapsInTracks (M, matchedLostCells, lostCells, currentMEntry, ...
                                                        startFrame, increment, clusterDir);   
               end
            end     % ~isempty (lostCellsToMatch)
         end    % ~isempty (lostCells) & ~isempty (newCells)
      end   % if imageCount > startFrame

      % Calculate single cell and cluster properties and generate binary and labeled images
      [cellProp, clusterProp] = ptCalculateCellArea (labeledCellImage, newCoord, distanceToCellArea, minSizeNuc);
      
      % Accumulate the cell properties
      cellProps (1 : size (emptyCell, 1), 1 : size (emptyCell, 2), imageCount) = emptyCell;
      cellProps (1 : size (cellProp, 1), 1 : size (cellProp, 2), imageCount) = cellProp;
   
      % Accumulate the cluster properties
      clusterProps (1 : size (emptyCluster, 1), 1 : size (emptyCluster, 2), imageCount) = emptyCluster;
      clusterProps (1 : size (clusterProp, 1), 1 : size (clusterProp, 2), imageCount) = clusterProp;
                                                  
      % Write segmented and cluster binary images, M and cell coordinates to disk
      if ~exist (saveDirectory, 'dir')
         fprintf (1, 'ptTrackCells: Save directory %s does not exist. Exiting...\n', saveDirectory);
         return;
      else
         % Go to the save directory
         cd (saveDirectory); 
         
         % Save M matrix (tracks)
         if exist ('M', 'var')
            save('M.mat', 'M');
         end

         % Save cell properties
         if exist ('cellProps', 'var')
            save ('cellProps.mat', 'cellProps');
         end

         % Save cluster properties
         if exist ('clusterProps', 'var')
            save ('clusterProps.mat', 'clusterProps');
         end
 
         % Go to the directory where the segmentation, cell and cluster info
         % is going to be saved. Create it first if it doesn't exist and cd into it
         infoDir = [saveDirectory filesep 'info'];
         if ~exist (infoDir)
            mkdir (saveDirectory, 'info');
         end
         cd (infoDir);

         % Create numerical index to number the files
         s=3; strg = sprintf('%%.%dd',s);
         indxStr = sprintf (strg, imageCount);
         clusterFile = ['clusters' indxStr];
         segmentFile = ['segmentedImg' indxStr];
         nucleiFile = ['imgNuclei' indxStr];
         
         % Save as tiff as well
         if exist ('clusterImage', 'var')
            save (clusterFile, 'clusterImage');
            imwrite (clusterImage, [saveDirectory filesep 'body' filesep clusterFile '.tif']);
         end

         % Do some processing on the segmented image first
         if exist ('segmentedImage', 'var')
            %save (segmentFile, 'segmentedImage');
            procSegmImage = zeros (size (segmentedImage));
            procSegmImage (find (segmentedImage == 2)) = 0.5; procSegmImage (find (segmentedImage == 3)) = 1;
            imwrite (procSegmImage, [saveDirectory filesep 'body' filesep segmentFile '.tif']);
         end

         % Save the nuclei image
         if exist ('imgNuclei', 'var')
            imwrite (imgNuclei, [saveDirectory filesep 'body' filesep nucleiFile '.tif']);
         end

         % Save coordinates in ascii file
         coordinatesFile = ['coordinates' indxStr];
         save ([saveDirectory filesep 'body' filesep coordinatesFile '.txt'], 'newCoord', '-ASCII');
      end   % if ~exist (saveDirector
      
      % Save the previous coordinates and image for the next run through the loop
      previousCoord = newCoord;
      previousImage = newImage;
   
   end   % for imageCount
end % if ~(lastImaNum > firstImaNum)

% Generate the cell and cluster properties. This function also takes
% coordinates into account that have been generated to close gaps
clusterDirectory = [saveDirectory filesep 'info'];
[cellProps, clusterProps] = ptCalculateCellAreaWithM (M, distanceToCellArea, minSizeNuc, clusterDirectory, ...
                                                      startFrame, endFrame, increment);
                                                  
% Generate the MPM matrix as well so that it can be used from disk if needed
MPM = ptTrackLinker (M);

% Go to the save directory
cd (saveDirectory);

% Save M matrix 
if exist ('M', 'var')
   save('M.mat', 'M');
end

% Save MPM matrix
if exist ('MPM', 'var')
   save ('MPM.mat', 'MPM');
end

% Save cell properties
if exist ('cellProps', 'var')
   save ('cellProps.mat', 'cellProps');
end

% Save cluster properties
if exist ('clusterProps', 'var')
   save ('clusterProps.mat', 'clusterProps');
end

% Tell the user we've finished
fprintf (1, '\nFinished analysis of job number %d.\n', jobNumber);
