function ptTrackCells (ptJob, jobNumber)
% ptTrackCells finds and links coordinates in a serie of phase contrast images 
%
% SYNOPSIS       ptTrackCells (ptJob, jobNumber)
%
% INPUT          ptJob : a structure which contains the data for the job to process. See below 
%                        for the exact structure details
%                jobNumber : which job is currently being dealt with. This is used to print
%                            some status information on the matlab command line
%
% OUTPUT         All outputs are written directly to disk 
%                M : described in ptTrackLinker
%                cellProps : 
%                  cellProps(:,1) = coord(:,1);
%	 	   cellProps(:,2) = coord(:,2);
%		   cellProps(:,3) = belongsto(:);  (number of cluster - label)
%		   cellProps(:,4) = numberOfOccurences(:);  (how many cells in the cluster this cell is in)
%		   cellProps(:,5) = bodycount(:);  (area of the cluster with the number given in belongsto)
%		   cellProps(:,6) = perimdivare(:);  (cluster)
%                nameClust : the binary images of the areas occupied by cells (per frame)
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
firstImageNr         = ptJob.firstimage;
lastImageNr          = ptJob.lastimage;
imageNameList        = ptJob.imagenameslist;
maxSearch            = ptJob.maxsearch;
saveDirectory        = ptJob.savedirectory;
percentBackground    = ptJob.noiseparameter;
sizeTemplate         = ptJob.sizetemplate;
boxSizeImage         = ptJob.boxsize;
minSizeNuc           = ptJob.minsize;
maxSizeNuc           = ptJob.maxsize;
intensityMax    = ptJob.intensityMax;
minEdge              = ptJob.minedge;
minimalQualityCorr   = ptJob.mincorrqualtempl;
minTrackCorr         = ptJob.mintrackcorrqual;
minDistCellToCell    = ptJob.minsdist;
approxDistance       = round (minDistCellToCell / 2);
erodeDiskSize        = 15;
emptyM               = zeros (1,4);      % All zeros M matrix entry

% Check that the first and last image numbers are actually the right way around
if ~(lastImageNr > firstImageNr)
   fprintf (1, 'ptTrackCells: Last image # is smaller than first image # in job %d\n', jobNumber);
   return;
else		% lastImaNum > firstImaNum
    
   % Set up the loop counter
   loopCount = 0;

   % Initialize the storage for lost cells
   lostCells = [];

   % loop from the first frame to the last frame using the given increment
   for imageCount = firstImageNr : increment : lastImageNr
       
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
      if imageCount == firstImageNr
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

      if imageCount > firstImageNr
         % Match the current coordinates with the ones found in the previous frame
         matchedCells = fsmTrackTrackerBMTNN (previousCoord, newCoord, maxSearch, maxSearch);

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
                     if isempty (lostCells)
                        lostCells = [unmatchedCellCoord, imageCount-1];
                     else
                        lostCells = cat (1, lostCells, [unmatchedCellCoord, imageCount-1]);
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
         
         % Now that we have new newCoord and new PreviousCoord coordinates, we can do a renewed match
         matchedCells = fsmTrackTrackerBMTNN (previousCoord, newCoord, maxSearch, maxSearch);

         % Add a number of zero rows to the matchedCells matrix to make space for the coordinates
         % that we later will calculate to close tracks
         zeroRows = zeros (size (lostCells,1), 4);
         matchedCells = [matchedCells ; zeroRows];

         % Add the newly found matches to the M matrix and make sure all the rows and colums 
         % are correct for track linking
         M (1 : size (emptyM, 1), 1 : size (emptyM, 2), loopCount - 1) = emptyM;
         M (1 : size (matchedCells, 1), 1 : size (matchedCells, 2), loopCount - 1) = matchedCells;

         % There's one more thing to do: see if we can match lost cells in previous frames to any of
         % the new ones found in this frame. Allow this for a max history of 1 frames.
         newCells = matchedCells (find (matchedCells (:,1) == 0 & matchedCells (:,2) == 0 & ...
                                        matchedCells (:,3) ~= 0 & matchedCells (:,4) ~= 0),3:4);
         if ~isempty (lostCells) & ~isempty (newCells)        
            lostCellsToMatch = lostCells (find (lostCells (:,3) >= max ((imageCount - 1), 1)), :);
            if ~isempty (lostCellsToMatch)
               matchedLostCells = fsmTrackTrackerBMTNN (lostCellsToMatch (:,1:2), newCells, maxSearch, maxSearch);
            
               % Kick out all entries that have not been matched
               matchedLostCells (find ((matchedLostCells (:,1) == 0 & matchedLostCells (:,2) == 0) | ...
                                       (matchedLostCells (:,3) == 0 & matchedLostCells (:,4) == 0)),:) = [];
            
               % If we did find any matches it means we can close a gap in the track of a cell
               % in the M matrix
               if ~isempty (matchedLostCells)
                  clusterDir = [saveDirectory filesep 'info'];
                  [M, lostCells] = ptCloseGapsInTracks (M, matchedLostCells, lostCells, imageCount-1, clusterDir);   
               end
            end     % ~isempty (lostCellsToMatch)
         end    % ~isempty (lostCells) & ~isempty (newCells)
      end   % if imageCount > firstImageNr

      % Calculate single cell and cluster properties and generate binary and labeled images
      [cellProps, clusterProps, clusterImage, labeled] = ptCalculateCellArea ...
                                                        (segmentedImage, newCoord, imgNuclei, ...
                                                         imgHalo, approxDistance, minSizeNuc, maxSizeNuc);

      % Write segmented and cluster binary images, M and cell coordinates to disk
      if ~exist (saveDirectory, 'dir')
         fprintf (1, 'ptTrackCells: Save directory %s does not exist. Exiting...\n', saveDirectory);
         return;
      else
         cd (saveDirectory);

         % Save M matrix (tracks), cell and cluster properties
         if exist ('M', 'var')
	        save('M.mat', 'M');
	     end
	     if exist ('cellProps', 'var')
            save ('cellProps.mat', 'cellProps')
	     end
	     if exist ('clusterProps', 'var')
            save ('clusterProps.mat', 'clusterProps')
	     end

         % Now we go to the directory where the segmentation and cluster info
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
            save (segmentFile, 'segmentedImage');
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

% Tell the user we've finished
fprintf (1, 'Finished analysis of job number %d.\n', jobNumber);
