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
%				     ptFindTemplateCells }
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
% 			fi_nucleus / la_nucleus : nucloi intensity first / last image
% 			fi_background / la_background : background intensity first / last image
% 			fi_halolevel / la_halolevel : halo intensity first image
% 			leveladjust : factor to adjust intensity difference nucloi/ background
% 			minsize : minimal size of nucloi 
% 			maxsize : maximal size of nucloi
% 			minsdist : minimal distance between two cells
% 			minmaxthresh : onoff - should minima and segmentation be used 
% 			clustering : onoff - should clustering be used
% 			increment : image to image increment (imagenameslist)
% 			noiseparameter : used to calculate threshold within ptFindTemplateCells for ignoring certain pixels
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

% Assign all the job data to local variables
imageDirectory       = ptJob.imagedirectory;
imageName            = ptJob.imagename;
increment            = ptJob.increment;
firstImaNum          = ptJob.firstimage;
lastImaNum           = ptJob.lastimage;
imageNamesList       = ptJob.imagenameslist;
levNucFirst          = ptJob.fi_nucleus;
levBackFirst         = ptJob.fi_background;
levHaloFirst         = ptJob.fi_halolevel;      
levNucLast           = ptJob.la_nucleus;
levBackLast          = ptJob.la_background;
levHaloLast          = ptJob.la_halolevel;      

% Range in which direct assignement of found coordinates from frame to frame is accepted
maxSearch            = ptJob.maxsearch;
saveDirectory        = ptJob.savedirectory;
percentBackground    = ptJob.noiseparameter;
sizeTemplate         = ptJob.sizetemplate;
boxSizeImage         = ptJob.boxsize;
 
% Minimal/maximal size of the black spot in the cells
minSizeNuc           = ptJob.minsize;
maxSizeNuc           = ptJob.maxsize;

% Get the maximum image intensity of the images in the job
maxImageIntensity    = ptJob.intensityMax;

% An educated guess of the minimal distance between neighbouring cells
% (better to big than to small, for this value is used for the static
% search. Searches with templates are more tolerant)
minDistCellToCell    = ptJob.minsdist; 

% This value influences the level between the nuclei and the background, as
% calculated from input (clicking in the pictures, which pop up shortly
% after the program starts
levelAdjust          = ptJob.leveladjust;

% Minimum 4
howManyTimeStepSlide = ptJob.timestepslide;

% Minimal distance to edge for tracking with template
minEdge              = ptJob.minedge;
minimalQualityCorr   = ptJob.mincorrqualtempl;
minTrackCorr         = ptJob.mintrackcorrqual;

% Segmentation method used
if ptJob.clustering
   method = 1;
elseif ptJob.minmaxthresh
   method = 2;
elseif ptJob.emclustering
   method = 3;
else
   fprintf(1, 'ptTrackCells: Wrong method. A valid method (1,2,3) has to be selected.\n');
   return;
end

% How much the blobs found in ptFindHalos shall be eroded. This is an indirect
% size criteria for ptFindHalos. Increase - minimal size of halos will be
% increased, decrease - ... decreased
erodeDiskSize = round ((sqrt (minSizeNuc)) / 2) ;

% Range within which to look for correlations of tracks over several pics. 
distanceCorrel = maxSearch * 1.5; 

% Within which distance shall the program body look for an area, to which
% it can allocate a cell (coordinates of a cell)
approxDistance = round(minDistCellToCell / 2);

% The following markers identify the origins of coordinates. They get added to the
% values of coordinates. Depending on how the coordinate was found, a certain marker 
% will be added. In this way, the coordinates themselves always carry all information.

% If you wish, you can add other markers, to distinguish between other classes of cells.
% DO NOT USE 0.9 AS A MARKER!!! (x <= 0.8 is ok), because the rounding will change the 
% coordinates themselves (which is not good)
% If you wish to mark cells with more then one property, markers placed at the second 
% digit (0.01 , 0.02 ...) could be used.
goodCellMarker = 0.1;
goodCell = round (goodCellMarker * 10);

templateCellMarker = 0.2;
templateCell = round (templateCellMarker * 10);

newCellMarker = 0.3;
newCell = round (newCellMarker * 10);

newTemplateCellMarker = 0.4;
newTemplateCell = round (newTemplateCellMarker * 10);

fileExtension = imageName(end-3:end);
bodyFileName = imageName(1:end-7);

howManyImg = lastImaNum - firstImaNum + 1;

levDiffFirst = abs (levNucFirst - levBackFirst) * levelAdjust;
levDiffLast = abs (levNucLast - levBackLast) * levelAdjust;

incrDiff = (levDiffLast - levDiffFirst) / (howManyImg - 1);
incrBack = (levBackLast - levBackFirst) / (howManyImg - 1);
incrHalo = (levHaloLast - levHaloFirst) / (howManyImg - 1);
incrNuc  = (levNucLast  - levNucFirst)  / (howManyImg - 1);

i = 0;
k = 0;
a = 0;
lostOnEdge = 0;

newImage = [];
MPM = [];
cellProps = [];
M = [];
coord = [];

emptyM = zeros(1,4);

% Go to the directory where the images are
if exist (imageDirectory, 'dir')
    cd (imageDirectory);
else
    fprintf (1, 'Image directory %s does not exist. Exiting...\n', imageDirectory);
    return;
end

% Check that the first and last image numbers are actually the right way around
if ~ (lastImaNum > firstImaNum)
   fprintf (1, 'Last image # is smaller than first image # in job %s\n', jobNumber);
   return
else    
   % Index is equal to the number of the image, countLoops keeps
   % track of the loops
   countLoops = 0;

   % Generate the initial mu0 values for the segmentation
   mu0Init = [levNucFirst ; levBackFirst ; levHaloFirst];
   
   % Loop through all the images of the job using the specified increment
   for jImageNum = firstImaNum:increment:lastImaNum
        
      % Let the user know where he/she is in the sequence by printing the number
      countLoops = countLoops + 1;
      fprintf (1, '    Image number = %d\n', jImageNum);
        
      % Go to the directory where the images are
      if exist (imageDirectory, 'dir')
          cd (imageDirectory); 
      else
          fprintf (1, 'Image directory %s does not exist. Exiting...\n', imageDirectory);
          return;
      end

      % Get the filename of the image with number jImageNum
      fileName = char(imageNamesList(jImageNum));
        
      % Read the current image and normalize the intensity values to [0..1]
      tempImage = imreadnd2 (fileName, 0, maxImageIntensity);

      % Subtract the background from the image (this is done using robustfit)
      [newImage, background] = ptGetProcessedImage (tempImage, 20);
      clear tempImage;
      
      % Calculate the size of the normalized image
      [img_h,img_w] = size(newImage);
        
      % If image is the first one in the sequence do the following:
      if jImageNum == firstImaNum
         % Prepare the space for the binary image
         binaryImage = zeros(img_h, img_w);
                      
         % Determine the intensity level of the halos
         haloLevel = (levHaloFirst - levBackFirst) * 2 / 3 + levBackFirst;
            
         % In case the clustering algorithm has been selected use the kmeans method
         % to segment the image
         if method == 1
            % Here's the call to the image segmentation function; we also
            % get mu0 back which we use to segment the following images
            [segmentedImage, dummy, mu0] = imClusterSeg (newImage, 0, 'method', 'kmeans', 'k_cluster', 3, ...
                                                         'mu0', mu0Init);
            %[segmentedImage, dummy, mu0] = imClusterSeg (newImage, 0, 'method', 'kmeans', 'k_cluster', 3);
            inputImage = segmentedImage;
         elseif method == 2
            inputImage = newImage;
         elseif method == 3
            [segmentedImage, dummy, p, mu0] = imClusterSeg (newImage, 0, 'method', 'em', 'k_min', 2, 'k_max', 2);
            inputImage = newImage;
         end

         [nucCoord, imgNuclei] = ptFindNuclei (inputImage, levDiffFirst, minSizeNuc, maxSizeNuc, method);
         [haloCoord, imgHalo] = ptFindHalos (inputImage, erodeDiskSize, haloLevel, method);
                               
         if ~isempty (ptJob.coordinatespicone);
            newCoord = ptJob.coordinatespicone;
            
            % Just in case there are any zeros in there, let's remove them
            % because they will only cause trouble later on
            tempIndex = unique (find (newCoord(:,1) ~= 0));
            newCoord = newCoord (tempIndex,:);
         else 
            newCoord = ptCheckMinimalCellDistance (nucCoord, haloCoord, minDistCellToCell);           
         end

         % Calculate cell area and other properties and label the found areas
         [cellProps, binaryImage, labeled] = ptCalculateCellArea (inputImage, newCoord, imgNuclei, imgHalo, approxDistance, method);

         initProps = zeros (1, size (cellProps, 2));
                  
         % Now we mark these coordinates as good coordinates, meaning their
         % existence is credible. Futhermore cells marked as good cells will
         % be treated with priority, when it comes to allocating cells found in
         % a new picture
         previousCoord = newCoord + goodCellMarker;
                  
         emptyM = zeros (1,4);
                      
      % The following is what we do with images that are not the first
      else
          
         % Adjust level via increment
         levDiffInterpol = levDiffFirst + (countLoops - 1) * increment * incrDiff;
         levelBackInterpol = levBackFirst + (countLoops - 1) * increment * incrBack;
         levelHaloInterpol = levHaloFirst + (countLoops - 1) * increment * incrHalo;
         %levelNucInterpol = levNucFirst + (countLoops - 1) * increment * incrNuc;
         haloLevel = (levelHaloInterpol - levelBackInterpol) * 2 / 3 + levelBackInterpol;
         
         if method == 1
            % Segment the image; if necessary multiple times
            nucleiArea = 1;
            escapeCount = 0; maxCount = 5;
            while (nucleiArea > 0.5) && (escapeCount < maxCount)
                
                % Calculate the variation for mu0 in case we segment more than once
                % The second time start with the initial mu0
                if escapeCount == 1
                    mu0 = mu0Init;
                elseif escapeCount > 1
                    mu0 = mu0 + [0.005 ; 0.005 ; 0.005];
                end
                
                % Segment the image again
                [segmentedImage, dummy, mu0] = imClusterSeg (newImage, 0, 'method', 'kmeans', 'k_cluster', 3, 'mu0', mu0);
                inputImage = segmentedImage;

                % The number of pixels in the nuclei should be relatively small compared to all the pixels in the image.
                % If large, something went wrong and we'll have to do the clustering again. Build in some sort
                % of mechanism that we don't hang in here forever.
                bwNuclei = zeros (size (segmentedImage));
                bwNuclei (find (segmentedImage == 1)) = 1;
                nucleiArea = bwarea (bwNuclei) / (size (segmentedImage, 1) * size (segmentedImage, 2));

                % Increase the counter so that we eventually come out of the loop
                escapeCount = escapeCount + 1;
           end

            % Did we run against the max count and not segment correctly?: Stop
            if escapeCount == maxCount
               error ('ptTrackCells: error, could not segment image correctly.');
               return;
            end

         elseif method == 2
            inputImage = newImage;
            haloLevel = (levHaloFirst - levBackFirst) * 2 / 3 + levBackFirst;
         elseif method == 3
            [segmentedImage, dummy, p0, mu0] = imClusterSeg (newImage, 0, 'method', 'em', 'k_min', 2, 'k_max', 2, 'mu0', mu0, 'p0', p0);
         end

         % Find the coordinates of all the nuclei
         [nucCoord, imgNuclei] = ptFindNuclei (inputImage, levDiffInterpol, minSizeNuc, maxSizeNuc, method);

         % Find the coordinates of all the halos
         [haloCoord, imgHalo] = ptFindHalos (inputImage, erodeDiskSize, haloLevel, method);
         
         % Ensure a minimal distance between cells
         newCoord = ptCheckMinimalCellDistance (nucCoord, haloCoord, minDistCellToCell);  
         
         % Mark coordinates as high quality coordinates. These cells are found in the usual way, 
         % so they are superior to cells found for instance by templates
         newCoord = round (newCoord) + goodCellMarker;

         % Initialize some variables
         prevCoordNoMarker = [];
         coordMarkerInfo = [];
         coordIndex = [];
         previousCoordNewCells = [];
         newCellCoordIndex = [];
         newCellCoord = [];
         notAllocatedCellCoord = [];

         % Here we seperate the cells found already long ago from the ones found 
         % later (we believe in the existence of the first more than in the 
         % existence of the later)
         prevCoordNoMarker = floor (previousCoord + 0.1);
         
         % Save the marker information
         coordMarkerInfo = round (10 * (previousCoord - prevCoordNoMarker));
         
         % Find all the coordinates which belong to template cells or good
         % cells
         coordIndex =  find (coordMarkerInfo(:,1) == templateCell | coordMarkerInfo(:,1) == goodCell);

         % Old cells
         previousCoordNewCells = previousCoord (coordIndex,:);

         clear coordIndex;
         clear prevCoordNoMarker;

         %new cells
         newCellCoordIndex =  find (coordMarkerInfo(:,1) == newCell | coordMarkerInfo(:,1) == newTemplateCell);
         newCellCoord = previousCoord (newCellCoordIndex,:);
                       
         clear coordMarkerInfo;
         clear newCellCoordIndex;

         % First we try to allocate the old cells to any set of coordinates found in the new picture
         matchedCells = fsmTrackTrackerBMTNN (previousCoordNewCells, newCoord, maxSearch, maxSearch);
                      
         clear previousCoordNewCells;
                  
         % If there are new cells, now is the time to try to allocate them
         % to any set of coordinates not yet allocated to an old cell
         if size (find (matchedCells(:,1) == 0 & matchedCells(:,2) == 0), 1) > 0 & size (newCellCoord, 1) > 0
            foundNewCells = find (matchedCells(:,1) == 0 & matchedCells(:,2) == 0);
            % In any case the sets coordinates not allocated to an old cell
            % are marked as a new cells, for they will either be allocated
            % to a new cell (found a few pictures ago), or actually become a new cell
            notAllocatedCellCoord(:,1) = round (matchedCells (foundNewCells,3)) + newCellMarker;
            notAllocatedCellCoord(:,2) = round (matchedCells (foundNewCells,4)) + newCellMarker;

            % Erase the unallocated set of coordinates from matchedCells,
            % because they will turn up again in newMatchedCells no matter what happens
            matchedCells (foundNewCells,:) = [];
            clear foundNewCells

            newMatchedCells = [];
            newMatchedCells = fsmTrackTrackerBMTNN (newCellCoord, notAllocatedCellCoord, maxSearch, maxSearch);
                                  
            % Combine the two allocation matrices to one
            revisedMatchedCells = cat (1, matchedCells, newMatchedCells);
            clear newMatchedCells;
         else 
            % Re-mark the newly found cells (the ones that have no corresponding coords in column 
            % 1 and 2 of matchedCells)as new cells with the newCell marker
            foundNewCells = find (matchedCells(:,1) == 0);
            matchedCells (foundNewCells,3:4) = round (matchedCells (foundNewCells,3:4)) + newCellMarker;
                                  
            clear foundNewCells;
                                  
            % If we have any new cell coordinates ....
            if size (newCellCoord, 1) > 0
               tempZeros = zeros (size (newCellCoord,1), 4);
               tempZeros(:,1:2) = newCellCoord(:,1:2);
               revisedMatchedCells = cat (1, matchedCells, tempZeros);
               clear tempZeros;
            else
               % Since their are no new cells newMatchedCells is empty, so matchedCells is complete
               revisedMatchedCells = matchedCells;
            end
         end

         if length (find ((revisedMatchedCells(:,3) ~= 0) & (revisedMatchedCells(:,4) ~= 0))) ~= size (newCoord, 1)
            lostNewCell = 1;
         end

         if length( find( (revisedMatchedCells(:,1)~=0) & (revisedMatchedCells(:,2)~=0))) ~= size(previousCoord,1)
            lostOldCell = 1
         end
                       
         clear matchedCells;
         clear newCellCoord;
         clear notAllocatedCellCoord;
                         
         % So far so good, but we still have cells that weren't found in the new picture. Instead of just giving up, we now
         % try to track them with templates. Of course these cells will be marked (within the tracking via templates routine - called 
         % ptFindTemplateCells). If they are old cells they will be marked as "old cells tracked via template", if they are new
         % cells as "new cells tracked via template". The actual marker is a digit (fraction) behind the dot.
         unmatchedCells = [];
         unmatchedCells = find (revisedMatchedCells (:,3) == 0 & revisedMatchedCells (:,4) == 0);
                           
         if ~isempty (unmatchedCells)
            for jCount =  1 : size (unmatchedCells, 1)
               unmatchedCellsCoord = [0,0];
               unmatchedCellsCoord = [revisedMatchedCells(unmatchedCells(jCount),1), revisedMatchedCells(unmatchedCells(jCount),2)];
				
               % first we check if the lost cell was located
               % near the edge (distance to edge smaller than
               % minEdge). If so, we asume the cell has wandered
               % out of the picture and we stop following it
               if abs(unmatchedCellsCoord(1,1)) < minEdge | abs(unmatchedCellsCoord(1,1)-img_w) < minEdge | abs(unmatchedCellsCoord(1,2)) < minEdge | abs(unmatchedCellsCoord(1,2)-img_h) < minEdge
                  % one never knows what kind of data
                  % somebody may wish
                  % to possess
                  lostOnEdge = lostOnEdge + 1;
               else
                  correlation = [];
                  templateCellCoord = [];
                    
                  % this is the finding via templates
                  % routine. It will return the found
                  % coordinates and the value of the
                  % corralation
                  [templateCellCoord, correlation] = ptFindTemplateCells (unmatchedCellsCoord, previousImage, newImage, levelBackInterpol, templateCellMarker, ...
                                           newCell, newTemplateCell, newTemplateCellMarker, percentBackground, sizeTemplate, boxSizeImage);
                                                      
                  if sqrt ((unmatchedCellsCoord(1,1) - templateCellCoord(1,1)).^2 + (unmatchedCellsCoord(1,2) - templateCellCoord(1,2)).^2) > 2 * maxSearch
                     % Stop: the program made a mistake
                     break;
                  end
                                                      
                  % Of course we only accept the correlation, if it of a minimal
                  % quality. NOTE: that within the routine the background (everything 
                  % that doesn't belong to a cell) of the template and of the searcharea
                  % is replaced by random noise. So a background - background correlation
                  % will average zero
                            
                  if correlation > minimalQualityCorr
                     % Now we calculate the minimal distance between the cell found by correlation and any other cell
                     cellDistance = [];
                     minDistIndex = [];
                     
                     [cellDistance, minDistIndex] =  min (sqrt ((revisedMatchedCells(:,3) - templateCellCoord(1,1)).^2 + (revisedMatchedCells(:,4) - templateCellCoord(1,2)).^2));
                               
                     if cellDistance < round (minDistCellToCell / 1.5)
                                             
                        % The following is just extracting the information of the coordinates found via template 
                        % and the coordinates that are very close to it
                        tmpCoordNoMarker = [];
                        tmpCoordMarkerInfo = [];
                        newCoordNoMarker = [];
                                    
                        tmpCoordNoMarker = floor (revisedMatchedCells (minDistIndex, 3) + 0.1);
                        tmpCoordMarkerInfo = round (10 * (revisedMatchedCells (minDistIndex, 3) - tmpCoordNoMarker));
                        
                        clear tmpCoordNoMarker;
                       
                        newCoordNoMarker = floor (templateCellCoord(1, 1) + 0.1);
                        newCoordMarkerInfo = round (10 * (templateCellCoord(1,1) - newCoordNoMarker));
                        
                        clear newCoordNoMarker;
                         
                        % If the coordinates found via template belong to an old cell, we trust them, otherwise 
                        % the cell won't be propagated
                        if newCoordMarkerInfo == templateCell
                                                   
                           % If the cell close by is a new cell, we allocate it's coordinates to the 
                           % cell we try to track by template (which is an old cell - we've checked that above) and
                           % erase the new cell
                           if tmpCoordMarkerInfo == newCell
                              revisedMatchedCells (unmatchedCells(jCount), 3) = round (revisedMatchedCells (minDistIndex, 3)) + goodCellMarker;
                              revisedMatchedCells (unmatchedCells(jCount), 4) = round (revisedMatchedCells (minDistIndex, 4)) + goodCellMarker;
                              revisedMatchedCells (minDistIndex, 3:4) = 0;
                                        
                               % If it is a new cell, itself propagated via templates, we erase it and allocate the
                               % coordinates we have (template) found to the old cell
                           elseif tmpCoordMarkerInfo == newTemplateCell
                              revisedMatchedCells (minDistIndex,3:4) = 0;
                              revisedMatchedCells (unmatchedCells(jCount),3) = templateCellCoord (1,1);
                              revisedMatchedCells (unmatchedCells(jCount),4) = templateCellCoord (1,2);
                           end
                        end

                        clear newCoordMarkerInfo;
                        clear tmpCoordMarkerInfo;
                         
                        % If the minimal distance to the nearest cell is big enough,we allocate the coordinates we 
                        % have (template) found to the old cell
                        revisedMatchedCells (unmatchedCells(jCount), 3) = templateCellCoord (1,1);
                        revisedMatchedCells (unmatchedCells(jCount), 4) = templateCellCoord (1,2);
                     end
                     clear cellDistance;
                     clear minDistIndex;
                  end
                  clear templateCellCoord;
                  clear correlation;
               end
               clear unmatchedCellsCoord;
            end
         end
                       
         % I know all these 'end' are confusing, so for orientation: we 
         % have now got the matrice (revisedMatchedCells), which defines the allocations between the
         % the cells of the last picture and the cells of the picture
         % we are looking at at the moment.
                       
         clear unmatchedCells;
                          
         % Cosmetics
         keptCells = find (ismember (revisedMatchedCells (:,1:4), [0 0 0 0], 'rows') - 1);
         revisedMatchedCells = revisedMatchedCells (keptCells,:);
                      
         clear keptCells;
                     
         M (1 : size (emptyM, 1), 1 : size (emptyM, 2), countLoops - 1, 1) = emptyM;
         %stuff this information into a stack
         M (1 : size (revisedMatchedCells, 1), 1 : size (revisedMatchedCells, 2), countLoops - 1, 1) = revisedMatchedCells;
                          
         if length (find ((revisedMatchedCells (:,1) ~= 0) & (revisedMatchedCells(:,2) ~= 0))) ~= size (previousCoord, 1)
            lostOldCell = 1;
         end 
                       
         clear revisedMatchedCells;
         revisedMatchedCells = emptyM;
                       
         % This is a quite elaborate and (if I may say so) sophisticated
         % part of the programm, which basically changes the tracks of
         % certain cells, based on an analysis of the various tracks of
         % cells over the last 'howManyTimeStepSlide' frames (defined in PolyTrack (GUI)). 
                   
         if countLoops > howManyTimeStepSlide + 1
                           
            templateRow = [];
            MPMslide = [];
            MPMflo = [];
            xtempl = [];
            ytempl = [];

            % MPMslide is a matrix which gives the tracks of all cells over
            % the last 'howManyTimeStepSlide' frames. It gets updated after every new frame.
            MPMslide = ptTrackLinker (M (:, :, (countLoops - howManyTimeStepSlide + 1) : (countLoops - 1)));

            MPMflo = floor (MPMslide + 0.1);
            MPMIdentCell = round (10 * (MPMslide - MPMflo));

            clear MPMflo;

            % After a few frames, new cells will be accepted as old cells, but only if, 
            % after their initial finding, they weren't propagated by templates 
            foundNewCells = find (MPMIdentCell (:,3) == newCell);
            
            if ~isempty (foundNewCells)
               for jnewCell = 1 : length (foundNewCells)
                  % Make sure a new cell wasn't propagated via templates
                  if isempty (find (MPMIdentCell (foundNewCells (jnewCell),:) == newTemplateCell)) & ...
                                    MPMIdentCell (foundNewCells (jnewCell), end) ~= 0;
                     markOld = [];
                     markOld = find (M (:, 3, countLoops - 1) == MPMslide (foundNewCells (jnewCell), end-1) & ...
                                     M(:, 4, countLoops - 1) == MPMslide (foundNewCells (jnewCell), end));
                     M (markOld, 3:4, countLoops - 1) = round (M (markOld, 3:4, countLoops - 1)) + goodCellMarker;
                     MPMslide (foundNewCells (jnewCell), end-1) = M (markOld, 3, countLoops - 1);
                     MPMslide (foundNewCells (jnewCell), end) = M (markOld, 4, countLoops - 1);
                                          
                     clear markOld;
                                     
                  elseif  MPMIdentCell (foundNewCells (jnewCell), end) == newTemplateCell & ...
                          MPMIdentCell (foundNewCells (jnewCell), end-2) == newTemplateCell
                     delCell = [];
	             delCell = find (M (:, 3, countLoops - 1) == MPMslide (foundNewCells (jnewCell), end-1) & ...
                                     M (:, 4, countLoops - 1) == MPMslide (foundNewCells (jnewCell), end));
	             M (delCell, 3:4, countLoops - 1) = 0;
	             MPMslide (foundNewCells (jnewCell), end-1) = 0;
                     MPMslide (foundNewCells (jnewCell), end) = 0;

                     clear delCell
                  end
               end
            end
                  
            clear foundNewCells;

            % Here we search for old cells, that have been propagated via
            % template search over the last five frames. Afterwards we will
            % compare their tracks with tracks of newly found cells. If they
            % corralate strongly, we replace the templatestuff (old cell) by the newly
            % found cell
            templateRow = find (ismember (MPMIdentCell (:, 3:end), templateCell * ones (1, length (MPMIdentCell (1, 3:end))), 'rows'));
                  
            if ~isempty (templateRow)
               % Newly found cells
               newCellRow = find (MPMIdentCell (:,3) == newCell & MPMIdentCell (:,end) == newCell);
                      
               if ~isempty (newCellRow)
                           
                  if ~isempty (find (MPMIdentCell (newCellRow, 3:end)))
                     BIGISHTROUBLE = 1;
                  end
                                
                  xtempl = [];
                  ytempl = [];
                             
                  % Extract the x and y coordinates of the cells         
                  xtempl (:,:) = (MPMslide (templateRow(:), 3:2:end-1))'; 
                  ytempl (:,:) = (MPMslide (templateRow(:), 4:2:end))';
                         
                  % We want to correlate the displacement vectors (displacement of a cell from one frame to the next), 
                  % so we need to calculate it
                  moveXXTemplate = diff (xtempl);
                  moveYYTemplate = diff (ytempl);
                                
                  for indexnewCell = 1 : length (newCellRow)
                                                      
                     xnewone = [];
                     ynewone = [];
                           
                     moveXXnewCell = [];
                     moveYYnewCell = [];
                                                
                     xnewone = (MPMslide (newCellRow (indexnewCell), 3:2:end-1))';
                     ynewone = (MPMslide (newCellRow (indexnewCell), 4:2:end))';
                     
                     moveXXnewCell = diff (xnewone);
                     moveYYnewCell = diff (ynewone);
             
                     % Make sure we correlate two things of the same length
                     if size (xtempl,1) == length (xnewone)
                        corrTempRow = [];
                        nomore = [];
                        for indexTempRow = 1 : length (templateRow)
                                      
                           lostfound = [];
                           tempfound = [];
                           writeHereRow = [];
                                          
                           distance = min (sqrt ((xtempl (:,indexTempRow) - xnewone(:)).^2 + (ytempl (:,indexTempRow) - ynewone(:)).^2));
                                           
                           % Now, if these points at some point are really close - exchange them
                           if distance < round (minDistCellToCell / 2)
                                                         
                              % this is the substitution of the cell tracked by templates by the newly found cell. 
                              % It's kind of tricky, because we have to change the values not in MPMslide, but in M. 
                              % Since M is not sorted, we have to always find the right spot
                              templateCoord = [];
                              newCeCoord = [];
                                                       
                              % First we copy out the coordinates of the old cell, tracked by templates
                              templateCoord = MPMslide (templateRow (indexTempRow),:);
                              % Now we erase the values from MPMslide, so that we won't by accident use them again
                              MPMslide (templateRow (indexTempRow),:) = 0;
                                                     
                              % Same procedure for the new cell we want to link to the old cell
                              newCeCoord = MPMslide (newCellRow (indexnewCell),:);
                              MPMslide (newCellRow (indexnewCell),:) = 0;
                                                  
                              for jFrameCount = 1 : ((size (MPMslide,2) - 2) / 2)
                                                          
                                 writeHereRow = [];
                                 eraseHereRow = [];
                                         
                                 % Ok,now: in templateCoord is the track of the old cell.newCeCoord is the track of the cell 
                                 % we wish to put at it's place. Both are of the form: [x1,y1,x2,y2,x3,y3,x4,y4,x5,y5] 
                                 % So always two numbers correspond to the coordinates in one frame(1,2,3...) mean different 
                                 % frames. So we find the location of the coordinates of the old cell, by searching elements 
                                 % of templateCoord in M. With every increase in jFrameCount we increase the index of 
                                 % templateCoord by two (because two values (x and y)per picture). M has three dimensions. 
                                 % Every element of the third dimension is a matrice(other two dimensions)of the form [x1,y1,x2,y2],
                                 % thus giving the connection between two frames, BUUUUTTTT: the next element will be 
                                 % [x2,y2,x3,y3]!!!!!! This is very important!!!!!!.  So with every frame (jFrameCount) we 
                                 % make one step into the third dimesion in M.% with each set of x and y of templateCoord we 
                                 % locate the rigth row in M. 
                                 writeHereRow = find (M (:, 1, countLoops - 5 + jFrameCount) == templateCoord (2 * jFrameCount - 1) & ...
                                                      M (:, 2 ,countLoops - 5 + jFrameCount) == templateCoord ( 2 * jFrameCount));
                                                           
                                 % Since we overwrite the old cell with the new, we have to eliminate the new, 
                                 % in order not to have it twice
                                 eraseHereRow = find (M (:, 1, countLoops - 5 + jFrameCount) == newCeCoord (2 * jFrameCount - 1) & ...
                                                      M (:, 2, countLoops - 5 + jFrameCount) == newCeCoord(2 * jFrameCount));
                                  
                                 % Now we juggle the information into the right place. Since M is [x1,y1,x2,y2] , 
                                 % next layer [x2,y2,x3,y3] and newCeCoord is [x1,y1,x2,y2,x3,y3,x4,y4,x5,y5] therefor:  
                                 % M (writeHereRow, :, countLoops-4) replaced by newCeCoord (1:4) and M (writeHereRow, :, countLoops-3) 
                                 % replaced by newCeCoord(3:6) and so on.
                                 M (writeHereRow, :, countLoops - 5 + jFrameCount) = newCeCoord (2 * jFrameCount - 1 : 2 * jFrameCount + 2);
                                                       
                                 % Erase redundant info
                                 M (eraseHereRow, :, countLoops - 5 + jFrameCount) = 0;
                              end

                              clear writeHereRow;
                              clear eraseHereRow;
                              clear templateCoord;
                              clear newCeCoord;

                              nomore = 1;
                              break;
                                                    
                           % Here too, we have to stay within a maximal distance
                           elseif distance < distanceCorrel 
                              % We correlate the displacement vectors of the old cell with the ones of the new cell
                              correlx = xcorr (moveXXTemplate (:, indexTempRow), moveXXnewCell(:), 'coeff');
                              correly = xcorr (moveYYTemplate (:, indexTempRow), moveYYnewCell(:), 'coeff');
                              
                              indCorr = round (length (correlx) / 2);
                                                            
                              % One correlation value for each row of templates
                              corrTempRow (indexTempRow) = correlx (indCorr) + correly (indCorr);
                           end
                        end
                              
                        % Nomore says if this new cell was already asigned to an old cell by distance criteria
                        % if this isn't the case, it's time to see, if there is an old cell (nearby as checked above)
                        % that correlates nicely
                        if isempty (nomore)
                           if ~isempty (corrTempRow);
                              indBC = [];
                              [bestCorr, indBC] = max (corrTempRow);
                                                   
                              if bestCorr > minTrackCorr
                                 % this is exactly the same procedure as above. I know it would be clever to make a function out of it
                                 % but there would be so many BIG matrices involved I decided to have it in here twice.
                                 templateCoord = [];
                                 newCeCoord = [];
                                                                      
                                 templateCoord = MPMslide (templateRow (indBC), :);
                                 MPMslide (templateRow (indBC), :) = 0;
                                                                        
                                 newCeCoord = MPMslide (newCellRow (indexnewCell), :);
                                 MPMslide (newCellRow (indexnewCell), :) = 0;
                                                                      
                                 for jFrameCount = 1 : (size (MPMslide, 2) - 2) / 2
                                    writeHereRow = [];
                                    eraseHereRow = [];
                                  
                                    writeHereRow = find (M (:, 1, countLoops - 5 + jFrameCount) == templateCoord (2 * jFrameCount - 1) & ...
                                                         M (:, 2, countLoops - 5 + jFrameCount) == templateCoord (2 * jFrameCount));
                                    eraseHereRow = find (M (:, 1, countLoops - 5 + jFrameCount) == newCeCoord (2 * jFrameCount - 1) & ...
                                                         M (:, 2, countLoops - 5 + jFrameCount) == newCeCoord (2 * jFrameCount));
                                          
                                    M (writeHereRow, :, countLoops - 5 + jFrameCount) = (newCeCoord (2 * jFrameCount - 1:2 * jFrameCount + 2));
                                                                 
                                    M (eraseHereRow, :, countLoops - 5 + jFrameCount) = 0;
                                 end
                                 clear writeHereRow;
                                 clear eraseHereRow;
                                 clear templateCoord;
                                 clear newCeCoord;
                              end
                              clear indBC;
                           end    
                        end
                        clear corrTempRow;
                        clear nomore;
                     else
                        % The displacement vector of old cell and new cell (moveXXTemplate,moveXXnewCell...)
                        % aren't of the same length, so check what's wrong.
                        error = 1;
                     end
                     clear xnewone;
                     clear ynewone;
                   
                     clear moveXXnewCell;
                     clear moveYYnewCell;
                  end
               end
            end
            clear templateRow;
            clear MPMslide;
            clear xtempl;
            clear ytempl;
         end
       
         %clear newCoord;
          
         previousCoord = [];
         coordIndex = [];
         coordIndex = find (M (:, 3, countLoops - 1) | M (:, 4, countLoops - 1));
         
         previousCoord (:,1) = M (coordIndex, 3,countLoops - 1);
         previousCoord (:,2) = M (coordIndex, 4,countLoops - 1);
         clear coordIndex;
          
         tempProps = [];
                           
         if method == 1
            inputImage = segmentedImage;
         elseif method == 2
            inputImage = newImage;
         elseif method == 3
            inputImage = newImage;
         end
                 
         [tempProps, binaryImage, labeled] = ptCalculateCellArea (inputImage, previousCoord, imgNuclei, imgHalo, approxDistance, method);

         cellProps (1 : size (initProps, 1), 1 : size (initProps, 2), countLoops, 1) = initProps;
         cellProps (1 : size (tempProps, 1), 1 : size (tempProps, 2), countLoops, 1) = tempProps;
                 
         clear tempProps;
      end
      
      % Keep the previous image for the next loop
      previousImage = newImage;
          
      % We save all the data we found to the results directory.
      if exist (saveDirectory, 'dir')
          cd (saveDirectory);
      else
          fprintf (1, 'ptTrackCells: directory %s does not exist. Exiting...\n', saveDirectory);
          return;
      end
		
      save('M.mat', 'M')
      save ('cellprops.mat', 'cellProps')
      
      % Now we go to the directory where the segmentation and cluster info
      % is saved. Create it first if it doesn't exist
      bodyDirectory = [saveDirectory filesep 'body'];
      if exist (bodyDirectory)
          cd ('body');
      else
          mkdir (saveDirectory, 'body');
          cd ('body');
      end
      
      % Format
      s=3;
      strg = sprintf('%%.%dd',s);

      % Create numerical index
      indxStr = sprintf (strg, jImageNum);
      nameClust = ['clusters' indxStr]; 
      segmented = ['segmentedImg' indxStr];
      save (nameClust, 'binaryImage');
      save (segmented, 'segmentedImage');

      % Save as tiff as well
      imwrite (binaryImage, [saveDirectory filesep 'body' filesep nameClust '.tif']);
      
      % Do some processing on the segm image first
      if exist ('segmentedImage', 'var')
          procSegmImage = zeros (size (segmentedImage));
          procSegmImage (find (segmentedImage == 2)) = 0.5; procSegmImage (find (segmentedImage == 3)) = 1;
          imwrite (procSegmImage, [saveDirectory filesep 'body' filesep segmented '.tif']);
      end
            
      % Save coordinates in ascii file for Yang
      coordinates = ['coordinates' indxStr];
      saveCoord = floor (newCoord + 0.1);
      save ([saveDirectory filesep 'body' filesep coordinates '.txt'], 'saveCoord', '-ASCII');
   end
end

% Tell the user we've finished
fprintf (1, 'Finished analysis of job number %d.\n', jobNumber);
         


% AK: the poem below is kept because Colin has worked so hard on it :-)

%PLEASE NOTE WELL:
%by changing these values, 
%you have nothing to gain.
%no profit, no revenues,
%only TROUBLE that will rain down on your brain.
%By now, my friend, it should be quite plain:
%If you change these values ... you are insane.
