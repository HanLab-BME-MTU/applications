function [M, clusterProps, cellProps, frameProps, imageCount, validFrames] = ptTrackCells (ptJob, jobNumber, saveIntResults)
% ptTrackCells finds and links coordinates in a serie of phase contrast images 
%
% SYNOPSIS       [M, clusterProps, cellProps, imageCount] = ptTrackCells (ptJob, jobNumber)
%
% INPUT          ptJob : a structure which contains the data for the job to process. See below 
%                        for the exact structure details
%                jobNumber : which job is currently being dealt with. This is used to print
%                            some status information on the matlab command line
%                saveIntResults : 1 if intermediate results have to be
%                                 saved; 0 if not saved
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
%                frameProps :
%                  frameProps (:,1) = average area all cells/clusters
%                  frameProps (:,2) = average area convex hull around clusters
%                  frameProps (:,3) = average ratio area/convex-hull-area of clusters
%                imageCount : the last image that was processed (also needed in case of crashes)
%             imreadnd2   validFrames: an array showing which frames were processed (1) and which frames 
%                             were bad and therefore not processed (0)
%
% DEPENDENCIES   ptTrackCells uses { imClusterSeg 
%				     ptTrackLinker
%				     ptCheckMinimalCelptTrackCellslDistance
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
% 			mintracklength : minimal length of track in MPM
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source, make function independent of GUI handle
% Andre Kerstens        Jun 04          Modification necessary because of extra parameter ptGetProcessedImage
% Andre Kerstens        Jul 04          In case of too many bad segmentations, the function will return 
%                                       and give up on the job
% Andre Kerstens        Jul 04         ptTrackCells Added frame properties to functions
% Andre Kerstens        Jul 04          Segmentation function now uses precalculated mu0 values
% Johan de Rooij        May 05          added binning to imClusterSeg
%                                       needs to be added to GUI!!
%                                       use ptImClusterSeg now!!
%                                       set at line 118 now, default = 2.
% Johan de Rooij        Jul 05          bugfix: restart muO calc after bad
%                                       frames. needs revision, but works.
% Johan de Rooij        Jul 05          Improve segmentation by using
%                                       edgeImage for background loss.


% Get a pointer to the chromdynMain gui
hPtMain = findall(0,'Tag','polyTrack_mainwindow','Name','PolyTrack');
if ~isempty(hPtMain)
    % Get the currently selected project and experiment data
    handles = guidata(hPtMain);
end 

% Tell the user that we've started
fprintf (1, 'ptTrackCells: Starting analysis of job number %d:\n', jobNumber);

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
minTrackLength       = ptJob.mintracklength;
minDistCellToCell    = ptJob.minsdist;
timeStepSlide        = ptJob.timestepslide;
distanceToCellArea   = round (minDistCellToCell / 2);
timePerFrameGUI      = ptJob.timeperframe;
erodeDiskSize        = 15;
bgKernel             = 21;
edgeKernel           = 5;
emptyM               = zeros (1,4);      % All zeros M matrix entry
emptyCell            = zeros (1,3);
emptyCluster         = zeros (1,5);
emptyFrame           = zeros (1,5);
clustSize            = str2num(get(handles.GUI_st_txtNumberKclustere,'String'));
clustBinSize         = 1/str2num(get(handles.GUI_st_txtBinSize_txNew,'String')); %1/binSize, because the image must be binSize smaller

% Created info directory if needed
if saveIntResults
    mkdir(saveDirectory,'info');
end

% Initialize the vector that holds the valid frame info and the
% accompanying frame rates (taken from the image file header)
% We'll start with zeros so that in case a whole movie breaks of, we still
% have some frames to postprocess
frameCount = ceil((endFrame - startFrame + 1) / increment);
validFrames = zeros(2,frameCount);

% Check that the first and last image numbers are actually the right way around
if ~(endFrame > startFrame)
   fprintf (1, '     ptTrackCells: Last image # is smaller than first image # in job %d\n', jobNumber);
   return;
else		% lastImaNum > firstImaNum
    
   % Set up the loop and M counters
   loopCount = 0;
   mCount = 0;
     
   % Initialize the storage for lost cells
   lostCells = [];
   
   % Initialize bad frame counter
   badFrameCounter = 0;

   % loop from the first frame to the last frame using the given increment
   for imageCount = startFrame : increment : endFrame
       
      % Print the image number that is being processed
      fprintf (1, '     ptTrackCells: Processing frame # %d on %s\n', imageCount, datestr(now));
      
      % Increase the loop counter
      loopCount = loopCount + 1;
      
          
      % Read the image to be processed from disk using ptGetProcessedImage
      % Go to the directory where the images are
      if ~exist (imageDirectory, 'dir')
         fprintf (1, '     ptTrackCells: Image directory %s does not exist. Exiting...\n', imageDirectory);
         return;
      end
      cd (imageDirectory);

      % Get the filename of the image with number jImageNum %JvR: do you
      % mean imageCount instead of jImageNum?
      imageFilename = char(imageNameList(imageCount));
      
      % Get the information on the image
      imageInfo = imfinfo (imageFilename, 'tif');
               
      % Read the current image and normalize the intensity values to [0..1]
      tempImage = imreadnd2 (imageFilename, 0, intensityMax);
       
      %JvR; Binning if necessary
      if clustBinSize<1
        tempImage = imresize(tempImage, clustBinSize, 'bicubic');
      end
      
      % Check the variance of the whole image. If it's a very small value
      % (< 0.005 %JvR: Johan changed in to <0.001 because in our images we get small variance as well) 
      %it will be a bad frame (lots of black pixels when the
      % shutter was completely closed)
      varImage = max(var(tempImage)); 
      
      if varImage >0.0001    % 0.001
          % Process the image: subtract background, increase contrast between halos and nuclei
          [newImage, backgroundLevel, edgeImage] = ptGetProcessedImage (tempImage, intensityMax, ...
                                                                        bgKernel, edgeKernel);
            
          tempImage=[]; %JvR Delete tempImage, this matrix will not be used anymore
               
          % Calculate the size of the image
          [imgHeight, imgWidth] = size (newImage);

          % Segment the image; if necessary multiple times
          nucleiArea = 1;
          equality = 0.1;
          escapeCount = 0; maxCount = 3;
         %%%% while ((nucleiArea > 0.5) || (haloArea > 0.5) || (equality < 0.8)) && (escapeCount < maxCount) 
                
             % Segment the image (again)
             % JvR:
             % First we would like to divide the image in cell and
             % background. Then the dust and noise pixels the background as well
             % as for the cells must be deleted which is done by the bwareaopen function is used.
              
             % first substract background from image:
             % calculate background:
             edgeImageOpened = imerode (edgeImage, strel ('disk', 1));
             edgeImageClosed = imfill (edgeImageOpened, 'holes');
             edgeImageSeg = imerode (edgeImageClosed, strel ('disk', 5));
             edgeImageSeg = bwareaopen (edgeImageSeg, minSizeNuc);
             
             
             %JvR code: should work but the minSizeNuc is not a good
             %criterea
             %edgeImage=imerode(edgeImage,strel('disk',1));
             %edgeImageSeg=(edgeImage-1).^2; %reverse the 0 and 1 in the matrix
             %edgeImageSeg = bwareaopen (edgeImageSeg, round(((minSizeNuc)^2)*3.14)); %Delete small dust from the cells
             %edgeImageSeg=(edgeImageSeg-1).^2; %reverse the 0 and 1 in the matrix
             %edgeImageSeg = bwareaopen (edgeImageSeg, round(((minSizeNuc)^2)*3.14)); %Delete small dust from background 
             
             H = fspecial('disk',2);
             BlurImage = imfilter(newImage,H,'replicate');

             % noBlurImagew determine the input vector for segmentation:
             segVec = BlurImage(find(edgeImageSeg)); 
             newImage=[]; 
             % and do the segmentation on just this vector!!
             % JvR: the number of segmentations is now the variable
             % clustersize which is set by a slider of the Number of
             % K-clusters
             % JvR; Ive added the binning option, which is set by the
             % Binning slider
             [segmentedVec, dummy, mu0Calc] = ptImClusterSeg (segVec, 0, 'k_cluster', clustSize);
             if isempty(segmentedVec) % if the segmentation did not work then go out of the loop and make the frame invalid
               varImage=1;
             else
                 segmentedImage = zeros(size(edgeImageSeg)); 
                 segmentedImage(find(edgeImageSeg)) = segmentedVec;
                 DumbackgroundNumber=ceil((clustSize/2)+0.5);
                 segmentedImage(find(edgeImageSeg == 0)) = DumbackgroundNumber;
                 edgeImageSeg=[]; % JvR; Delete the image out of the memory   
                 segVec=[]; % JvR; Delete the array out of the memory 
                 segmentedVec=[]; % JvR; Delete the array out of the memory

                 % we may need a little check here, to make sure that we are
                 % not finding bulshit further down the line, because this
                 % resulted in black dots at the wrong spots.
                 % therefore we quickly determine the nuclei it finds in the
                 % segmented image: 
                 [nucCoord, imgNuclei] = ptFindNuclei (segmentedImage, minSizeNuc, maxSizeNuc); 

                 % now we compare this frame to the previous (if there was one)
                 if mCount > 0
                     testImage = imgNuclei - previmgNuclei;
                     % if they are equal, than the difference between the two
                     % should be quite small, so many zeros in the resulting
                     % matrix.
                     equality = length(find(testImage == 0))/numel(segmentedImage);
                     % I think it should be over 80% equal..this stringency can
                     % be increased later!
                 else
                     equality = 1;
                 end
                 % and store for comparison with the next frame:
                 testImage=[]; %JvR: Delete testImage out of memory
                 previmgNuclei=[]; %JvR: Delete testImage out of memory

                 % The number of pixels in the nuclei should be relatively small compared to all the pixels in the image.
                 % If large, something went wrong and we'll have to do the clustering again. Build in some sort
                 % of mechanism that we don't hang in here forever.
                 bwNuclei = zeros (size (segmentedImage));
                 bwNuclei (find (segmentedImage == 1)) = 1;
                 nucleiArea = bwarea (bwNuclei) / (size (segmentedImage, 1) * size (segmentedImage, 2));
                 bwNuclei=[]; %jvR: Delete bwNuclei matrix out of memory
                 % Do the same for the halos since the total area shouldn't be to big as well
                 bwHalo = zeros (size (segmentedImage));
                 bwHalo (find (segmentedImage == clustSize )) = 1; %standart clusterSize=5
                 haloArea = bwarea (bwHalo) / (size (segmentedImage, 1) * size (segmentedImage, 2));
                 bwHalo=[]; %jvR: Delete bwNuclei matrix out of memory

                 % Increase the counter so that we eventually come out of the loop
                 escapeCount = escapeCount + 1;
             %%%% end   % while              
             end           
      end  % if varImage > 0.001
       
      % In case we found a bad frame, mark it in the valid frame array
      if varImage <= 0.0001 | escapeCount >= maxCount
         fprintf (1, '       ptTrackCells: Bad frame %d: continue with next.\n', imageCount);
         
         % Increase bad frame counter
         badFrameCounter = badFrameCounter + 1;
         
         if badFrameCounter > 10   % Stop since we wouldn't be able to track anymore
            fprintf (1, '     ptTrackCells: Too many (>10) bad frames in a row.\nCheck the bitdepth.\nExiting...\n');
            M = [];
            clusterProps = [];
            cellProps = [];
            frameProps = [];
            imageCount = [];
            validFrames = [];
            return;
         end
         
      else   % A good frame was found         
          % Reset the bad frame counter
          badFrameCounter = 0;
             
          % Increase the M counter
          mCount = mCount + 1;
          
          % We just processed a valid frame
          validFrames(1,mCount) = loopCount;
      
          % Calculate the amount of time (in secs) between this and the previous frame
          % datenum gives back the difference in days (time is the
          % fraction) which we have to convert back to seconds
          try
              DateTime = imageInfo.DateTime;
          catch
              % Tell the user there is no time-info in files
              fprintf (1, 'these tiffs have no timestamp, job: %d.\n', jobNumber);
          end
          if mCount > 1
              %some images don't come with this piece of information..
             if exist ('DateTime') == 1;
                 % There might be several formats of storing the date; let's
                 % try a couple of them
                 julianTime = datenum(imageInfo.DateTime, 'yyyy:mm:dd HH:MM:SS');
                %julianTime = datenum(imageInfo.DateTime, 'mmm dd yy HH:MM:SS');
                frameTime = round(abs(julianTime - prevJulianTime) * 60 * 60 * 24);
         
                if frameTime >= 1  % 1 sec is the minimum we accept
                validFrames(2,mCount) = frameTime;
                else
                    validFrames(2,mCount) = timePerFrameGUI;
                end
         
             prevJulianTime = julianTime;
             else
                 validFrames(2,mCount) = timePerFrameGUI;
             end
          else
              if exist ('DateTime') == 1;
                 prevJulianTime = datenum(imageInfo.DateTime, 'yyyy:mm:dd HH:MM:SS');
                 %prevJulianTime = datenum(imageInfo.DateTime, 'mmm dd yy HH:MM:SS');
              end
          end 
          
          % Find all the cell nuclei coordinates
          % we did this already!!
          % [nucCoord, imgNuclei] = ptFindNuclei (segmentedImage, minSizeNuc, maxSizeNuc);

          % Store the nuclei coords as intermediate result
          if saveIntResults
              save([saveDirectory filesep 'info' filesep 'nucleiCoord_pass1_' num2str(imageCount) '.mat'],'nucCoord');
          end
          
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
          [haloCoord, wouldBeNucCoord, imgHalo] = ptFindHalos (segmentedImage, minSizeNuc, maxSizeNuc,clustSize);

          % Make sure the minimum cell to cell distance is valid
          newCoord = ptCheckMinimalCellDistance (nucCoord, wouldBeNucCoord, minDistCellToCell);
 
          % Store the nuclei coords incl wouldbe ones as intermediate result
          if saveIntResults
              save([saveDirectory filesep 'info' filesep 'nucleiCoord_pass2_' num2str(imageCount) '.mat'],'newCoord');
          end
          
          % Check whether we can find some more cells based on average cell size and cluster areas
          % that have no coordinates in them: if these are big enough we label them as cells
          [avgCoord, clusterImage, labeledCellImage] = ptFindCoordFromClusters (edgeImage, newCoord, ...
                                                                                minSizeNuc, edgeKernel);
          
          edgeImage=[]; %JvR: Delete edgeImage out of memory
          % If we found any new cells add them to the already found ones
          if ~isempty (avgCoord)
             newAvgCoord = cat (1, newCoord, avgCoord);
          else
             newAvgCoord = newCoord;
          end

          % Again make sure the minimum cell to cell distance is valid
          newCoord = ptCheckMinimalCellDistance (newAvgCoord, [], minDistCellToCell);

          % Store the nuclei coords incl new ones as intermediate result
          if saveIntResults
              save([saveDirectory filesep 'info' filesep 'nucleiCoord_pass3_' num2str(imageCount) '.mat'],'newCoord');
          end
          
          % From frame 2 we should start matching coordinates
          if mCount > 1

             % Match the current coordinates with the ones found in the previous frame
             if isempty(newCoord) %JvR if newCoord is empty, it should have a dimention of 2!!
               newCoord=[0 0]; 
             end
             matchedCells = ptTracker (previousCoord, newCoord, maxSearch, maxSearch);
             
             % Use template matching to find cells that where not found in the BMTNN match
             if ~isempty(matchedCells)
                 unmatchedCells = find (matchedCells (:,3) == 0 & matchedCells (:,4) == 0);
                 unmatchedCellsCoords = matchedCells (unmatchedCells,1:2);
                 
                 % Make sure the minimum cell to cell distance is valid for
                 % unmatched cells as well
                 unmatchedCellsCoord = ptCheckMinimalCellDistance (unmatchedCellsCoords, [], minDistCellToCell);
             else
                 unmatchedCellsCoord = [];
                 unmatchedCells = [];
             end

             % Keep lost cells for later
             % AK: added this to test how algoritm performs without
             % template matching
             if ~isempty(unmatchedCellsCoord)
                 currentMEntry = mCount-1.*ones(size(unmatchedCellsCoord,1),1);
                 if isempty (lostCells)
                    lostCells = [unmatchedCellsCoord, currentMEntry];
                 else
                    lostCells = cat (1, lostCells, [unmatchedCellsCoord, currentMEntry]);
                 end   
             end
             
             % Add a number of zero rows to the matchedCells matrix to make
             % space for the coordinates that we later will calculate to close
             % tracks and to add any previous coords we found
             %zeroRows = zeros (size (lostCells,1) + size (foundPrevCells,1), 4);
             % AK: the above line was replaced with the next one to test
             % algorithm without templ matching
             zeroRows = zeros (size (lostCells,1) + size (unmatchedCellsCoord,1), 4);
             matchedCells = [matchedCells ; zeroRows];

             % Add the newly found matches to the M matrix and make sure all the rows and colums 
             % are correct for track linking
             M (1 : size (emptyM, 1), 1 : size (emptyM, 2), mCount - 1) = emptyM;
             M (1 : size (matchedCells, 1), 1 : size (matchedCells, 2), mCount - 1) = matchedCells;

             % See if we can match lost cells in previous frames to any of
             % the new ones found in this frame. Allow this for a max history of timeStepSlide frames.
             newCells = matchedCells (find (matchedCells (:,1) == 0 & matchedCells (:,2) == 0 & ...
                                            matchedCells (:,3) ~= 0 & matchedCells (:,4) ~= 0),3:4);
             if ~isempty (lostCells) & ~isempty (newCells)        
                lostCellsToMatch = lostCells (find (lostCells (:,3) >= max ((mCount - timeStepSlide), 1)), :);
                if ~isempty (lostCellsToMatch)
                   if size(lostCellsToMatch,1) < size(newCells,1)
                       matchedLostCells = ptTracker (newCells, lostCellsToMatch(:,1:2), maxSearch, maxSearch);
                       matchedLostCellsTmp = matchedLostCells(:,3:4);
                       matchedLostCells(:,3:4) = matchedLostCells(:,1:2);
                       matchedLostCells(:,1:2) = matchedLostCellsTmp;
                   else
                       matchedLostCells = ptTracker (lostCellsToMatch(:,1:2), newCells, maxSearch, maxSearch);
                   end

                   % Kick out all entries that have not been matched (they get a chance later)
                   matchedLostCells (find ((matchedLostCells (:,1) == 0 & matchedLostCells (:,2) == 0) | ...
                                           (matchedLostCells (:,3) == 0 & matchedLostCells (:,4) == 0)),:) = [];

                   % If we did find any matches it means we can close a gap in the track of a cell
                   % in the M matrix
                   if ~isempty (matchedLostCells)
                      clusterDir = [saveDirectory filesep 'info'];
                      % The current M entry is different from the current frame nr: recalculate
                      %currentMEntry = ceil ((imageCount - startFrame) / increment);
                      currentMEntry = mCount-1;
                      [M, lostCells] = ptCloseGapsInTracks (M, matchedLostCells, lostCells, currentMEntry, ...
                                                            startFrame, increment, clusterDir, validFrames);   
                   end
                end     % ~isempty (lostCellsToMatch)
             end    % ~isempty (lostCells) & ~isempty (newCells)

% AK: commented out the next block to test algorithm without templ matching
%              % In case we found previous coords by template, we should modify a
%              % zero row in the previous M entry as well, otherwise we will get
%              % problems linking tracks later on 
%              if mCount > 2
%                 foundPrevCells (1,:) = [];
%                 if ~isempty (foundPrevCells)
%                    zeroInd = find (M (:,1,mCount - 2) == 0 & M (:,2,mCount - 2) == 0 & ...
%                                    M (:,3,mCount - 2) == 0 & M (:,4,mCount - 2) == 0);
%                    M (zeroInd(1:size (foundPrevCells,1)), 1 : size (foundPrevCells, 2), mCount - 2) = ...
%                       foundPrevCells;
%                 end
%              end
% 
%              % Empty foundPrevCells for the next loop run
%              clear foundPrevCells; clear zeroInd;      

          end   % if imageCount > startFrame

          % Calculate single cell and cluster properties and generate binary and labeled images
          [cellProp, clusterProp, frameProp] = ptCalculateCellAreaUsingVariance (labeledCellImage, newCoord, ...
                                                                   distanceToCellArea, minSizeNuc, edgeKernel);
          labeledCellImage=[]; % JvR: Delete matrix out of memory
          % Accumulate the cell properties
          cellProps (1 : size (emptyCell, 1), 1 : size (emptyCell, 2), imageCount) = emptyCell;
          cellProps (1 : size (cellProp, 1), 1 : size (cellProp, 2), imageCount) = cellProp;
   
           % Accumulate the cluster properties
          clusterProps (1 : size (emptyCluster, 1), 1 : size (emptyCluster, 2), imageCount) = emptyCluster;
          clusterProps (1 : size (clusterProp, 1), 1 : size (clusterProp, 2), imageCount) = clusterProp;

          % Accumulate the frame properties
          frameProps (1 : size (emptyFrame, 1), 1 : size (emptyFrame, 2), imageCount) = emptyFrame;
          frameProps (1 : size (frameProp, 1), 1 : size (frameProp, 2), imageCount) = frameProp;

          % Write segmented and cluster binary images, M and cell coordinates to disk
          if ~exist (saveDirectory, 'dir')
             fprintf (1, '     ptTrackCells: Save directory %s does not exist. Exiting...\n', saveDirectory);
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
             haloFile = ['imgHalo' indxStr];

             % Save as tiff as well
             if exist ('clusterImage', 'var')
                save (clusterFile, 'clusterImage');
                imwrite (clusterImage, [saveDirectory filesep 'body' filesep clusterFile '.tif']);
             end

              % Save intermediate results if needed    
              if saveIntResults    
                 % Do some processing on the segmented image first
                 if exist ('segmentedImage', 'var')
                    %save (segmentFile, 'segmentedImage');
                    procSegmImage = segmentedImage/clustSize ;
                    imwrite (procSegmImage, [saveDirectory filesep 'body' filesep segmentFile '.tif']);
                 end

                 % Save the nuclei image
                 if exist ('imgNuclei', 'var')
                    imwrite (imgNuclei, [saveDirectory filesep 'body' filesep nucleiFile '.tif']);
                 end

                 % Save Halo image
                 if exist ('imgHalo', 'var')
                    imwrite (imgHalo, [saveDirectory filesep 'body' filesep haloFile '.tif']);
                 end
                 
                 
                 % Save coordinates in ascii file
                 coordinatesFile = ['coordinates' indxStr];
                 save ([saveDirectory filesep 'body' filesep coordinatesFile '.txt'], 'newCoord', '-ASCII');
             end
          end   % if ~exist (saveDirector

          % Save the previous coordinates and image for the next run through the loop
          previousCoord = newCoord;
          previmgNuclei= imgNuclei;
          
      end  % if varImage <= 0.001 | escapeCount >= maxCount
   end   % for imageCount
end % if ~(lastImaNum > firstImaNum)

% Remove the zero entries from the validFrames array
tempValidFrames = validFrames (:, find (validFrames(1,:) ~= 0));
validFrames = tempValidFrames;
                                                  
MPM = ptTrackLinker (M);

save ('MPMBeforeProcessing.mat', 'MPM');

% Make sure the tracks have a minimum track length (specified on GUI)
MPM = ptMinTrackLength (MPM, minTrackLength);

% Generate the cell and cluster properties. This function also takes
% coordinates into account that have been generated to close gaps
clusterDirectory = [saveDirectory filesep 'info'];
[cellProps, clusterProps, frameProps] = ptCalculateCellAreaWithMPM (MPM, distanceToCellArea, minSizeNuc, clusterDirectory, ...
                                                                    startFrame, endFrame, increment, validFrames);

% Go to the save directory
cd (saveDirectory); 

% Save M matrix 
if exist ('M', 'var')
   save('M.mat', 'M');
end

% Save MPM matrix
if exist ('MPM', 'var')
   MPM=MPM./clustBinSize; %JvR and Johan: to resize according to binning
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

% Save cluster properties
if exist ('frameProps', 'var')
   save ('frameProps.mat', 'frameProps');
end

% Save validFrames array
if exist ('validFrames', 'var')
   save ('validFrames.mat', 'validFrames');
end

% Tell the user we've finished
fprintf (1, '\nptTrackCells: Finished analysis of job number %d.\n', jobNumber);
