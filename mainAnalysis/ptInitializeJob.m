function [newCoord] = ptInitializeJob (ptJob, jobNumber)
% ptInitializeJob finds coordinates in an image
%
% SYNOPSIS       ptInitializeJob (ptJob, jobNumber)
%
% INPUT          ptJob : a structure which contains the data for the job to process. See below
%                        for the exact structure details
%                jobNumber : which job is currently being dealt with. This is used to print
%                            some status information on the matlab command line
%                
% OUTPUT         newCoord : list of coordinates manually checked and edited by user.
%
% DEPENDENCIES   ptInitializeJob uses { imClusterSeg
%				        ptTrackLinker
%				        ptCheckMinimalCellDistance
%				        ptFindNuclei
%				        ptUserCellProcessing
%				        ptFindHalos }
%                                  
%                ptInitializeJob is used by { PolyTrack }
%	
% REMARK         the ptJobs structure looks as follows:
%                       imagedirectory : where are the images located
%                       imagename : what are the images called (sort of a template)
%                       imagenameslist : list of images within imagedirectory with imagename
%                       firstimage : which images shall we start with (refers to imagenameslist)
%                       lastimage : which images shall be the last one (refers to imagenameslist)
%                       intensitymax : highest value image can have (calc from bitdepth)
%                       fi_nucleus / la_nucleus : nucloi intensity first / last image
%                       fi_background / la_background : background intensity first / last image
%                       fi_halolevel / la_halolevel : halo intensity first image
%                       leveladjust : factor to adjust intensity difference nucloi/ background
%                       minsize : minimal size of nucloi
%                       maxsize : maximal size of nucloi
%                       minsdist : minimal distance between two cells
%                       minmaxthresh : onoff - should minima and segmentation be used
%                       clustering : on off - should clustering be used
%                       intensityMax : max intensity of the image
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source, make function independent of GUI handle

% First tell the user we're busy initializing
fprintf (1, 'Initializing phase has started for job number %d...\n', jobNumber);

% Assign all the job data to local variables
imageDirectory    = ptJob.imagedirectory;
imageNameList     = ptJob.imagenameslist;
firstImageNr      = ptJob.firstimage;
lastImageNr       = ptJob.lastimage;
levNucFirst       = ptJob.fi_nucleus;
levBackFirst      = ptJob.fi_background;
levHaloFirst      = ptJob.fi_halolevel;
levelAdjust       = ptJob.leveladjust;

% Minimal/maximal size of the nuclei
minSizeNuc        = ptJob.minsize;
maxSizeNuc        = ptJob.maxsize;

% An educated guess of the minimal distance between neighbouring cells
% (better to big than to small, for this value is used for the static
% search. Searches with templates are more tolerant.)
MinDistCellCell   = ptJob.minsdist;

% Get the maximum image intensity of the images in the job
maxImageIntensity = ptJob.intensityMax;

% User selectable clustering method
if ptJob.clustering
   method = 1;
elseif ptJob.minmaxthresh
   method = 2;
elseif ptJob.emclustering
   method = 3;
else
   fprintf(1, 'ptInitializeJob: Wrong method. A valid method (1,2,3) has to be selected.\n');
   return;
end

% Get the difference between nuclei and background intensity and adjust this
levDiffFirst = abs (levNucFirst - levBackFirst) * levelAdjust;

% How much the blobs found in ptFindHalos shall be eroded. This is an indirect
% size criteria for ptFindHalos. Increase - minimal size of halos will be
% increased, decrease - ... decreased
ErodeDiskSize = round ((sqrt (minSizeNuc)) / 2) ;

% Go to the directory where the images are stored
cd (imageDirectory)

% Make sure the images are in correct sequence. If not throw an error
if ~lastImageNr > firstImageNr
   fprintf (1, 'ptInitializeJob: Error: the last image # is before the first image # in job number %d\n', jobNumber);
   return;
else
   % Get the name of the first image in the list
   imageName = char (imageNameList (firstImageNr));

   % Read that image from disk using the given min and max intensity values and do
   % background subtraction as well
   tempFirstImage = imreadnd2 (imageName, 0, maxImageIntensity);
   [firstImage, background] = ptGetProcessedImage (tempFirstImage, 20);
   clear tempFirstImage;

   % Get the size of the image
   [img_h, img_w] = size (firstImage);

   % Determine the intensity level of the halos
   HaloLevel = (levHaloFirst - levBackFirst) * 2 / 3 + levBackFirst;

   % Depending on the method selected by the user do the appropriate segmentation
   if method == 1
      % Segment the image using the k-means method
      [segmentedImage, dummy, mu0] = imClusterSeg (firstImage, 0, 'method', 'kmeans', 'k_cluster', 3, 'mu0', ...
                                                   [levNucFirst ; levBackFirst ; levHaloFirst]);
      %[segmentedImage, dummy, mu0] = imClusterSeg (firstImage, 0, 'method', 'kmeans', 'k_cluster', 3);
      % Copy the image for later use
      inputImage = segmentedImage;

   elseif method == 2 
      % Specify only input image since the image will be segmented in the following functions
      inputImage = firstImage;
   elseif method == 3
      % Segment the image using the EM method
      [segmentedImage, dummy, mu0] = imClusterSeg (firstImage, 0, 'method', 'em', 'k_min', 2, 'k_max', 2);

      inputImage = firstImage;
   end

   % Find areas that are really dark and match cells into them
   [coordNuc, nucImage] = ptFindNuclei (inputImage, levDiffFirst, minSizeNuc, maxSizeNuc, method);

   % Find areas that look like round, big spots of pure light. We do this because sometimes
   % the pictures are of poor quality and display huge halos around certain cells
   [coordHalo, haloImage] = ptFindHalos (inputImage, ErodeDiskSize, HaloLevel, method);

   % Ensure a minimal distance between two found cells
   newCoord =  ptCheckMinimalCellDistance (coordNuc, coordHalo, MinDistCellCell);           

   % Let the user manually fill in the missing cells and erase the wrong cells
   newCoord = ptUserCellProcessing (firstImage, newCoord);

   % Round coordinates and return to the calling function
   newCoord = round (newCoord);

   % Tell the user we've finished
   fprintf (1, 'ptInitializeJob: Finished Initializing phase for job number %d.\n', jobNumber);

   return;
end

