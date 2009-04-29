% This is a help page for analyzing EB3 comets.
% 2009 Kathryn Applegate

% HINT #1: for any function, you can type "edit functionName" at the command
% line to see the comments at the top.  these (hopefully!) provide more
% info about what the input parameters should be and what the output should
% look like.

% HINT #2: for each function below, change the relevant input parameters 
% (listed above the function). then copy the lines and function call as a 
% block and hit the "F9" key to execute.  alternatively, you can copy and
% paste onto the command line.


%% PRE-PROCESSING STEPS 

% the first time you begin work on a given project (i.e. movie), you need
% to move the raw images into a subdirectory under the parent, called
% "images".  if you have a bunch of movies for which you want to do this,
% just call the following function:

makeImageDirectory;


% next you should set up the "roi_x" directory where all analysis will be
% stored.  each roi_x directory will appear at the same level as "images".
% you can choose several (up to 9) per project.  if you don't want to spend
% the time to choose a region around the cell, you will be asked to select
% a background point.  choose a point in the center of a cell, like near
% the centrosome (if one is evident).  the detection function will use the
% ROI you select or the background point to determine some internal
% parameters.

selectROI     = 1; % 1 to select ROI(s), 0 to select background point (shortcut)
overwriteROIs = 1; % 1 if you want to overwrite previous ROIs
setupRoiDirectories(selectROI,overwriteROIs,[]);



%% BATCH MOVIE PROCESSING

% use getProj to set up a batch job. this function collects all the
% projects under a user-selected top-level directory.

queryString=[]; % enter a string like 'rac' or 'control' to find only projects which have this in the file path
[projList,projPathParsed]=getProj(queryString);


% to detect, track, and post-process for multiple movies, open the batch
% function and check the parameters.  then run the function by typing the
% function name at the command line.

edit batchEB3analysis


%% SINGLE MOVIE PROCESSING

% DETECTION AND TRACKING

% detection of features - results stored in roi_x/feat directory

timeRange = []; % frame range [start end], use [] to use all frames
bitDepth  = 12; % change according to your camera (should be 12, 14, or 16)
savePlots = 1;  % 1 to save overlay plots of detection results in subfolder; 0 if not (may run faster)
[movieInfo]=eb1SpotDetector([],timeRange,bitDepth,savePlots);


% tracking of features - open the script and check parameters there. save
% any changes.
edit scriptTrack_EB3.m

% now run the tracking.  results will be stored in roi_x/track directory
scriptTrack_EB3


% POST-PROCESSING

% the metaEB3analysis function extracts some useful information from the
% tracks - results stored in roi_x/meta directory. open the function and
% read the header to see what the various fields represent...we can add
% more functionality here for sure.
% (i.e. k: 2.0s,105nm;  y: 2.0s,84nm;  c: 0.8s,110nm)

secPerFrame = .8; % frame rate (seconds)
pixSizeNm   = 110; % real-space pixel size (nanometers)
[projData]=metaEB3analysis([],secPerFrame,pixSizeNm);


% use popHist to make histograms of growth/gap/shrinkage speeds and the
% distribution representing how much the speed changes between any two
% frames. these are saved both as matlab figures (can be opened from within
% matlab) and as tifs (poor quality but easy to view), in the meta folder

popHist(projData)


% once you have detected/tracked a region of interest (like a cell), you
% may want to choose sub-regions to further analyze. use subRoiConfig to
% choose up to 9 sub-regions.  the data from tracking will be pulled
% automatically and stored in the roi_x/subROIs/sub_x folders.

subRoiConfig

% after running subRoiConfig, if you want to see the tracks, go to the
% sub_x folder and load subIdx and the roiYX then go to the roi_x folder
% which you're using and load projData.  then run the following, and choose
% an image to use for the overlay according to the prompt.

timeRange=[];
img=[];
ask4sel=[];
[selectedTracks] = plusTipPlotTracks(projData,subIdx,timeRange,img,ask4sel,[],roiYX);

% you can access the data this way:
subTrackProfiles=projData.nTrack_start_end_velMicPerMin_class_lifetime(subIdx,:);
% and here we can make a histogram of growth speeds, for example...
growthSpeeds=subTrackProfiles(subTrackProfiles(:,5)==1,4);
figure; hist(growthSpeeds,25)



%% VISUALIZATION OF THE RESULTS

% use the GUI
plusTipTrackViz

plusTipPlotTracks

