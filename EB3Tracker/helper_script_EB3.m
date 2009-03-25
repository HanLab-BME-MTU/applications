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
% automatically and stored in the roi_x/subROIs/sub_x folders.  after
% running it, you can run getProj again to add these sub-projects into the
% list for batch processing, or you can run metaEB3analysis at the
% individual project level.

subRoiConfig



%% VISUALIZATION OF THE RESULTS

% the function featVelMovie makes movies of moving dots corresponding to tracks. 
% the color of the dots corresponds to speed (where cool colors represent
% shrinkage speeds, and warm colors represent growth and/or pause phase)

% projData is the output from the metaEB3analysis function, but you can give [] to select the file
timeRange = [1 30];   % frame range [start end], use [] to use all frames
velLimit  = [];       % max speed to use for color min/max (i.e. it will plot all tracks faster than velLimit as same shade of red), if [], will use full range
roiYX     = [];       % [] to choose a ROI on the fly; or load one by dragging it into the workspace
featVelMovie(projData,timeRange,velLimit,roiYX);


% using trackMovie, you can make movies of one or more individual tracks,
% or of all the tracks within the timeRange and within the input ROI. the
% movies will be stored in whatever folder you select when prompted, along
% with the ROI coordinates used to make that movie. the movie name is
% chosen automatically as either:
%
% "allTracks_startFrame_endFrame_01" (or 02, 03, etc, depending on how many
% allTracks movies have been produced in the past)
%
% -OR-
%
% "track_trackNumber_startFrame_endFrame_01" (or 02, 03, etc, for
% individual movies)
%
% movies are never overwritten so you can make variations on a theme, for
% example: you might want to make 3 movies using the same ROI and timeRange
%           - one where you just show the raw images
%           - one where you show just the tracks
%           - one where you show the tracks with the feature centroids plotted
%             according to time

indivTrack = []; % track numbers if you want to make n movies of individual tracks. use [] for all tracks
timeRange  = [1 10]; % frame range [start end], use [] to use all frames
roiYX      = []; % [] to choose a ROI on the fly; or load one by dragging it into the workspace
magCoef    = []; % use [] to make largest possible movie based on screen size, otherwise try 3 or 4.
showTracks = 1; % 1 to show tracks, 0 if not
showDetect = 1; % 1 for feature positions on tracks only, 2 for all detected features, 3 for current frame (see function header)
trackMovie([],indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect)


% with plotTracks2D_EB3 you can overlay the tracks onto a single image and select them
% selectedTracks are not saved anywhere but will appear in the work space

timeRange       = [1 10]; % frame range over which to plot
img             = [];     % [] to select the image
ask4sel         = 'y'     % 'y' to select tracks, 'n' otherwise
plotCurrentOnly = []      % [] for all tracks in the frame range, a number in timeRange if you only want to show the tracks active in that frame
roiYX           = []      % [] to use whole image, or load one by dragging it into the workspace
[selectedTracks] = plotTracks2D_EB3([],timeRange,img,ask4sel,plotCurrentOnly,roiYX,[]);




