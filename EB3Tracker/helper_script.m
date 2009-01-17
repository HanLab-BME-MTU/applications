
% HINT #1: for any function, you can type "edit functionName" at the command
% line to see the comments at the top.  these (hopefully!) provide more
% info about what the input parameters should be and what the output should
% look like.

% HINT #2: for each function below, change the relevant input parameters 
% (listed above the function). then copy the lines and function call as a 
% block and hit the "F9" key to execute.  alternatively, you can copy and
% paste onto the command line.


%%%%%%% PRE-PROCESSING STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
overwriteROIs = 1; % 1 if you want to 
setupRoiDirectories(selectROI,overwriteROIs,[]);




%%%%%%% DETECTION AND TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% detection of features - results stored in roi_x/feat directory
timeRange = []; % [] for all frames; change to [1 100], for example, to detect features in frames 1-100
bitDepth  = 16; % change according to your camera (should be 12, 14, or 16)
savePlots = 1;  % 1 to save overlay plots of detection results in subfolder; 0 if not (may run faster)
[movieInfo]=eb1SpotDetector([],timeRange,bitDepth,savePlots);

% tracking of features - open the script and check parameters there. save
% any changes.
edit scriptTrack_EB3.m

% now run the tracking.  results will be stored in roi_x/track directory
scriptTrack_EB3




%%%%%%% POST-PROCESSING STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function extracts some useful information from the tracks - results
% stored in roi_x/meta directory. open the function and read the header to
% see what the various fields represent...we can add more functionality
% here for sure.
secPerFrame = 2;    % frame rate (seconds)
pixSizeNm   = 105;  % real-space pixel size (nanometers)
[projData]=metaEB3analysis([],secPerFrame,pixSizeNm);

% this function (formerly called popVelScatter) makes movies of moving dots
% corresponding to tracks. the color of the dots corresponds to speed
% (where cool colors represent shrinkage speeds, and warm colors represent 
% growth and/or pause phase)
projData  = projData; % use output from metaEB3analysis function, or give [] to select the file
timeRange = [];       % [] for all frames; change to [1 10], for example, to make movie over just frames 1-10
velLimit  = [];       % max speed to use for color min/max (i.e. it will plot all tracks faster than velLimit as same shade of red), if [], will use full range
roiYX     = [];       % [] to choose a ROI on the fly; or load one by dragging it into the workspace
featVelMovie(projData,timeRange,velLimit,roiYX);

% here we can make movies of one or more individual tracks, or of all the
% tracks within the timeRange and within the input ROI
% the movies will be stored in whatever folder you select when prompted,
% along with the ROI coordinates used to make that movie. the movie name is
% chosen automatically as either:
% allTracks_startFrame_endFrame_01 (or 02, 03...depending on how many allTracks movies have been produced in the past)
% OR
% track_trackNumber_startFrame_endFrame_01 (or 02, 03...etc)
% movies are never overwritten so you can make variations on a theme, for
% example: you might want to make 3 movies using the same ROI and timeRange
%           - one where you just show the raw images
%           - one where you show just the tracks
%           - one where you show the tracks with the feature centroids plotted
% according to time
indivTrack = []; % [] to show all tracks; [1 5 44] to make 3 movies of just tracks 1, 5, and 44, etc.
timeRange  = []; % [] for all frames; change to [1 10], for example, to make movie over just frames 1-10 (ignored if using indivTrack - then whole track is plotted)
roiYX      = []; % [] to choose a ROI on the fly; or load one by dragging it into the workspace
magCoef    = 3;  % [] for largest movie possible on your screen; something like 3 should be fine if smaller movie desired
showTracks = 1;  % 1 to show the tracks, 0 to not show them
showDetect = 1;  % 0 to not show time-color-coded detection, 1 to show detection just on tracks, 2 to show ALL detected features (including ones not selected in tracking step)
trackMovie([],indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect);

% here we can overlay the tracks onto an image and select them
% selectedTracks are not saved anywhere but will appear in the work space
timeRange       = [1 10]; % frame range over which to plot
img             = [];     % [] to select the image
ask4sel         = 'y'     % 'y' to select tracks, 'n' otherwise
plotCurrentOnly = []      % [] for all tracks in the frame range, a number in timeRange if you only want to show the tracks active in that frame
roiYX           = []      % [] to use whole image, or load one by dragging it into the workspace
[selectedTracks] = plotTracks2D_EB3([],timeRange,img,ask4sel,plotCurrentOnly,roiYX,[]);

% make histograms of seg/gap speeds and the distribution representing how
% much the speed changes between any two frames. these are saved both as
% matlab figures (can be opened from within matlab) and as tifs (poor
% quality but easy to view), in the meta
projData  = projData; % use output from metaEB3analysis function, or give [] to select the file
popHist(projData)










