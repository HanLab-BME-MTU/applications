% This is the sequence of files that should be performed.
% The file settings.ini has to be placed in the same folder as the movie
% files. The movie must have the format file0001.tif, file0002.tif
% file0003.tif etc.

% An example of a good settings.ini file is as follows:

% [Image Parameters]
% fileName=Five
% frames=150
% thresholdFactor=1.9
% pixelsDilate=2
% skelFill=1
% diskSize=40
% intensityBased=1
% closeSize=6
% 
% 
% [Actions]
% saveTracks=1
% plotSkeletonizationSteps=0
% saveMovie=1
% 
% [Tracking parameters]
% maxDisplacement=30
% minLenght=30
% lostFrames=5
% dimensions=2
% 
% [Selected Phyllopodia]
% idAndOrigin=[11 146 108; 5 146 108; 9 146 70; 10 148 78]


%% Test Single frames

% First evaluate the skeletonization of single frames to make sure the
% basic shape is correct and the endpoints are accurate.
% Settings.ini can be opened with this same editor
% The frame number and the image folder are written in
% singleImageSkeletonAnalysis.m

singleImageSkeletonAnalysis

%% Create all skeletons

% Done. Now we need to create skeletons of all frames and locate endpoints
% and joints.
% analyzeSkeletonMovie.m is the file
% Remember to set in settings.ini the option plotSkeletonizationSteps to 0
% and the number of frames the movie has.

analyzeSkeletonMovie

%% Track Endpoints

% Now we have to track the endpoints and the maximal displacement and
% minimum number of tracked frames is important.
% The results will be placed in the folder Tracks and the movie will show
% the trajectories of the tracked endpoints
% Now we'll put numbers to how much the phyllopodia moved and a prliminary
% spreadsheet will be generateed.

trackSkeletons

%% Create segments

% This routine creates a movie with segments from all tracked endopints to
% the center of mass. It may be become useless, but now is the one that
% placed id numbers in the movie to select which ones to track.

makeSegmentsMovie

%% Create Segments to selected points

% Now, instead of having a common origin in the center off mass and all the
% endpoints to be tracked, we'll arbitrarily say which ones.
% The last line of Settings.ini needs tha info. We have to provide the ID
% of the endpoint and an approximate location of the origin. The routine
% will find the closest joint to that approximate location.
% The format MUST be [ID1 X1 Y1; ID2 X1 Y2; etc] as many as you want.
% The analysis includes a new spreadsheet with the lenght of the segment
% from the origin to the end of each selected endpoint at every frame (ArbitrarySegments.dat).

makeSegmentsWithArbitraryPoints

%% Follow the phyllopodia

% Finally, instead of drawing a straight line to the origin, it routine
% will follow the pillopodia. Uf, this one was hard!
% The new length in a spreadsheet is still missing, but I will do it.

followSkeletonAnalysis
