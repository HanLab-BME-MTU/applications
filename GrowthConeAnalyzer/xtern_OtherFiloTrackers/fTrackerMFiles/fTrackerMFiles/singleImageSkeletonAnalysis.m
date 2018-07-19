% Serves as a test tool in order to determine the good settings before
% runing the full analysis. The frame to be tested and the folder where the
% files are have to be properly set.
% Now that we have the gui, it became useless.
% The parameters are varied in the settings.ini file, which should be
% places in the same folder as the movie.

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName = getFromINI(settingsRead, 'fileName');
frames   = getFromINI(settingsRead, 'frames');
thresholdFactor=getFromINI(settingsRead, 'thresholdFactor');
pixelsDilate=getFromINI(settingsRead, 'pixelsDilate');
skelFill=getFromINI(settingsRead, 'skelFill');
diskSize=getFromINI(settingsRead, 'diskSize');
intensityBased=getFromINI(settingsRead, 'intensityBased');
closeSize=getFromINI(settingsRead, 'closeSize');
EdgeStrongness = getFromINI(settingsRead, 'EdgeStrongness'); 

plotSkeletonizationSteps = 0;

skeletonEnds   = [];
skeletonJoints = [];

%% Read frame and skeletonize
intwarning('off')
thisFrame = uint8(imread([imagesFolder fileName num2str(frameToTest, '%04.f') '.tif']));
[thisFrameSkeleton, thisCentroid,singleCone,closedCone,imageBackground,dilatedCone,OtherImage] = getSkeleton(thisFrame, thresholdFactor, ...
    pixelsDilate, plotSkeletonizationSteps, skelFill, diskSize, closeSize, intensityBased,EdgeStrongness);

if intensityBased == 0
    edgesImage = OtherImage{1};
    filledCone = OtherImage{2};
else
    binaryImage = OtherImage{1};
    filteredCone = OtherImage{2};    
end

%% Find locations of joints and ends
thisSkeletonEnds   = find_skel_ends(thisFrameSkeleton);
thisSkeletonJoints = find_skel_intersection(thisFrameSkeleton);
   
%% Show an image with the skeleton superimposed with the joinst and ends as
%% colored circles

rgbFrame = showSkeletonImage(thisFrame, thisFrameSkeleton, thisSkeletonEnds, thisSkeletonJoints,plotSkeletonizationSteps);