% Runs the skeletonization of all frames and saves all locations of ends
% and joints
% The file Settings.ini has to be placed in the same foder as the movie

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');
fileName=getFromINI(settingsRead, 'fileName');
frames=getFromINI(settingsRead, 'frames');
thresholdFactor=getFromINI(settingsRead, 'thresholdFactor');
pixelsDilate=getFromINI(settingsRead, 'pixelsDilate');
skelFill=getFromINI(settingsRead, 'skelFill');
closeSize=getFromINI(settingsRead, 'closeSize');
diskSize=getFromINI(settingsRead, 'diskSize');
plotSkeletonizationSteps=getFromINI(settingsRead, 'plotSkeletonizationSteps');
intensityBased=getFromINI(settingsRead, 'intensityBased');
diskSize=getFromINI(settingsRead, 'diskSize');
EdgeStrongness = getFromINI(settingsRead, 'EdgeStrongness'); 
saveMovie=getFromINI(settingsRead, 'saveMovie');
saveTracks=getFromINI(settingsRead, 'saveTracks');

skeletonEnds=[];
skeletonJoints=[];
centroids=[];

%% Creates the folder were the folowing analyses will place the files

warning off
mkdir([imagesFolder '\Skeletons']);
mkdir([imagesFolder '\nakedSkeletons']);
mkdir([imagesFolder '\Tracks']);
mkdir([imagesFolder '\Segments']);
mkdir([imagesFolder '\Arbitrary']);
mkdir([imagesFolder '\FollowSkeleton']);


%% Read movie and skeletonize
itwaitbar = waitbar(0/frames,['Analyzing frame ',num2str(0),' of ',num2str(frames)]);
for itFrame=1:frames
    thisFrame = imread([imagesFolder fileName num2str(itFrame-1, '%04.f') '.TIF']);
    [thisFrameSkeleton, thisCentroid] = getSkeleton(thisFrame, thresholdFactor, ...
    pixelsDilate, plotSkeletonizationSteps, skelFill, diskSize, closeSize, intensityBased,EdgeStrongness);
    
    thisSkeletonEnds   = find_skel_ends(thisFrameSkeleton);
    thisSkeletonJoints = find_skel_intersection(thisFrameSkeleton);
    thisCentroid       = [thisCentroid(1, 2) thisCentroid(1, 1)];
    
    %% Saves the frames with their skeletons
    if saveMovie == 1
        saveSkeletonImage(thisFrame, thisFrameSkeleton, thisSkeletonEnds,...
            thisSkeletonJoints, thisCentroid, [fileName num2str(itFrame-1, '%04.f')] , imagesFolder);
        imwrite(thisFrameSkeleton, [imagesFolder 'nakedSkeletons\' fileName num2str(itFrame-1, '%04.f') '.tif'], 'TIF')
    end    
    thisSkeletonEnds   = [thisSkeletonEnds, ones(size(thisSkeletonEnds, 1), 1)*itFrame];
    thisSkeletonJoints = [thisSkeletonJoints, ones(size(thisSkeletonJoints, 1), 1)*itFrame];
    thisCentroid       = [thisCentroid itFrame];
    skeletonEnds       = [skeletonEnds; thisSkeletonEnds];
    skeletonJoints     = [skeletonJoints; thisSkeletonJoints];
    centroids          = [centroids; thisCentroid];
    itFrame
    itwaitbar = waitbar(itFrame/frames,itwaitbar,['Analyzing frame ',num2str(itFrame),' of ',num2str(frames)]);
end

%% Saves locations

itwaitbar = waitbar(.5,itwaitbar,['Saving Tracks']);
if saveTracks==1
    save([imagesFolder 'Tracks\' fileName 'Ends.dat'], 'skeletonEnds', '-ASCII', '-DOUBLE');
    save([imagesFolder 'Tracks\' fileName 'Joints.dat'], 'skeletonJoints', '-ASCII', '-DOUBLE');
    save([imagesFolder 'Tracks\' fileName 'CenterOfMass.dat'], 'centroids', '-ASCII', '-DOUBLE');
end
itwaitbar = waitbar(.5,itwaitbar,['Saving log of the settings']);   
%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'Skeletons\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');
close(itwaitbar)