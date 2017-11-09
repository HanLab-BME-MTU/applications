% Saves a movie with a composition of the algorithm steps for presentations

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


clear

imagesFolder='Images\Chris\052407 R18 movies\R18 RFP\21_4 hrs pi\';

%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName=getFromINI(settingsRead, 'fileName');
frames=getFromINI(settingsRead, 'frames');

warning off
mkdir([imagesFolder '\Presentation']);

if exist([imagesFolder 'Presentation\' fileName '.TIF']) == 2
        delete([imagesFolder 'Presentation\' fileName '.TIF'])
end

%% Make movie images
for itFrame=1:frames
    thisSkeleton=imread([imagesFolder 'Skeletons\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    thisTracks=imread([imagesFolder 'Tracks\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    thisSegments=imread([imagesFolder 'Segments\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    thisArbitrary=imread([imagesFolder 'Arbitrary\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    thisFollowSkeleton=imread([imagesFolder 'FollowSkeleton\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    thisNakedSkeleton=imread([imagesFolder 'nakedSkeletons\' fileName num2str(itFrame-1, '%04.f') '.tif']);
    thisNakedSkeletonRGB=cat(3,thisNakedSkeleton,thisNakedSkeleton,thisNakedSkeleton)*255;
    
    thisFrameComposed=[thisSkeleton thisTracks thisSegments;...
        thisArbitrary thisFollowSkeleton thisNakedSkeletonRGB];
    
    %imwrite(thisFrameComposed, [imagesFolder 'Presentation\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'JPG')
    imwrite(thisFrameComposed,[imagesFolder 'Presentation\' fileName '.TIF'],...
        'TIF','Compression','none','WriteMode','append','ColorSpace','rgb')
end
    