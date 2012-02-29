% Script to test the consistency between the old plusTipTracker interface
% and its integration into the movie management platform
% Assuming test images are located on the server, create a projData and a
% movie object using these images. Perform sucessfully detection, tracking
% and classification on these images using both interfaces and check
% results validity.
%
% Sebastien Besson, Feb 2012

% Location to the test images
mainPath = fullfile(getenv('HOME'),'files','LCCB','comet','plusTipTrackerTest');

frameRate=1;
pixelSize=66;


%% projData construction

projData.imDir = fullfile(mainPath,'images');
projData.anDir = fullfile(mainPath,'roi_1');

if ~exist(projData.anDir,'dir'),mkdir(projData.anDir); end
imFiles=imDir(projData.imDir);
I = imread(fullfile(projData.imDir,imFiles(1).name));

imwrite(true(size(I)),fullfile(projData.anDir,'roiMask.tif'));
roiYX = [1 size(I,1); 350 size(I,1); size(I,2) 350; size(I,2) 1; 400 1; 1 400; 1 size(I,1)];
roiMask = poly2mask(roiYX(:,1),roiYX(:,2),size(I,1),size(I,2));
save(fullfile(projData.anDir,'roiYX.mat'),'roiYX');
imwrite(roiMask,fullfile(projData.anDir,'roiMask.tif'));

%% Movie object creation
MD=MovieData(Channel(projData.imDir),mainPath);
MD.setPath(mainPath);
MD.setFilename('movieData.mat');
MD.sanityCheck;
MD.camBitdepth_=16;
MD.timeInterval_=frameRate;
MD.pixelSize_=pixelSize;

% Create Roi
roiPath = fullfile(mainPath,'movieROI');
if ~exist(roiPath,'dir'),mkdir(roiPath); end
roiMaskPath = fullfile(roiPath,'mask.tif');
imwrite(roiMask,roiMaskPath);
MD.addROI(roiMaskPath,roiPath);
movie = MD.rois_(1);
movie.setPath(roiPath);
movie.setFilename('movieData.mat');
movie.save;

% Add PlusTipTracker package
movie.addPackage(PlusTipTrackerPackage(movie));
procConstr= PlusTipTrackerPackage.getDefaultProcessConstructors;

%% Comet detection 
% 1 - plusTipTacker interface
plusTipCometDetector(projData,[],16,1)

% 2 - PlusTipTrackerPackage
if isempty(movie.packages_{1}.processes_{1})
    movie.packages_{1}.createDefaultProcess(1)
end   
detProc=movie.packages_{1}.processes_{1};
detProc.run()

% 3- Compare output
s=load(movie.processes_{1}.outFilePaths_{1});
s2=load(fullfile(projData.anDir,'feat','movieInfo.mat'));

assertEqual(s.movieInfo,s2.movieInfo)

%% Comet Tracking 
% 1- plusTipTracker
timeWindow=5;
minTrackLen=3;
minRadius=2;
maxRadius=10;
maxFAngle=30;
maxBAngle=10;
maxShrinkFactor=1.5;
fluctRad=1.0;
timeRange=[];
diagnostics=false;
plusTipCometTracker(projData,timeWindow,minTrackLen,minRadius,maxRadius,...
    maxFAngle,maxBAngle,maxShrinkFactor,fluctRad,timeRange,diagnostics)

% 2- PlusTipTracker package
if isempty(movie.packages_{1}.processes_{2})
    movie.packages_{1}.createDefaultProcess(2)
end   
trackProc=movie.packages_{1}.processes_{2};
trackProc.run();

% 3- Compare output
s=load(trackProc.outFilePaths_{1});
s2=load(fullfile(projData.anDir,'track','trackResults.mat'));

assertEqual(s.tracksFinal,s2.tracksFinal);

%% Post-processing
% 1- plusTipTracker
makeHist=true;
plusTipPostTracking(projData,frameRate,pixelSize,[],makeHist)

% 2- PlusTipTracker package
if isempty(movie.packages_{1}.processes_{3})
    movie.packages_{1}.createDefaultProcess(3)
end   
postProc=movie.packages_{1}.processes_{2};  
postProc.run();

% 3- Compare output
s=load(postProc.outFilePaths_{1});
s2=load(fullfile(projData.anDir,'meta','projData.mat'));

assertEqual(s.projData,s2.projData);