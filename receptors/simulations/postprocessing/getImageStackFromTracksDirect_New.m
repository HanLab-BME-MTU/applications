function getImageStackFromTracksDirect_New(tracksSim,timeInfo,intensityInfo,spaceInfo,saveInfo)
% GETIMAGESTACKFROMTRACKSDIRECT get images from compTracks of the simulated
% data. The images have the format as the experimental data, mimic the
% fluorescence for each recepptor. The outputs are used to calculate tracks
% in the packageGUI function uTrackPackageGUI
%
%   SYNPOSIS: getImageStackFromTracksDirect(tracksSim,saveDir)
%
%   Input:
% 
%       tracksSim     : is the compTrack,tracks structure array in default
%       format.
% 
%       timeInfo      : Structure with fields:
%           
%           .simTimeStep : Time between frames/time points (s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%           .cutOffTime     : Total time of simulation (s) to be considered,
%                             it can be less than simTime.
%
%       intensityInfo : Structure with fields:
%         
%           .bgav           : average background level (can be either single
%                             value or, for background inequality, an image
%                             of the same size as imsize). In counts, e.g.
%                             assuming a 16-bit camera.
%           .bgnoise        : std of backgound noise. In counts, e.g.
%                             assuming a 16-bit camera.
%           .scaleFactor    : is the scale factor for the intensities  
%
%        spaceInfo    : Structure with fields:
%         
%           .pixelSize      : pixel size in microns
%           .psfSigma        : width of point-spread function (in pixels),
%                             where sigma = (1/3)*(radius of Airy disc)
%           .imsize         : size of image [sx,sy]
%          
%         saveInfo    : Structure with fields:
%         
%           .saveVar      : variable that indicates whether or not tif files
%                           are saved to file (1/0)       
%           .saveDir      : directory where the images will be saved.
%          
% 
%
% Modified September 2016, Luciana de Oliveira 
%% Input

%call imput variables from structure

% simulation tracks
 tracksSim=tracksSim.compTracks;

% time info
simTimeStep = timeInfo.simTimeStep;
sampleStep = timeInfo.sampleStep;
cutOffTime = timeInfo.cutOffTime;

%intensity info
bgav=intensityInfo.bgav;
bgnoise=intensityInfo.bgnoise;
scaleFactor=intensityInfo.scaleFactor;

%space info
pixelSize =spaceInfo.pixelSize;
psfSigma =spaceInfo.psfSigma;
imsize =spaceInfo.imsize;

%save info
saveVar =saveInfo.saveVar;
saveDir =saveInfo.saveDir;

%% Calculations

%convert tracks structure to matrix
trackedFeatureInfo = convStruct2MatIgnoreMS(tracksSim,0); 

%Use a portion of the simulation if desired (truncating) given by the value
% of endIter,i.e, only the first total endIter=cutOffTime*simTimeStep will
% be considered for the analysis

endIter=cutOffTime/simTimeStep;


trackedFeatureInfo = trackedFeatureInfo(:,1:endIter);

%Scale up intensities
trackedFeatureInfo(:,4:8:end) = trackedFeatureInfo(:,4:8:end) * scaleFactor;

%rescaling in function of pixel size
trackedFeatureInfo(:,1:8:end) = trackedFeatureInfo(:,1:8:end)/pixelSize;
trackedFeatureInfo(:,2:8:end) = trackedFeatureInfo(:,2:8:end)/pixelSize;

%remove 0 columns
trackedFeatureInfo(:,3:8:end) = [];
trackedFeatureInfo(:,4:7:end) = [];
trackedFeatureInfo(:,4:6:end) = [];
trackedFeatureInfo(:,4:5:end) = [];
trackedFeatureInfo(:,4:4:end) = [];

%Replace zeros with NaNs
trackedFeatureInfo(trackedFeatureInfo == 0) = NaN;

%Subsample whole simulation every convStep iterations
convStep = 3*round(sampleStep/simTimeStep);
xCoord = trackedFeatureInfo(:,1:convStep:end);
yCoord = trackedFeatureInfo(:,2:convStep:end);
amps = trackedFeatureInfo(:,3:convStep:end);
[numTracks,numIters] = size(xCoord);

% Alleviate boundary effects on detection: add 10 pixels*psfSigma to
% particle positions; a few lines down image size will be expanded as
% well
xCoord = xCoord + 10*psfSigma;
yCoord = yCoord + 10*psfSigma;

%reformat track information following function input
trackInfo = reshape([xCoord;yCoord;amps],numTracks,3*numIters);

%Convert image size from micron square to pixels also add 10*psfSigma pixels
%to each image size to alleviate boundary effects on detection

imsize(1:2) = round(imsize/pixelSize) + 20*psfSigma;

% rad is the radius used for the generation of the Airy disc
% for each object, in increments of psfSigma 

rad=3*psfSigma;

%call the function that will make the image and save it in saveDir
makeAiryImageFromMPM(trackInfo,bgav,bgnoise,psfSigma,imsize,rad,saveVar,[],saveDir);

end



