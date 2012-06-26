load /home/mv89/matlab/my_stuff/Data/2012_data/Vinay/2012_02_22/Cell2/Analysis/movieData.mat

iChannel     = 2;
pathForImag  = MD.channels_(iChannel).channelPath_;
pathForMask  = MD.processes_{5}.outFilePaths_{iChannel};
sigmaPSF     = 1;
KSigma       = 4;
searchRadius = 5;
nFrames      = MD.nFrames_;
imSize       = MD.imSize_;
pixelSize    = 100;
%All the output files for each function will be stored here
outputPath   = '/home/mv89/matlab/my_stuff/Data/2012_data/Vinay/2012_02_22/Cell2/Analysis';
% featuresInfo.mat - 
% tracks.mat       - 
% 

%'mode', 'xyArtc' 'xyarstc' (r = sx, s = sy)
%  Symbols: xp : x-position
%              yp : y-position
%               A : amplitude
%              sx : standard deviation along the rotated x axis
%              sy : standard deviation along the rotated y axis
%               t : angle
%               c : background
%'alpha', 0.05 - some test
%'kSigma', 4 - cutoff in number of standard deviations
%'minDist', .25 - minimum distance between 2 detected objects
%'filterSigma',psfSigma*sqrt(2)


%This function detects the FA's
featuresInfo = getMovieParticleDetection1(pathForImag, pathForMask, sigmaPSF, 'outputPath',outputPath);
% This function will write featuresInfo.mat


%This function tracks the detected FA's 

[tracksFinal] = getMovieParticleTracking1(featuresInfo,searchRadius, 'outputPath',outputPath);
%Include the outputPath - used to write the tracks.mat file


%This function will write a distanceTransform_iFrame.mat for each frame
% Better give it a separate folder as input 
distTransPath = [outputPath filesep 'distanceTransform'];
getMovieDistanceTransform1(pathForMask, distTransPath, nFrames)
 
 
 %
 getMoviePairTracks1(featuresInfo, tracksFinal, sigmaPSF, KSigma, nFrames,imSize,pixelSize, distTransPath, outputPath)