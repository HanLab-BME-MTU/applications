%% Initialization
mainPath = fullfile(getenv('HOME'),'Desktop','052710_cre_CSUX_2');
chDirs= dir(fullfile(mainPath,'ch*'));
nChan = numel(chDirs);

% Create channels
C(nChan,1)=Channel();
for i=1:nChan
    C(i)=Channel(fullfile(mainPath,chDirs(i).name,'roi'));
    C(i).emissionWavelength_=str2double(chDirs(i).name(3:end));
end

% Create movie
MD = MovieData(C,mainPath);
MD.setPath(mainPath);
MD.setFilename('movieData.mat');
MD.sanityCheck

%%
pathForImag  = MD.channels_(1).channelPath_;
sigmaPSF     = 1;
KSigma       = 4;
searchRadius = 5;
nFrames      = MD.nFrames_;
imSize       = MD.imSize_;
pixelSize    = 100;
MD.pixelSize_ = pixelSize;
%All the output files for each function will be stored here
outputPath   = fullfile(mainPath,'FADetection');
if ~isdir(outputPath), mkdir(outputPath); end

%% Package initialization
MD.reset
MD.addPackage(FocalAdhesionPackage(MD));

%% Segmentation
MD.packages_{1}.createDefaultProcess(1);
MD.packages_{1}.createDefaultProcess(2);
cellfun(@(x) x.run(struct('ChannelIndex',1)),MD.packages_{1}.processes_(1:2));

pathForMask  = MD.packages_{1}.processes_{2}.outFilePaths_{1,1};

%% Particle detection
% Use process anisotropic Gaussian detection
MD.packages_{1}.createDefaultProcess(3);
MD.packages_{1}.processes_{3}.run(struct('ChannelIndex',1,'MaskChannelIndex',1));

% Use Sylvain's function
featuresInfo = getMovieParticleDetection1(pathForImag, pathForMask, sigmaPSF, 'outputPath',outputPath);

movieInfo = MD.packages_{1}.processes_{3}.loadChannelOutput(1);

%% Particle tracking
% Run tracking
MD.packages_{1}.createDefaultProcess(4);
MD.packages_{1}.processes_{4}.run(struct('ChannelIndex',1));

% Use Sylvain's function
tracksFinal = getMovieParticleTracking1(featuresInfo,searchRadius, 'outputPath',outputPath);

tracks = MD.packages_{1}.processes_{4}.loadChannelOutput(1);

%% Track grouping
% Use process anisotropic Gaussian detection
MD.packages_{1}.createDefaultProcess(5);
MD.packages_{1}.processes_{5}.run(struct('ChannelIndex',1));

pathForMask =MD.processes_{2}.outFilePaths_{1};
distTransPath = [outputPath filesep 'distanceTransform'];
if ~isdir(distTransPath), mkdir(distTransPath); end
getMovieDistanceTransform1(pathForMask, distTransPath, nFrames)
  
%
getMoviePairTracks1(movieInfo, tracks, sigmaPSF, KSigma, nFrames,imSize,pixelSize, distTransPath, outputPath)
