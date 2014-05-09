% U2OSstack4.m is a script for Kris's 2014-04-09/DsRed Dual-cube 561stack3
% typical L1 force reconstruction
function [] = forceReconstructionL1UsingTFMPackage(beadFolder,cellFolder,refImgPath,NA,pixelSize,timeInterval,refROI)
% forceReconstructionL1UsingTFMPackage runs TFM package all the way to
% force reconsturction without any stopping using L1 reconstruction and
% using L-curve.

% Example
% beadFolder = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Kris
% Kubow/2014-04-09/DsRed Dual-cube 561stack2/Beads';
% cellFolder = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Kris
% Kubow/2014-04-09/DsRed Dual-cube 561stack2/Pax'
% refImgPath = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Kris Kubow/2014-04-09/DsRed Dual-cube 561stack2/Ref/DsRed Dual-cube 561stack2 trypsin.tif';
% NA = 1.4;
% pixelSize = 71;
% timeInterval = 3;
% Get refROI mask with:
% sampleCellImg = imread('Cell_11_w1561 TMR_t001.tif');
% figure, imshow(sampleCellImg,[]), hold on
% h=imrect;
% refROI = wait(h);
% Also make sure you make ROI tif using manualRectangle
% forceReconstructionL1UsingTFMPackage(beadFolder,cellFolder,refImgPath,NA,pixelSize,timeInterval,refROI)

%% Now force reconstruction via movieData (non-GUI mode)
% Retrieve current location
%% Channel creation
% Create a channels object
channel = Channel(beadFolder);
channel.fluorophore_='alexa647';
channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
channel.imageType_='TIRF';

channel(2) = Channel(cellFolder);
channel(2).fluorophore_='tmr';
channel(2).emissionWavelength_=name2wavelength('tmr')*1e9;
channel(2).imageType_='TIRF';
%% MovieData creation
% Constructor needs an array of channels and an output directory (for analysis)
outputDir = fileparts(beadFolder);
MD = MovieData(channel,outputDir);

% Set the path where to store the MovieData object.
MD.setPath(outputDir);
MD.setFilename('movieData.mat');

% Set some additional movie properties
MD.numAperture_= NA; %1.49;
MD.pixelSize_= pixelSize;%71;
MD.camBitdepth_=16;
MD.timeInterval_ = timeInterval;%5;
MD.notes_= 'Typical TFM analysis with L1 '; 

% Run sanityCheck on MovieData. 
% Check image size and number of frames are consistent. 
% Save the movie if successfull
MD.sanityCheck;
% Save the movie
MD.save;

%% Load the movie
clear MD
MD=MovieData.load(fullfile(outputDir,'movieData.mat'));

roiDir = [outputDir filesep 'ROI'];
if ~exist(roiDir,'dir')
    mkdir(roiDir)
end

% use manualRenctangle for making roi mask, rename it to roiMask.tif and
% place it to roiDir.
%% Apply ROI
roiMD = MD.addROI([roiDir filesep 'roiMask.tif'],roiDir);
roiMD.setFilename('movieData.mat');
roiMD.save;

%% Create TFM package and retrieve package index
roiMD.addPackage(TFMPackage(MD));
iPack=  roiMD.getPackageIndex('TFMPackage');

%% Create first process
roiMD.getPackage(iPack).createDefaultProcess(1)
params = roiMD.getPackage(iPack).getProcess(1).funParams_;
params.referenceFramePath = refImgPath;
% Get refROI mask with:
% sampleCellImg = imread('Cell_11_w1561 TMR_t001.tif');
% figure, imshow(sampleCellImg,[]), hold on
% h=imrect;
% refROI = wait(h);
params.cropROI = refROI;
params.ChannelIndex=1;
params.minCorLength=21;
params.maxFlowSpeed=5;
params.doPreReg=true;
roiMD.getPackage(iPack).getProcess(1).setPara(params);
%% Run the stage drift correction
roiMD.getPackage(iPack).getProcess(1).run();

%% Create second process
roiMD.getPackage(iPack).createDefaultProcess(2)
params = roiMD.getPackage(iPack).getProcess(2).funParams_;

%% Parameters in displacement field tracking

params.referenceFramePath = refImgPath;
params.maxFlowSpeed = 40;
params.alpha = 0.05;
params.minCorLength = 15;
params.maxFlowSpeed = 40;
params.highRes = true;
params.mode = 'accurate';
roiMD.getPackage(iPack).getProcess(2).setPara(params);
%% Run the displacement field tracking
roiMD.getPackage(iPack).getProcess(2).run();
%% Create third process and run
roiMD.getPackage(iPack).createDefaultProcess(3)
params = roiMD.getPackage(iPack).getProcess(3).funParams_;
roiMD.getPackage(iPack).getProcess(3).setPara(params);
roiMD.getPackage(iPack).getProcess(3).run();

%% Create force reconstruction process and run
roiMD.getPackage(iPack).createDefaultProcess(4)
params = roiMD.getPackage(iPack).getProcess(4).funParams_;

params.YoungModulus = 3500;
params.regParam = 1e-4;
params.method = 'FastBEM';
params.solMethodBEM = '1NormReg';
params.useLcurve = true;
params.basisClassTblPath = '/project/....mat';

params.solMethodBEM = solMethodBEM;
params.basisClassTblPath = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM basis functions/basisClass90x Kubow.mat';

roiMD.getPackage(iPack).getProcess(4).setPara(params);
roiMD.getPackage(iPack).getProcess(4).run();

roiMD.save;
end