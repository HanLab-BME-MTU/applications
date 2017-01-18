% testing a typical L1 force reconstruction

% forceReconstructionL1UsingTFMPackage runs TFM package all the way to
% force reconsturction without any stopping using L1 reconstruction and
% using L-curve.
% previously: forceReconstructionL1UsingTFMPackage(beadFolder,cellFolder,refImgPath,NA,pixelSize,timeInterval,refROI)


% Analysis Output Directory
analysis_dir = '/work/bioinformatics/s170480/Analysis/TFM';
% ROI_dir = '/work/bioinformatics/s170480/Analysis/TFM/ROI'

% Raw Input Image paths
% (pulled from: )
beadFolder = '/work/bioinformatics/s170480/Data/TFM/Real/Cell5/Beads';
cellFolder = '/work/bioinformatics/s170480/Data/TFM/Real/Cell5/Cell';
refImgPath = '/work/bioinformatics/s170480/Data/TFM/Real/Cell5/Avg_Ref.tif';
roiMask = '/work/bioinformatics/s170480/Data/TFM/Real/Cell5/roiMask2.tif';

% Imaging/Microscope Parameters
NA = 1.49;
pixelSize = 60;
timeInterval = 1;

% (Some) Algorithm Parama
solMethodBEM='1NormReg';
useLcurve = true;
YoungModulus = 16000; %3500;
regParam = 0.0036;


% Get refROI mask with:
% sampleCellImg = imread('Cell_11_w1561 TMR_t001.tif');
% figure, imshow(sampleCellImg,[]), hold on
% h=imrect;
% refROI = wait(h);

% Also make sure you make ROI tif using manualRectangle

% Retrieve current location
% Channel creation
% Create a channels object
channel = Channel(beadFolder);
channel.fluorophore_='alexa647';
channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
channel.imageType_='TIRF';

channel(2) = Channel(cellFolder);
channel(2).fluorophore_='egfp';
channel(2).emissionWavelength_=name2wavelength('egfp')*1e9;
channel(2).imageType_='TIRF';


% MovieData creation
% Constructor needs an array of channels and an output directory (for analysis)
outputDir = analysis_dir; %fileparts(analysis_dir);
MD = MovieData(channel, outputDir);

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

% Load the movie
clear MD;
MD = MovieData.load(fullfile(outputDir,'movieData.mat'));

roiDir = [outputDir filesep 'ROI'];
if ~exist(roiDir,'dir')
    mkdir(roiDir);
end

% use manualRenctangle for making roi mask, rename it to roiMask.tif and
% place it to roiDir.
% Apply ROI
roiMD = MD.addROI(roiMask, roiDir);
roiMD.setPath(roiDir);
roiMD.setFilename('movieData.mat');
roiMD.save;
MD.save;

% Create TFM package and retrieve package index
roiMD.addPackage(TFMPackage(roiMD));
iPack =  roiMD.getPackageIndex('TFMPackage');


%% Create first process
disp('===================================================================');
disp('Running (1st) TFM Process');
disp('===================================================================');
roiMD.getPackage(iPack).createDefaultProcess(1)
params = roiMD.getPackage(iPack).getProcess(1).funParams_;
params.referenceFramePath = refImgPath;
% Get refROI mask with:
% sampleCellImg = imread('Cell_11_w1561 TMR_t001.tif');
% figure, imshow(sampleCellImg,[]), hold on
% h=imrect;
% refROI = wait(h);
params.cropROI = [1249, 1159, 194, 162]; %; refROI;
params.ChannelIndex=[1 2];
params.minCorLength=21;
params.maxFlowSpeed=5;
params.doPreReg=true;
roiMD.getPackage(iPack).getProcess(1).setPara(params);
% Run the stage drift correction
roiMD.getPackage(iPack).getProcess(1).run();


%% Create second process
disp('===================================================================');
disp('Running (2nd) TFM Process');
disp('===================================================================');
roiMD.getPackage(iPack).createDefaultProcess(2)
params = roiMD.getPackage(iPack).getProcess(2).funParams_;
% Parameters in displacement field tracking
params.referenceFramePath = refImgPath;
params.alpha = 0.05;
params.minCorLength = 19;
params.maxFlowSpeed = 40;
params.highRes = true;
params.mode = 'accurate';
roiMD.getPackage(iPack).getProcess(2).setPara(params);
roiMD.save;
% Run the displacement field tracking
roiMD.getPackage(iPack).getProcess(2).run();


%% Create third process and run
disp('===================================================================');
disp('Running (3rd) TFM Process');
disp('===================================================================');
roiMD.getPackage(iPack).createDefaultProcess(3)
params = roiMD.getPackage(iPack).getProcess(3).funParams_;
roiMD.getPackage(iPack).getProcess(3).setPara(params);
roiMD.getPackage(iPack).getProcess(3).run();


%% Create force reconstruction process and run
disp('===================================================================');
disp('Running (4th) TFM Process');
disp('===================================================================');
roiMD.getPackage(iPack).createDefaultProcess(4)
params = roiMD.getPackage(iPack).getProcess(4).funParams_;
params.YoungModulus = YoungModulus;
params.regParam = regParam;
params.method = 'FastBEM';
params.solMethodBEM = solMethodBEM; %'1NormRegLaplacian';
params.useLcurve = useLcurve;
params.basisClassTblPath = '/work/bioinformatics/s170480/Data/TFM/BasisFunctions/basisClass16kPa12pix.mat';
% Run
roiMD.getPackage(iPack).getProcess(4).setPara(params);
roiMD.getPackage(iPack).getProcess(4).run();
roiMD.save;