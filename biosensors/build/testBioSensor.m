% saveFolder : where to save movieData file
% imgFolder : where your original image data is stored

%% Channel creation : you first always need to create the channels
% typically the biosensor data has at least two channels.
% a FRET channel (Acceptor) and a Donor channel
% the mCherry here is a third channel on which I can generally get better
% signal to noise and typically perform the segmentation.
channels = {'C1_mCherry','C2_Donor','C3_FRET'};
for iCh = 1:3
    % Create a channels object
    imgFolderC = [imgFolder filesep channels{iCh}];
    channel(iCh) = Channel(imgFolderC);
end
% Constructor needs an array of channels and an output directory (for analysis)
MD = MovieData(channel,saveFolder);


% Set the path where to store the MovieData object.
MD.setPath(saveFolder);
MD.setFilename('movieDataBiosensors.mat');

MD.sanityCheck;

% Save the movie
MD.save;
%% Threshold : Create a binary mask by thresholding
% This is actually the stage in movieData I replace
% as I use different functions for segmentation not in process format quite yet.
% This is the traditional way to go however. You will see for this movie
% how the binarization is crap using this method.

% Create the Threshold Process
thresProc = ThresholdProcess(MD);

% Save the process in the movie object
MD.addProcess(thresProc);

% Create a segmentation package
segPackage = SegmentationPackage(MD);

%  Save the package in the movie object
MD.addPackage(segPackage);

% Get the threshold parameters so you can modify them
params = MD.processes_{1}.funParams_;

%params.ChannelIndex = [1,2,3];

% Set the Thresholding Parameters
params.MethodIndx = 2;
% Here are the methods for the thresholdProcess
% MD.processes_{1}.getMethods
%  1 = 'thresholdFluorescenceImage'; (ie minMax)
%  2 = 'thresholdOtsu';
%  3 = 'thresholdRosin';
%  4 = 'intensityBinnedGradientThreshold'

params.ChannelIndex = [2,3];
params.GaussFilterSigma = 1;

% Just as a side note if you ever need the function name the process calls
% you can find it via
% funcHandle = MD.processes_{idxProc}.funName_
% where idxProc is the index of the process.
% and funcHandle is the function handle called by the process
% For instance the function called in this case is @thresholdMovie
% this function should give you documentation regarding the function
% parameters

%Resave the parameters
parseProcessParams(MD.processes_{end},params);

MD.processes_{1}.run(params); % run the process



MD.packages_{1}.setProcess(1,thresProc);

%% START Biosensors Package
%% These are the associated steps in the Biosensors (ie FRET processing GUI)
% These numbers are a bit hardwired which makes it a bit cumbersome
% Note if they are optional I do not run them here. But the you would
% run them in the same fashion.
% Note these are the same order as in the biosensor GUI via
% movieSelectorGUI:
% STEPS:
% 1) Dark Current Correction (Optional)
% 2) Shade Correction
% 3) Segmentation (performed above)
% 4) Background Mask
% 5) Mask Refinement (Optional)
% 6) Background Subtraction
% 7) Transformation (Optional)
% 8) Bleedthrough/crosstalk correction
% 9) Ratioing
% 10) Photobleach correction (Optional)
% 11) Ratio Output (I use my own functions for this step for more
% flexibility- just makes a series of .tifs and movies)

%% Shade correct : ie flat field correction - Step II in GUI

sCP = ShadeCorrectionProcess(MD);
MD.addProcess(sCP);
shadeParams = MD.processes_{2}.funParams_;
shadeParams.ChannelIndex = [2,3]; % only have for donor and acceptor channel

shadeDir = [imgFolder filesep 'ShadeCorrectCrop'];
shadeDirs{1} = [shadeDir filesep 'Donor'];
shadeDirs{2} = [shadeDir filesep 'FRET'];

% change params
shadeParams.ShadeImageDirectories = shadeDirs;
% set params
parseProcessParams(MD.processes_{2},shadeParams);
% run
MD.processes_{end}.run(shadeParams);

% associate to package
bioPack = BiosensorsPackage(MD);
MD.addPackage(bioPack);
MD.packages_{2}.setProcess(2,sCP); %  number of setProcess refers to step number as indicated above
MD.packages_{2}.setProcess(3,thresProc); % number of set process refers to step number as indicated above
% note a thresh process can got to multiple packages
%% Run Background Estimation
bmP = BackgroundMasksProcess(MD);
MD.addProcess(bmP);

% Get the threshold parameters so you can modify them
bkParams = MD.processes_{end}.funParams_;

bkParams.ChannelIndex = [2,3];
%bkParams.SegProcessIndex = [1];
% set params
parseProcessParams(MD.processes_{end},bkParams);
MD.processes_{end}.run(bkParams);


MD.packages_{2}.setProcess(4,bmP);
MD.save;
%% Run Background Subtraction
bsP = BackgroundSubtractionProcess(MD);
MD.addProcess(bsP);
bsParams= MD.processes_{end}.funParams_;
bsParams.ChannelIndex = [ 2, 3];
bsParams.MaskChannelIndex = [2,3];
parseProcessParams(MD.processes_{end},bsParams);
MD.processes_{end}.run(bsParams);

MD.packages_{2}.setProcess(6,bsP);

%% Run transformation if applicable - this is to align the channels
% load the transformation matrix

% trans = TransformationProcess(MD);
% MD.addProcess(trans);
% idxTrans =  find(cellfun(@(x) sum(strcmpi( x.name_,'transformation')), MD.processes_));
% transParams = MD.processes_{idxTrans}.funParams_;
% transParams.ChannelIndex = 2;
% transDir =  upDirectory(imgFolder,5);
% transParams.TransformFilePaths{1} = [transDir filesep 'transformation_matrix.mat'];
% transParams.TransformMasks = 0;
% parseProcessParams(MD.processes_{idxTrans},transParams);
% MD.processes_{idxTrans}.run(1);
%
%% Run Ratio
ratProcess = RatioProcess(MD);
MD.addProcess(ratProcess);
idxRat = find(cellfun(@(x) sum(strcmpi(x.name_,'Ratioing')),MD.processes_));
%
ratioParams=MD.processes_{idxRat}.funParams_;
ratioParams.SegProcessIndex = 1;
ratioParams.MaskChannelIndex = [2,2];
ratioParams.ChannelIndex = [3,2];
%ratioParams.SegProcessIndex = 1;

% note the option not to apply the masks I don't believe is available in
% the GUI : this is helpful if you want to apply masks to the ratio image
% using another function not in the  MD processes.
%ratioParams.ApplyMasks = 0; %
parseProcessParams(MD.processes_{idxRat},ratioParams);
MD.processes_{idxRat}.run(ratioParams);

MD.packages_{2}.setProcess(9,ratProcess);
% eventually want shade correction process...
% bioProcess = BiosensorPackage(

MD.save

% I use my own visualization to have more flexibility
% so I don't usually don't run the last step.

% you should be able to visualize everything you have run here via the
% command by loading the movieData file in movieSelectorGUI.
% The
%% Note you should now be able to load this into the movieSelectorGUI
% It will show you some exclamation point flags for steps 3 to end : as it
% is worried that
end

