packageUnzipLocation = "<location where you unzip package file>"

addpath(genpath(packageUnzipLocation)); % add the build package path
disp('Added packageUnzipLocation to path');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For GUI based workflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


movieSelectorGUI;
disp('load the movies input movieSelector');
disp('then select the segmentation package');
disp('and follow the steps');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For script based workflow see below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define output directories
out_dir = 'where you want your output analysis to go'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define input image locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create MovieData
saveFolder = out_dir;
ch1 = Channel('mDiaActin_mini/images_mDia1');  %% <<< define input image dir locations
ch2 = Channel('mDiaActin_mini/images_actin'); %% <<< define input image dir locations
% Constructor needs an array of channels and an output directory (for analysis)
MD = MovieData([ch1 ch2], saveFolder);
MD.setPath(saveFolder);
MD.setFilename('movieData.mat');
% Check image size and number of frames are consistent.
% Save the movie if successfull
MD.sanityCheck; % 
% Set some additional movie properties
MD.numAperture_=1.4;
MD.pixelSize_=215; % in nm after binning
MD.timeInterval_=5;% in sec
MD.camBitdepth_=16;
MD.notes_='Created for test purposes';
% Save the movieData
MD.save;


% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note you can use the GUI after creating the MovieData too

% movieSelectorGUI('MD', MD);
% % or 
% SegmentationPackage(MD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Run Threshold
thresProc = ThresholdProcess(MD);
MD.addProcess(thresProc);
% Create a segmentation package
segPackage = SegmentationPackage(MD);
MD.addPackage(segPackage);
% % Associate the threshold process to the package
MD.packages_{1}.setProcess(1,thresProc);
% Get the threshold parameters so you can modify them
params = MD.processes_{1}.funParams_;
% Set the Thresholding Parameters
params.MethodIndx = 1;
params.GaussFilterSigma = 1;
% Resave the parameters
parseProcessParams(MD.processes_{end}, params);
% Run the process
MD.processes_{1}.run(); %
MD.save

%% Refine the segmentation mask
MD.addProcess(MaskRefinementProcess(MD))
MD.packages_{1}.setProcess(2,MaskRefinementProcess(MD));   
MD.processes_{2}.run()

disp('Finish Windowing test run script successfully');