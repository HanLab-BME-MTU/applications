% MATLAB-based testing/performance suite for BioSensors Package
% Andrew R. Jamieson 2017
% Test BioSensors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Preconditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['PWD:', pwd]);
% Dump path for debugging
s_path = strsplit(path,':');
s_match = (cellfun(@(x) regexp(x,'toolbox'), s_path, 'UniformOutput', false))';
matlab_paths = s_path(cellfun(@isempty, s_match))';
disp('    [MATLAB] current top-level paths....');
% disp(matlab_paths);

disp(['Java heap max: ' num2str(java.lang.Runtime.getRuntime.maxMemory/1e9) 'GB'])
disp('Starting biosensor script');

%----Initialization of temp dir
package_name = 'BioSensor';
t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
mkdir(tmpdir);

cd(tmpdir);

% Download test data for BioSensor (Use BioHPC internal cloud account danuserweb)
url = 'https://lamella.biohpc.swmed.edu/index.php/s/HOq0LUsukadEx1l/download';
zipPath = fullfile(tmpdir, 'Neurite.zip');
urlwrite(url, zipPath);
unzip(zipPath, tmpdir);

saveFolder = [tmpdir filesep 'Neurite' filesep 'Analysis'];
imgFolder = [tmpdir filesep 'Neurite' filesep 'Channels'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channels = {'C1_mCherry','C2_Donor','C3_FRET'};
for iCh = 1:3
    % Create a channels object
    imgFolderC = [imgFolder filesep channels{iCh}];
    channel(iCh) = Channel(imgFolderC);
end

MD = MovieData(channel, saveFolder);
MD.setPath(saveFolder);
MD.setFilename('movieDataBiosensors.mat');
MD.sanityCheck;
MD.save;
MD.reset();


% Generate Package
Package_ = BiosensorsPackage(MD);
% Set-up analysis infrastructure via command line interace
MD.addPackage(Package_);

% Specify processes to include in test.
stepNames = Package_.getProcessClassNames;
disp('===============================');
disp('Available Package Process Steps');
disp('===============================');
disp(stepNames');
steps2Test = [3, 4, 5, 6, 9];
assert(length(Package_.processes_) >= length(steps2Test));
assert(length(Package_.processes_) >= max(steps2Test));
disp('Selected Package Process Steps');

for i=steps2Test
	disp(['Step ' num2str(i) ': ' stepNames{i}]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START Biosensors Package
% These are the associated steps in the Biosensors (ie FRET processing GUI)
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 2: Shade Correction
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 2: Shade Correction');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 2;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;

funParams.ChannelIndex = [2,3]; % only have for donor and acceptor channel
shadeDir = [imgFolder filesep 'ShadeCorrectCrop'];
shadeDirs{1} = [shadeDir filesep 'Donor'];
shadeDirs{2} = [shadeDir filesep 'FRET'];

funParams.ShadeImageDirectories = shadeDirs;
MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();


%% Step 3: Threshold & Segmentation (otsu)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 3: Threshold & Segmentation');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 3;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;

disp('Available methods');
disp({MD.processes_{2}.getMethods.name}')

% Set the Thresholding Parameters
funParams.MethodIndx = 1;
% Here are the methods for the thresholdProcess
% MD.processes_{1}.getMethods
%  1 = 'thresholdFluorescenceImage'; (ie minMax)
%  2 = 'thresholdOtsu';
%  3 = 'thresholdRosin';
%  4 = 'intensityBinnedGradientThreshold'

funParams.ChannelIndex = [2,3];
funParams.GaussFilterSigma = 1;

%Resave the parameters
MD.getPackage(1).getProcess(step_).setPara(funParams);

%% Step 3a: MinMax Threshold
disp('% ------- Step 3a: MinMax --------------');
step_ = 3;
MD.getPackage(1).processes_{step_}.run();

%% Step 3b: Rosin Threshold
disp('% ------- Step 3b: Rosin --------------');
step_ = 3;
funParams = MD.getPackage(1).processes_{step_}.funParams_;
funParams.MethodIndx = 3;
MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();

%% Step 3c: Gradient Threshold
disp('% ------- Step 3c: BinnedGradient --------------');
step_ = 3;
funParams = MD.getPackage(1).processes_{step_}.funParams_;
funParams.MethodIndx = 4;
MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();

%% Step 3d: Otsu Threshold
disp('% ------- Step 3d: Otsu --------------');
step_ = 3;
funParams = MD.getPackage(1).processes_{step_}.funParams_;
funParams.MethodIndx = 2;
funParams.GaussFilterSigma = 2;
MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();


%% Step 4: Background mask
% note a thresh process can got to multiple packages
% Run Background Estimation
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 4: Background mask');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 4;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;
funParams.ChannelIndex = [2, 3];

MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();


%% Step 5: Mask Refinement
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 5: % Mask Refinement');% Should select automaticall last one
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 5;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;
funParams.ChannelIndex = [2, 3];

MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();


%% Step 6: Background Subtraction
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 6: Background Subtraction');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 6;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;

funParams.ChannelIndex = [2, 3];
funParams.MaskChannelIndex = [2, 3];

MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();
MD.save

%% Step 9: Ratioing
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Step 9: Ratioing');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

step_ = 9;
if isempty(MD.getPackage(1).processes_{step_})
    MD.getPackage(1).createDefaultProcess(step_);
end
funParams = MD.getPackage(1).processes_{step_}.funParams_;

funParams.ChannelIndex = [3, 2];
funParams.MaskChannelIndex = [3, 2];

MD.getPackage(1).getProcess(step_).setPara(funParams);
MD.getPackage(1).processes_{step_}.run();
MD.save