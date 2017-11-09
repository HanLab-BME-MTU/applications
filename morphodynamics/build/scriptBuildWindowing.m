% Simple script for building application packages

% Build on BioHPC - Linux

if strcmp(computer('arch'),'win64')
	matlab_repo_root = 'C:\Users\Andrew\GIT\matlab\';
else
	matlab_repo_root = [getenv('HOME') filesep 'matlab'];
end

package_name = 'WindowingPackage';
institution_name = 'Danuser Lab - UTSouthwestern';
t_stamp = datestr(now,'ddmmmyyyyHHMM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set output path for package build files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build_dir = fullfile(matlab_repo_root, 'builds');
out_dir = fullfile(matlab_repo_root, ['builds' filesep package_name filesep t_stamp]);
zip_file = fullfile(matlab_repo_root, ['builds' filesep package_name t_stamp '.zip']);
start_paths = path;
start_dir = pwd;
mkdir(out_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;  % Clears out repo paths
% Add necesary repos for building 
repo_dirs = fullfile(matlab_repo_root, {'extern'; 'applications/morphodynamics/'; 'common'});
for i = 1:length(repo_dirs)
    cur_dir = repo_dirs{i};
    disp(['Adding ' cur_dir]);
    addir(cur_dir);
end
ScriptIn(1) = {'WindowingPackage.m'};
ScriptIn(2) = {'SegmentationPackage.m'};
buildPackage(ScriptIn, out_dir);
cd(out_dir); % check results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add license 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Be sure to modify the statement as needed.
addCopyingStatement = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'addCopyingStatement']);
CopyingStatement = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'CopyingStatement']);
GPL_license = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'GPL-License.txt']);
copyfile(addCopyingStatement, out_dir);
copyfile(CopyingStatement, out_dir);
copyfile(GPL_license, out_dir);

if strcmp(computer('arch'),'win64')
    disp('Win64');
    disp(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']);
    winopen(out_dir);
    uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']));
    system('C:\Program Files\Git\git-bash.exe');
else
	disp(['Run this code in the pop up: "bash addCopyingStatement ' package_name '  "' institution_name '"']);
	disp(['in this directory          : ' out_dir])
	uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' "' institution_name '"']));
	% system(['bash addCopyingStatement ' package_name ' ' institution_name]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zip up package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---> Basically Zip up this directory now for the package.
cd(build_dir);
if exist(zip_file, 'file')
	zip_file = [zip_file(1:end-4) '1_.zip'];
end
zip(zip_file, out_dir);
msgbox(['Package zipped here ' zip_file])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the package build works in tmp dir!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = questdlg(['Run test on zip package?'],'Question..','Yes','No','Yes');

if strcmp(choice, 'Yes')
	restoredefaultpath;  % Clears out repo paths
    t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
    tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
	mkdir(tmpdir);
	disp(['Created tmpdir ' tmpdir]);
	unzip(zip_file, tmpdir);
	addpath(genpath(tmpdir)); % add the build package path
	disp('Added tmpdir to path');
	out_dir = fullfile(tmpdir, 'analysis');
	mkdir(out_dir);
	disp(['Created output dir ' out_dir]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Gather Test image
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% move test image to test package dir.
	url = 'https://lamella.biohpc.swmed.edu/index.php/s/AfiyWzGV4bUYknV/download'
	zipPath = fullfile(tmpdir, 'mDiaActin_mini.zip');
	urlwrite(url, zipPath);
	unzip(zipPath, tmpdir);	
    delete(zipPath);
	cd(tmpdir);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create MovieData
	saveFolder = 'analysis';
	ch1 = Channel('mDiaActin_mini/images_mDia1');
	ch2 = Channel('mDiaActin_mini/images_actin');
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
    
    %% Run the protrusion vectors
	MD.addProcess(ProtrusionProcess(MD));
	% Get the Protrusion Parameters so you can modify them
	protParams = MD.processes_{end}.funParams_;
	% Set the Protrusion Parameters
	protParams.ChannelIndex = 1; % 
	protParams.SegProcessIndex = 2; % (refined mask)
	% Resave the parameters (if reconfigured)
	parseProcessParams(MD.processes_{end},protParams);
	% Run the process
	MD.processes_{end}.run();
	MD.save
	
    %% Run the Windowing
	MD.addProcess(WindowingProcess(MD));
	% Get the Windowing Parameters so you can modify them
	windParams = MD.processes_{end}.funParams_;
	parseProcessParams(MD.processes_{end},windParams);
	% Run
	MD.processes_{end}.run();
	MD.save;

	%% Run protrusion sampling process
	MD.addProcess(ProtrusionSamplingProcess(MD));
	% Get the Windowing Parameters so you can modify them
    protParams = MD.processes_{end}.funParams_;
% 	protParams.OutputDirectory = [MD.outputDirectory_ filesep 'protrusion_samples_' outName];
	parseProcessParams(MD.processes_{end},protParams);
	MD.processes_{end}.run();

	MD.save;

	disp('Finish Windowing test run script successfully');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	uiwait(msgbox(['Please check the analysis output...:' out_dir]));
	restoredefaultpath;	
end

addpath(start_paths);
cd(start_dir);


