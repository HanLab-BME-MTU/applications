% Simple script for testing application packages

if strcmp(computer('arch'),'win64')
	% matlab_repo_root = '/home2/s170480/matlab/'
	matlab_repo_root = 'C:\Users\Andrew\GIT\matlab\';
	% Build on Andrew's Windows 10 machine
else
	matlab_repo_root = [getenv('HOME') filesep 'matlab'];
end

package_name = 'GEFScreen';
institution_name = 'Danuser Lab - UTSouthwestern';
start_paths = path;
start_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set output path for package build files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice = questdlg(['Run test in TMP dir?'],'Question..','Yes','No','Yes');

if strcmp(choice, 'Yes')
	tmpdir = fullfile(tempname, package_name);
	mkdir(tmpdir);
	out_dir = tmpdir;
	disp(['Created tmpdir ' out_dir]);
elseif strcmp(choice, 'No')
	out_dir = uigetdir(matlab_repo_root);
else
	error('No Choice made!');		
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Test image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move test image to test package dir.
if strcmp(computer('arch'),'win64')
	test_img = 'C:\Users\Andrew\Data\raw\Assaf\Angeles_20150402_14hrs_5min_AA01_7.tif';
else
	test_img = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/testSW/Angeles_20150402_14hrs_5min_AA01_7.tif';
end

[a, b, c] = fileparts(test_img);
test_img_name = fullfile(out_dir, [b c]);

if ~strcmp(out_dir, a) 
	copyfile(test_img, out_dir); 	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add only appropriate paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;  % Clears out repo paths
% Add necesary repos for testing
repo_dirs = fullfile(matlab_repo_root, {'extern'; ['applications' filesep 'monolayer']; 'common'});
for i = 1:length(repo_dirs)
    cur_dir = repo_dirs{i};
    disp(['Adding ' cur_dir]);
    addir(cur_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run test script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testScript_assaf(test_img_name);

restoredefaultpath;
addpath(path);



