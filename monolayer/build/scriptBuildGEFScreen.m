% Simple script for building application packages

% Build on BioHPC - Linux

if strcmp(computer('arch'),'win64')
	matlab_repo_root = 'C:\Users\Andrew\GIT\matlab\';
else
	matlab_repo_root = [getenv('HOME') filesep 'matlab'];
end

package_name = 'GEFScreen';
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
repo_dirs = fullfile(matlab_repo_root, {'extern'; 'applications/monolayer/'; 'common'});
for i = 1:length(repo_dirs)
    cur_dir = repo_dirs{i};
    disp(['Adding ' cur_dir]);
    addir(cur_dir);
end
ScriptIn(1) = {'testScript_assaf.m'};
ScriptIn(2) = {'export_fig.m'};
buildPackage(ScriptIn, out_dir);
cd(out_dir); % check results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Add license 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Be sure to modify the statement as needed.
% %{
% addCopyingStatement = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'addCopyingStatement']);
% CopyingStatement = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'CopyingStatement']);
% GPL_license = fullfile(matlab_repo_root,['documentation' filesep 'license' filesep 'GPL-License.txt']);
% copyfile(addCopyingStatement, out_dir);
% copyfile(CopyingStatement, out_dir);
% copyfile(GPL_license, out_dir);

% if strcmp(computer('arch'),'win64')
%     disp('Win64');
%     disp(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']);
%     winopen(out_dir);
%     uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']));
%     system('C:\Program Files\Git\git-bash.exe');
% else
% 	disp(['Run this code in the pop up: "bash addCopyingStatement ' package_name '  "' institution_name '"']);
% 	disp(['in this directory          : ' out_dir])
% 	uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' "' institution_name '"']));
% 	% system(['bash addCopyingStatement ' package_name ' ' institution_name]);
% end
% %}
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
	tmpdir = fullfile(tempname, package_name);
	out_dir = fullfile(tmpdir, 'analysis');
	mkdir(tmpdir);
	disp(['Created tmpdir ' tmpdir]);
	unzip(zip_file, tmpdir);
	addpath(genpath(tmpdir)); % add the build package path
	disp('Added tmpdir to path');
	mkdir(out_dir);
	disp(['Created output dir ' out_dir]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Gather Test image
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% move test image to test package dir.
	if strcmp(computer('arch'),'win64')
		test_img = 'C:\Users\Andrew\Data\raw\Assaf\Angeles_20150402_14hrs_5min_AA01_7.tif';
	else
		test_img = '/work/bioinformatics/s170480/Data/Assaf/Angeles_20150402_14hrs_5min_AA01_7.tif';
	end

	[a, b, c] = fileparts(test_img);
	test_img_name = fullfile(out_dir, [b c]);

	if ~strcmp(out_dir, a) 
		copyfile(test_img, out_dir); 	
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	cd(out_dir);
	testScript_assaf(test_img_name);
	uiwait(msgbox(['Please check the analysis output...:' out_dir]));
	restoredefaultpath;	
end

addpath(start_paths);
cd(start_dir);


