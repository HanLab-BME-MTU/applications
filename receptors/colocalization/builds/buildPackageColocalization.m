% Simple script for building application packages

% Build on BioHPC - Linux
% matlab_repo_root = '/home2/s170480/matlab/'
% Build on Andrew's Windows 10 machine

matlab_repo_root = [getenv('HOME') filesep 'matlab'];
package_name = 'ColocP2C';
institution_name = 'Jaqaman & Danuser Labs - UTSouthwestern';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set output path for package build files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build_dir = fullfile(matlab_repo_root, 'builds');
out_dir = fullfile(matlab_repo_root, ['builds' filesep package_name]);
zip_file = fullfile(matlab_repo_root, ['builds' filesep package_name '.zip']);
start_paths = path;
start_dir = pwd;

if exist(out_dir, 'dir') == 7
    choice = questdlg(['Remove old build dir?'],'Question..','Yes','No','Yes');
    if strcmp(choice, 'Yes')
    	rmdir(out_dir, 's');
    	pause(.25);
    	rehash(); % weird issue with matlab rmdir?
    end
end
rehash();
mkdir(out_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;  % Clears out repo paths
% Add necesary repos for building 
repo_dirs = fullfile(matlab_repo_root, {'extern'; 'applications'; 'common'});
for i = 1:length(repo_dirs)
    cur_dir = repo_dirs{i};
    disp(['Adding ' cur_dir]);
    addir(cur_dir);
end
cellScript{1} = 'scriptGeneralColocalization.m';
cellScript{2} = 'scriptMultiChannelColocalization.m';
buildPackage(cellScript, out_dir);
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
    uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']));
    winopen(out_dir);
    system('C:\Program Files\Git\git-bash.exe');
else
	disp(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']);
	disp(['in this directory          : ' out_dir])
	uiwait(msgbox(['Run this code in the pop up: "bash addCopyingStatement ' package_name ' ' institution_name '"']));
	% system(['bash addCopyingStatement ' package_name ' ' institution_name]);
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Test image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move test image to build package dir.
test_img = '/home2/avega/Documents/test_0001.tif';
copyfile(test_img, out_dir); 
[a, b, c] = fileparts(test_img);
test_img_name = [b c];
%}
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
	java_tmpdir = char(java.lang.System.getProperty('java.io.tmpdir'));
	uuid = java.util.UUID.randomUUID();
	uuid = char(uuid.toString());
	tmpdir = fullfile(java_tmpdir, uuid);
	mkdir(tmpdir);
	% Unzip test imagesa
	unzip(zip_file, tmpdir);
	% Run
	restoredefaultpath;  % Clears out repo paths
	addpath(genpath(tmpdir)); % add the build package path
	cd([tmpdir filesep package_name]);
	scriptGeneralColocalization
    scriptMultiChannelColocalization
    cd(start_dir);
	restoredefaultpath;
	try
        rmdir(tmpdir, 's');
    	pause(.25);
    	rehash(); % weird issue with matlab rmdir?
    catch ee
         disp(['unable to delete tmp dir' tmpdir]);
    end
end

addpath(start_paths);
cd(start_dir);


