% Simple script for building application packages

% Build on BioHPC - Linux
% matlab_repo_root = '/home2/s170480/matlab/'
% Build on Andrew's Windows 10 machine
matlab_repo_root = 'C:\Users\Andrew\GIT\matlab\';
package_name = 'colocalization';
institution_name = 'UTSouthwestern';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set output path for package build files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = fullfile(matlab_repo_root, ['builds' filesep 'colocalization'])
if isdir(out_dir)
    rmdir(out_dir, 's');
end
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

buildPackage('scriptGeneralColocalization.m', out_dir);
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
system(['bash addCopyingStatement ' package_name ' ' institution_name]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Test image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move test image to build package dir.
test_img = 'C:\Users\Andrew\Data\raw\Tony\Colocalization\test_0001.tif';
copyfile(test_img, out_dir); 
[a, b, c] = fileparts(test_img);
test_img = [b c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zip up package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---> Basically Zip up this directory now for the package.
cd('../');
zip('colocalization.zip', 'colocalization');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the package build works in tmp dir!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
java_tmpdir = char(java.lang.System.getProperty('java.io.tmpdir'));
uuid = java.util.UUID.randomUUID();
uuid = char(uuid.toString());
tmpdir = fullfile(java_tmpdir, uuid);
mkdir(tmpdir);
% Unzip test imagesa
unzip('colocalization.zip', tmpdir);
% Run
restoredefaultpath;  % Clears out repo paths
cd([tmpdir filesep 'colocalization']);
addpath(genpath(tmpdir)); % add the build package path
scriptGeneralColocalization(test_img);
cd(out_dir);
rmdir(tmpdir, 's');



