
addLicense = false;

if strcmp(computer('arch'),'win64')
    matlab_repo_root = 'C:\Users\Andrew\GIT\matlab\';
    data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
else
    matlab_repo_root = [getenv('HOME') filesep 'matlab'];
    data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
end

package_name = 'CMEAnalysisPackage_alpha';
institution_name = 'Danuser Lab - UTSouthwestern';
t_stamp = datestr(now,'ddmmmyyyyHHMM');


cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);

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
repo_dirs = fullfile(matlab_repo_root, {'extern'; 'applications/endocytosis/'; 'common'});
for i = 1:length(repo_dirs)
    cur_dir = repo_dirs{i};
    disp(['Adding ' cur_dir]);
    addir(cur_dir);
end
cellScript{1} = 'testCMEanalysis.m';
cellScript{2} = 'sortTiffStacks.m'; % helper function
buildPackage(cellScript, out_dir);
cd(out_dir); % check results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add license 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if addLicense
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
    testScript_path = evalc('which testCMEanalysis')
    restoredefaultpath;  % Clears out repo paths
    t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
    tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
    mkdir(tmpdir);
    copyfile(testScript_path, tmpdir);     
    disp(['Created tmpdir ' tmpdir]);
    unzip(zip_file, tmpdir);
    addpath(genpath(tmpdir)); % add the build package path
    disp('Added tmpdir to path');
    out_dir = fullfile(tmpdir, 'analysis');
    mkdir(out_dir);
    disp(['Created output dir ' out_dir]);
%   cd(out_dir);
    results = runtests('testCMEanalysis.m');
    disp(results.table);
    uiwait(msgbox(['Please check the analysis output...:' out_dir]));
    restoredefaultpath; 
end

addpath(start_paths);
cd(start_dir);


