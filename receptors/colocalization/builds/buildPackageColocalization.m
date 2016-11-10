% Simple script for building application packages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set output path for package build files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build on BioHPC - Linux
% out_dir = '/home2/s170480/matlab/builds/colocalize/'
% Build on Andrew's Windows 10 machine
out_dir = 'C:\Users\Andrew\GIT\matlab\builds\colocalization';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buildPackage('scriptGeneralColocalization.m',out_dir);
cd(out_dir); % check results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Test image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move test image to build package dir.
test_img = 'C:\Users\Andrew\Data\raw\Tony\Colocalization\test_0001.tif';
copyfile(test_img, out_dir); 
[a, b, c] = fileparts(test_img);
test_img = [b c];

% ---> Basically Zip up this directory now for the package.
zip('../colocalization.zip', '../colocalization/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the package build works!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;  % Clears out repo paths
addpath(genpath(out_dir)); % add the build package path
scriptGeneralColocalization(test_img);




