function makeBiosensorsPackage(outDir,bioDir,comDir)
%MAKEBIOSENSORSPACKAGE builds the biosensor package
% 
% makeBiosensorsPackage(outDir,bioDir,comDir)
% 
% THis function copies all the files needed to run the biosensosrs data
% processing package into a single folder for upload to the website. It
% assumes you have checked out all the files from the SVN repository, and
% that they are all in your matlab path.
% 
% Input:
% 
%   outDir - The directory to copy all the package files to.
%
%   bioDir - the directory you have checked out the biosensors project to.
% 
%   comDir - The directory you have checked out the common project to.
% 
% Hunter Elliott
% 12/2010
%

if nargin < 1 || isempty(outDir)    
    outDir = uigetdir(pwd,'Select output dir:');
end
if nargin < 2 || isempty(bioDir)    
    bioDir = uigetdir(pwd,'Select biosensors dir:');
end
if nargin < 3 || isempty(comDir)    
    comDir = uigetdir(pwd,'Select common dir:');
end

disp('Getting file list...')

%Get m files from biosensors directory
bioFuns = dir([bioDir filesep '*.m']);
bioFuns = {bioFuns(:).name}';
%Remove this function from the list. Fucking recursion.
bioFuns(strcmp(bioFuns,mfilename)) = [];


%Get everything these depend on also
bioFuns = depfun_notoolbox(bioFuns);
%Get extra non-function files (figures for GUIs)
bioExtras = dir([bioDir filesep '*.fig']);
bioExtras = {bioExtras(:).name}';
bioExtras = cellfun(@(x)([bioDir filesep x]),bioExtras,'UniformOutput',false);

%Get the movie management functions
comDir = [comDir filesep 'toolfun' filesep 'movieManagement'];
comFuns = dir([comDir filesep '*.m']);
comFuns = {comFuns(:).name}';
%Add the few odd files that are in other areas of common
comFuns = vertcat(comFuns,{'bleedthroughCorrectMovie.m','backgroundSubtractMovie.m','refineMovieMasks.m','separateNumberedFiles.m'}');

%..and everything these functions depend on
comFuns = depfun_notoolbox(comFuns);

%Get fig files from common also
comExtras = dir([comDir filesep '*.fig']);
comExtras = {comExtras(:).name}';
comExtras{end+1} = 'lccbGuiIcons.mat';
comExtras = cellfun(@(x)([comDir filesep x]),comExtras,'UniformOutput',false);

%Check and display toolbox dependency
disp('Checking toolbox dependency...')
tbs = toolboxesUsed(vertcat(bioFuns,comFuns));

disp('The package uses the following toolboxes:')
disp(tbs)

if ~exist(outDir,'dir')
    mkdir(outDir)
end

allFiles = vertcat(bioFuns,bioExtras,comFuns,comExtras);
nFiles = numel(allFiles);

disp(['Copying all '  num2str(nFiles) ' files ...'])

for j = 1:nFiles

    iLFS = max(regexp(allFiles{j},filesep));
    copyfile(allFiles{j},[outDir filesep allFiles{j}(iLFS+1:end)]);

end

disp(['Finished. Wrote package to ' outDir])
