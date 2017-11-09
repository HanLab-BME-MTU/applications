function script_FilamentRun_perWell(rowInd,rootDataFolder)
% run the filament analysis package for vim screen plates per well version

% take care of Matlab path
matlabDir = '/home2/lding/Ding/GitMatlab'; 
addpath(genpath(matlabDir),'-end'); 

% Set the root folder, as input
if(isempty(rootDataFolder))
 rootDataFolder ='/project/cellbiology/gdanuser/vimentin/tonyzhang/live cell nikon data/well plate screening/screen_20141204/raw data_analysis3';
end
% set the data folder according to the row index
% 1 for row A, 2 for B, and so on
rowName = char('A' + rowInd - 1);
rowFolderName = [ rowName,'_Row_Analysis'];

% load the MD file
MD = MovieData.load([rootDataFolder,filesep,rowFolderName,'movieData.mat']);

% run the filament analysis package
vimscreen_load_MD_run_filament_analysis_package(MD,[],'run_with_new_param',1,'save_old_data_tag','time'); 

% Mark the Completion of the network analysis

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp(['Finished filament segmentation for row ',rowName]);
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');

% run network analysis
load_MD_network_for_analysis(MD,[],30, 1,1, ones(16,1),1);

% Mark the Completion of the network analysis

disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
disp(['Finished network analysis for row ',rowName]);
disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
