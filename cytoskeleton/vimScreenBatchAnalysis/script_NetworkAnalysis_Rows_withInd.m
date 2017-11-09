function script_NetworkAnalysis_Rows_withInd(rowInd,rootDataFolder)
% run the filament analysis package for vim screen plates

% take care of Matlab path
matlabDir = '/home2/zzhang/GitMatlabCode'; 
addpath(genpath(matlabDir),'-end'); 

% Set the root folder, as input
if(isempty(rootDataFolder))
 rootDataFolder ='/project/cellbiology/gdanuser/vimentin/tonyzhang/live cell nikon data/well plate screening/screen_20141204/raw data_analysis/';
end
% set the data folder according to the row index
% 1 for row A, 2 for B, and so on
rowName = char('A' + rowInd - 1);
rowFolderName = [ rowName,'_Row_Analysis'];

% load the MD file
MD = MovieData.load([rootDataFolder,filesep,rowFolderName,'movieData.mat']);

% run the filament analysis package
load_MD_network_for_analysis(MD,[],30, 1,1, [ones(6,1); zeros(10,1)],1)
