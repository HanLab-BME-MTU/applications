%%% Description
% This script find the tiff stack under the folder <dataPath> and write the
% associated movie management file MovieData used for pixel reading, unit
% and output path management
%
% The file tree under <dataPath> must adopt the following structure:
%   <dataPath> / conditionName / CellNumber / chNumber / *.tif
%
% See data in <../data> for an example
%
% To load the movies use <loadMoviesManagementFile.m>.

%% %% USER INPUT
dataPath='/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1';
timeInterval=1;         % in secondes
lateralPixelSize=100;   % in nm
axialPixelSize=235;     % in nm

%% %% MOVIE CREATION PROCESS
movieListArray=createMovieList(dataPath,'timeInterval',timeInterval,'lateralPixelSize',lateralPixelSize,'axialPixelSize',axialPixelSize);
[~,xpName]=fileparts(dataPath);
save([dataPath filesep xpName '.mat'],'movieListArray');
