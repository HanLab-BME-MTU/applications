%%% Description
% This script shows an example of the use of the indexLatticeData function 
% that creates MovieList file that links movies/cells together. For each 
% movies that are described by a generic path, a  movieData structure 
% is created. They are then listed in aMovieList file.
%
% The example below collect all the files that correspond to the path:
%
% '/work/gdanuser/proudot/project/EB3-3D-track/data-analysis/four-phases/Prometaphase_cell*_deskew/deskewed-cropped/ch{ch}/*.tif'
%
% with * start meaning any possible string of character and {ch} meaning a channel number.
%
% For each movie, an analysis folder is created with a MovieData:  
%
% /work/gdanuser/proudot/project/EB3-3D-track/data-analysis/four-phases/Prometaphase_cell*_deskew/deskewed-cropped/analysis/
%
% Every movie is listed in a movieList file called 'prometaphaseCells.mat'
% Saved in 
% '/work/gdanuser/proudot/project/EB3-3D-track/data-analysis/four-phases/analysis'

%% %% USER INPUT
dataPath='/project/bioinformatics/Danuser_lab/shared/proudot/3d-vis/utrackPackage/scriptBased/UTrack-QD-v1/data/';
timeInterval=0.28;                % in secondes
axialPixelSize=400;               % in nm
lateralPixelSize=160;             % in nm    

movieList=indexLSFMData([dataPath '/Cell12Toy2/ch{ch}/*.tif'],dataPath,'movieListName','singleCellList.mat','axialPixelSize',axialPixelSize,'lateralPixelSize',lateralPixelSize,'timeInterval',timeInterval);

