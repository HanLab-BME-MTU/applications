
%First step, create a MovieList with all the movie you want to group  
%Of course, it is expected by now that you already ran the windowing package

%Load the MovieList
ML = MovieList.load(pathForTheMovieListFile);

%Format the protrusion map
%You can:
%           remove outliers using the input outLevel
%           exclude windows by length using minLength
%           detrend using trendType
%           exclude pre-selected fucked windows using excludeWin
%           enforce units of nm/sec by setting scale to true

%Read the function's help
cellData = formatEdgeVelocity(ML,'minLength',minLen,'outLevel',outLevel,'scale',scale,'trendType',trendT);

%Again, read the edgeVelocityQuantification help
[cellData,dataSet] = edgeVelocityQuantification(ML);

%cellData is a structure array which contains all the information at the cell level and at the window level
%cellData(iCell)

%dataSet contains all the information for the whole data set

%After running these two functions, 2 folders are created.
%One of each cell, in the same folder where movieData is
%One for the movieList, in the same folder where the movieList is