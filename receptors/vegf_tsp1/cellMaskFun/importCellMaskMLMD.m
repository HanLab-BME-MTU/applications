function importCellMaskMLMD(MLMD)
%importCellMaskMLMD imports cell masks for a list of movies
%
%SYNOPSIS importCellMaskMLMD(MLMD)
%
%INPUT MLMD: MovieList or MovieData object for movie(s) to be analyzed
%
%OUPUT Output is saved in directory ImportedCellMask belonging to each
%      analyzed movie.
%
%Khuloud Jaqaman, March 2015

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--importCellMaskMLMD: Function needs at least 1 input argument!');
    return
end

%% Analysis

%determine if input is a MovieList or MovieData object
if isa(MLMD,'MovieList') %if it is a movieDist
    
    listFlag = 1;
    
    %rename to ML
    ML = MLMD;
    clear MLMD
    
    %get number of movies
    numMovies = length(ML.movieDataFile_);
    
else %if it is a movieData
    
    listFlag = 0;
    
    %rename to MD
    MD = MLMD;
    clear MLMD
    
    %assign number of movies to 1
    numMovies = 1;
    
end

%go over all movies and import their cell mask
for iM = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{iM});
    end
    
    %add cell mask import process if never run before
    iProc = MD.getProcessIndex('ImportCellMaskProcess',1,0);
    if isempty(iProc)
        iProc=numel(MD.processes_)+1;
        MD.addProcess(ImportCellMaskProcess(MD));
    end
    
    %define function parameters
    
    %general
    funParams = MD.getProcess(iProc).funParams_;
    funParams.ChannelIndex = 1;
    
    %function-specific
    %none for now
    
    %general again
    parseProcessParams(MD.getProcess(iProc),funParams);
    
    %% Run process
    cellfun(@run,MD.processes_(iProc));
    
end

