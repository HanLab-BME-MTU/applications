function [slist, dataProperties, optData] = detectSpots(rawMovieName, filteredMovieName, dataProperties, verbose, options)
%DETECTSPOTS is a wrapper for the static segmentation in the chromdyn package
%
% SYNOPSIS: [slist, optData] = detectSpots(rawMovieName, filteredMovieName, dataProperties, options)
%
% INPUT rawMovieName:       filename of raw movie or handle to an imaris
%                           session where the raw movie has been loaded
%                           already. Please load the movie into imaris
%                           outside of DETECTSPOTS if it can't be read by
%                           cdLoadMovie.
%		filteredMovieName:  (opt) filename of filtered movie. If
%                           filteredMovie ends in ".new", movie will be
%                           filtered and saved under the given name minus
%                           the ".new". (filtered movies should have the
%                           extension ".fim"!)
%		dataProperties:     dataProperties-structure describing the movie
%       verbose:            (opt) 0: none / {1} just waitbar / 2 more output
%		options:            optional input for debugging or extensions
%
% OUTPUT slist:             1-by-nTimepoints structure with coordinates,
%                           amplitudes, uncertainties of detected spots
%        dataProperties     dataProperties structure with added field
%                           "amplitudeCutoff"
%        optData:           additional data for debugging
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: Jonas Dorn
% DATE: 07-Feb-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================
% TEST INPUT
%=================


% check number of input arguments
if nargin < 3
    error('not enough input arguments!')
end
% set verbosity
if nargin < 4 || isempty(verbose)
    verbose = 1;
end
% make debug optional
if nargin < 5 || isempty(options)
    options = [];
end

% check for the existence of raw and filtered movies. Also check whether
% the rawMovieName is not actually an imaris handle
if strcmp(class(rawMovieName),'COM.Imaris_Application')
    % load movie via imaris
    movieLoader = 'imaris';
    imarisHandle = rawMovieName;
elseif ~exist(rawMovieName,'file')
    error('rawMovie (''%s'') not found',rawMovieName)
else
    movieLoader = 'cdLoadMovie';
end

% check for filteredMovie. If '.new', we want to filter
if isempty(filteredMovieName)
    error(['if there is no filtered movie yet, specify a moviename'...
        'with full extension ending in ''.new''!'])
end
if strcmp(filteredMovieName(end-3:end),'.new')
    % make a new filteredMovie
    doFilter = 1;
    filteredMovieName = filteredMovieName(1:end-4);
    % make sure we get the right extension
    if ~strcmp(filteredMovieName(end-3:end),'.fim')
        filteredMovieName = [filteredMovieName,'.fim'];
    end
elseif ~exist(filteredMovieName,'file')
    error('filteredMovie (''%s'') not found',filteredMovieName)
else
    doFilter = 0;
end


% check that dataProperties is a structure
if ~isstruct(dataProperties)
    error('dataProperties must be a structuer')
end

% assign optData
optData = [];

%====================

%=================================
% PREPARE MOVIES
%=================================

% find maximum allowed data size
if isfield(dataProperties,'maxSize')
    loadOptions.maxSize = dataProperties.maxSize;
else
    % set maximum size to 100 Mb
    loadOptions.maxSize = 100000000;
end

% check for path in rawMovieName
if strcmp(movieLoader,'cdLoadMovie')
[moviePath,movieName,extension] = fileparts(rawMovieName);
if isempty(moviePath)
    moviePath = pwd;
    rawMovieName = fullfile(moviePath, [movieName, extension]);
end
end
% check for pathName
[fmoviePath,fmovieName,fextension] = fileparts(filteredMovieName);
if isempty(fmoviePath)
    fmoviePath = pwd;
    filteredMovieName = fullfile(fmoviePath, [fmovieName, fextension]);
end


%=================================


%=================================
% FIND COORDINATES
%=================================

% preassign coordinates
coordinates(1:dataProperties.movieSize(4),1) = ...
    struct('sp',[],'COM',[]);

% load movie. If there isn't a filtered movie already, load the raw movie
if doFilter
    % load raw movie
    switch movieLoader
        case 'cdLoadMovie'
            [rawMovie, movieHeader, loadStruct] = ...
                cdLoadMovie({rawMovieName,'corr/raw'}, [], loadOptions);
            % check if there are any leading darkframes we need to subtract
        deltaFrames = loadStruct.loadedFrames(1) - 1;
        case 'imaris'
            [rawMovie,movieSize,movieName,...
                moviePath,movieHeader,imarisHandle,loadStruct] = ...
                imarisImread(imarisHandle,[],[],loadOptions.maxSize);
            deltaFrames = 0;
    end

else
    % load filtered movie
    [filteredMovie, movieHeader, loadStruct] = ...
        cdLoadMovie({filteredMovieName,'filtered'}, [], loadOptions);
    deltaFrames = 0;

end

% loop through the movie to collect local maxima
done = 0;
while ~done

    % filter first if necessary
    if doFilter
        % filter movie
        filteredMovie = filtermovie(rawMovie,dataProperties.FILTERPRM);
        clear rawMovie

        % save filtered movie to file
        writemat(filteredMovieName,...
            filteredMovie,1,5);
    end

    % find local maxima
    coordinates(loadStruct.loadedFrames - deltaFrames) = ...
        spotfind(filteredMovie,dataProperties);
    clear filteredMovie

    % load more
    if isempty(loadStruct.frames2load)
        done = 1;
    else
        if doFilter
            switch movieLoader
                case 'cdLoadMovie'
                    [rawMovie, movieHeader, loadStruct] = ...
                        cdLoadMovie(loadStruct.movieType, [], loadStruct);
                case 'imaris'
                    [rawMovie,dummy,dummy,...
                        dummy,dummy,dummy,loadStruct] = ...
                        imarisImread(loadStruct);
            end
        else
            [filteredMovie, movieHeader, loadStruct] = ...
                cdLoadMovie(loadStruct.movieType, [], loadStruct);
        end
    end % load more

end % while loop find local max

%===========================================



%===========================================
% MIXTURE MODEL FITTING
%===========================================

% first, we need to find the amplitude cutoff. Then we can go on to do the
% 'real' mixture model fitting. Of course, if there is already an
% amplitudeCutoff, we don't need to go over it again

if ~isfield(dataProperties,'amplitudeCutoff') || ...
        dataProperties.amplitudeCutoff == 0

    % if we load from imaris, rawMovieName is the imaris handle
    [dataProperties, testRatios] = ...
        detectSpots_MMF_findAmplitudeCutoff(...
        rawMovieName, coordinates, dataProperties, movieLoader, verbose);
else
    % set testRatios to []
    testRatios = [];
end

% now, we can get the slists. Loop through rawMovie
% load raw movie
switch movieLoader
    case 'cdLoadMovie'
        [rawMovie, movieHeader, loadStruct] = ...
            cdLoadMovie({rawMovieName,'corr/raw'}, [], loadOptions);
        % check if there are any leading darkframes we need to subtract
        deltaFrames = loadStruct.loadedFrames(1) - 1;
    case 'imaris'
        [rawMovie,movieSize,movieName,...
            moviePath,movieHeader,imarisHandle,loadStruct] = ...
            imarisImread(imarisHandle,[],[],loadOptions.maxSize);
        deltaFrames = 0;
end

% preassign slist
slist(1:movieHeader.numTimepoints) = ...
    struct('sp',[],...
    'statistics',[],...
    'parms',[],...
    'COM',[]);


done = 0;
while ~done

    % prepare data for mixture model fitting
    loadedFrames = loadStruct.loadedFrames-deltaFrames;

    % do fitting
    slist(loadedFrames) = ...
        detectSpots_MMF_main(rawMovie,coordinates(loadedFrames),...
        dataProperties,testRatios(loadedFrames),verbose);

    % load more
    if isempty(loadStruct.frames2load)
        done = 1;
    else
        
        switch movieLoader
            case 'cdLoadMovie'
                [rawMovie, movieHeader, loadStruct] = ...
                    cdLoadMovie(loadStruct.movieType, [], loadStruct);
            case 'imaris'
                [rawMovie,dummy,dummy,...
                    dummy,dummy,dummy,loadStruct] = ...
                    imarisImread(loadStruct);
        end
    
    end % load more

end % while loop find spots
