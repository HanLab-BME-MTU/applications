function [slist, dataProperties, testRatios, debugData] = detectSpots(rawMovieName, filteredMovieName, dataProperties, verbose, options)
%DETECTSPOTS is a wrapper for the static segmentation in the chromdyn package
%
% SYNOPSIS: [slist, dataProperties, testRatios, debugData] = detectSpots(rawMovieName, filteredMovieName, dataProperties, options)
%
% INPUT rawMovieName:       filename of raw movie or handle to an imaris
%                           session where the raw movie has been loaded
%                           already. Please load the movie into imaris
%                           outside of DETECTSPOTS if it can't be read by
%                           cdLoadMovie.
%		filteredMovieName:  filename of filtered movie.
%		dataProperties:     dataProperties-structure describing the movie
%       verbose:            (opt) 0: none / {1} just waitbar / 2 more output
%		options:            optional input for debugging or extensions
%
% OUTPUT slist:             1-by-nTimepoints structure with coordinates,
%                           amplitudes, uncertainties of detected spots
%        dataProperties     dataProperties structure with added field
%                           "amplitudeCutoff"
%        testRatios         nTimepoints-by-1 cell with the testRatios for
%                           the amplitudes of the fitted spots
%        debugData:         additional data for debugging
%
% REMARKS  It is also possible to pass the raw (and maybe filtered) movie
%          instead of the names. This is not recommended, however.
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

% debug vector with:
% 1: return F-test values for old and new degrees of freedom
% 2: return residual images
debug = 0;


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
else
    if isfield(options,'debug')
        debug = options.debug;
    end
end

% update dataProperties
dataProperties = defaultDataProperties(dataProperties);

% check for the existence of raw and filtered movies. Also check whether
% the rawMovieName is not actually an imaris handle
if isnumeric(rawMovieName)
    % set movieLoader to none. Don't assign rawMovie.
    movieLoader = 'none';
elseif strcmp(class(rawMovieName),'COM.Imaris_Application')
    % load movie via imaris
    movieLoader = 'imaris';
    imarisHandle = rawMovieName;
elseif ~exist(rawMovieName,'file')
    error('rawMovie (''%s'') not found',rawMovieName)
else
    if isfield(dataProperties,'movieType') && strcmp(dataProperties.movieType,'sedat')
        movieLoader = 'sedat';
    else
        movieLoader = 'cdLoadMovie';
    end
end

% check for filteredMovie.
if isempty(filteredMovieName)
    error('please supply a filtered movie')
end
if isnumeric(filteredMovieName)
    % for now, filteredMovie can only be numeric if the rawMovie is
    % numeric, too.
    if strcmp(movieLoader,'none')
        % all is good.
    else
        error('It is not permitted to pass the filtered movie if the raw movie is passed by filename')
    end
elseif ~exist(filteredMovieName,'file')
    error('filteredMovie (''%s'') not found',filteredMovieName)
end


% check that dataProperties is a structure
if ~isstruct(dataProperties)
    error('dataProperties must be a structuer')
end

% assign debugData
debugData = [];

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
if ~strcmp(movieLoader,'none')
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
end

%=================================


%=================================
% FIND COORDINATES
%=================================

% preassign coordinates
coordinates(1:dataProperties.movieSize(4),1) = ...
    struct('sp',[],'COM',[]);



% find local maxima
switch dataProperties.detector_spotfind
    case 1 % standard spotfind

        % load filtered movie
        if strcmp(movieLoader,'none')
            filteredMovie = filteredMovieName;
            loadStruct.loadedFrames = 1:size(filteredMovie,5);
            loadStruct.frames2load = [];
            deltaFrames = 0;
        else
            % load filtered movie
            [filteredMovie, movieHeader, loadStruct] = ...
                cdLoadMovie({filteredMovieName,'filtered'}, [], loadOptions);
            deltaFrames = 0;
        end

        % loop through the movie to collect local maxima
        done = 0;
        while ~done

            coordinates(loadStruct.loadedFrames - deltaFrames) = ...
                spotfind(filteredMovie,dataProperties,verbose,loadStruct.loadedFrames);

            clear filteredMovie

            % load more
            if isempty(loadStruct.frames2load)
                done = 1;
            else

                [filteredMovie, movieHeader, loadStruct] = ...
                    cdLoadMovie(loadStruct.movieType, [], loadStruct);
            end % load more

        end % while loop find local max

    case 2
        % mammalian spotfind with changed cutoff
        % need to load movie within function
        coordinates = ...
            spotfind_mammalian(filteredMovieName,dataProperties,verbose,movieLoader);
end

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
    [dataProperties, testRatios, dbTmp] = ...
        detectSpots_MMF_findAmplitudeCutoff(...
        rawMovieName, coordinates, dataProperties, movieLoader, verbose, any(debug==2));
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
        if ~isfield(dataProperties,'crop')
            dataProperties.crop = [];
        end
        [rawMovie,movieSize,movieName,...
            moviePath,movieHeader,imarisHandle,loadStruct] = ...
            imarisImread(imarisHandle,[],dataProperties.crop,loadOptions.maxSize);
        deltaFrames = 0;

    case 'sedat'
        % load one frame
        rawMovie = sedatLoadRaw(1,rawMovieName,dataProperties);
        loadStruct.loadedFrames = 1;
        loadStruct.frames2load = 2:dataProperties.movieSize(4);
        deltaFrames = 0;
        % set numTimepoints in movieHeader
        movieHeader.numTimepoints = dataProperties.movieSize(4);

    case 'none'
        % read raw movie, set loadStruct and deltaFrames
        rawMovie = rawMovieName;
        loadStruct.loadedFrames = 1:size(rawMovie,5);
        loadStruct.frames2load = [];
        deltaFrames = 0;
        % need movieHeader for timepoints
        movieHeader.numTimepoints = size(rawMovie,5);
end

% preassign slist
slist(1:movieHeader.numTimepoints) = ...
    struct('sp',[],...
    'statistics',[],...
    'parms',[],...
    'COM',[]);

% take care of debug1
if any(debug) == 1
    debugData.fStats = cell(movieHeader.numTimepoints);
end


done = 0;
while ~done

    % prepare data for mixture model fitting
    loadedFrames = loadStruct.loadedFrames-deltaFrames;

    if ~isempty(testRatios)
        tr = testRatios(loadedFrames);
    else
        tr = [];
    end

    % do fitting
    [slist(loadedFrames), dbTmp] = ...
        detectSpots_MMF_main(rawMovie,coordinates(loadedFrames),...
        dataProperties,tr,verbose,debug);

    if debug == 1
        debugData.fStats(loadedFrames) = dbTmp.fStats;
    end

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

            case 'sedat'
                % load one more frame
                rawMovie = sedatLoadRaw(loadStruct.frames2load(1),rawMovieName,dataProperties);
                loadStruct.loadedFrames = loadStruct.frames2load(1);
                loadStruct.frames2load = loadStruct.frames2load(2:end);
        end

    end % load more

end % while loop find spots
