function [movie, movieHeader, loadStruct] = cdLoadMovie(movieType, dirName, loadOpt)
%CDLOADMOVIE loads movies in the Chromdyn project
%
% SYNOPSIS [movie, movieHeader, loadStruct] = cdLoadMovie(movieType, dirName, loadOpt)
%
% INPUT     movieType    : "raw", "corrected", or "filtered". "latest" will
%                           try to load filtered/corrected/raw movie in
%                           this order, if it exists.
%                           "corr/raw" will do the same, but without trying
%                           to load the filtered movie.
%                           "ask" will let the user choose the movie
%                           For synthetic movies, specify "raw"
%                           
%                           Alternatively, if the filename is known,
%                           movieType can be a 1-by-2 cell array with 
%                           {movieName, movieType}, where movieType is
%                           either "raw", "corrected", "filtered", or
%                           "synth"
%
%           dirName (opt): directory name in which the movie can be found.
%                          If empty, current directory will be used
%           loadOpt (opt): Load options
%                            - array with a vector of timepoints to be
%                              loaded (e.g. [startFrame:endFrame])
%                            - jobProperties structure with field .maxSize
%                              indicating maximum array size in bytes. If
%                              the movie is bigger than maxSize, only the
%                              first few frames will be loaded. With the
%                              information in loadStruct, cdLoadMovie can
%                              be called in a loop so that eventually, all
%                              frames of the movies will have been
%                              analyzed.
%                              If the structure contains an additonal field
%                              .frames2load, the program will load the
%                              appropriate number of frames from the list
%                              in .frames2load{1}.
%                            - loadStruct structure with fields
%                              .frames2load, a cell array containing the
%                               vectors with the timepoints to load. If
%                               no further timepoints are to be loaded, the
%                               cell will be empty and can be used as the
%                               exit condition in the calling loop.
%                              .correctionData: data necessary for
%                               correction
%                              .movieName: filename to load
%                              .movieType: if movieType was "latest" or
%                               "corr/raw": type which was chosen.
%
% OUTPUT    movie        : 5D movie file (x,y,z,wavelenghth, time)
%           movieHeader  : header of movie
%           loadStruct   : see description of loadOpt(3)
%
% c: 4/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================
% TEST INPUT
%=========================
goodTypes = {'raw';'corrected';'filtered';'latest';'corr/raw';'ask'};
synthType = length(goodTypes) + 1; % *.r3c-movie. Still used for synth movie
knownType = length(goodTypes) + 2; % known movieName

if nargin < 1 || isempty(movieType) || ~(ischar(movieType) || iscell(movieType))
    error('cdLoadMovie needs a movieType as first input argument!')
end

if ischar(movieType)
    % find movieType
    type = find(strcmpi(goodTypes,movieType));
    if isempty(type)
        error('movieType not recognized. Possible types are: ''raw'', ''corrected'', ''filtered'', ''latest'', ''corr/raw'', or ''ask''!');
    end
else % movieType is a cell
    if length(movieType) ~= 2
        error('if known movieName: Please supply a cell array with {movieName, movieType}')
    end
    movieType = knownType; % assign everything later
end

if nargin < 2 || isempty(dirName)
    dirName = pwd;
end

if nargin < 3
    loadOpt = [];
end

%==========================


%==========================
% PREPARE LOADING
%==========================

% find movie file if necessary
if type == 6
    [movieName, dirName, chooseIdx] = ...
        uigetfile({'*.fim;moviedat*','filtered movie';...
        '*.r3d;*.r3c','corrected movie';...
        '*.r3d;*.r3c','raw movie'},...
        'Please choose movie type and location!');
    if movieName == 0
        % user aborted
        warning('No movie loaded');
        movie = 0;
        movieHeader = 0;
        loadStruct = 0;
        return
    else
        % translate type
        type = 4 - chooseIdx;
    end
else
    % we want to use movieName below - assing empty here
    movieName = [];
end




% goto movie directory
oldDir = cd(dirName);

% if we're in a loop, there is no point in loading the header all the time.
% Therefore, we return it empty
emptyHeader = 0;

% if there is a loadStruct, we know the movie already. Otherwise, use type
% to decide on filename (*.r3d or *.fim)
if isstruct(loadOpt) && isfield(loadOpt,'movieName')
    loadStruct = loadOpt;
    emptyHeader = 1;
else
    % generate loadStruct
    loadStruct(1) = struct('movieName',[],'frames2load',[],'correctionData',[], 'movieType',[]);

    % get all files in the directory. If we know the file already, this
    % is a lot faster!
    if isempty(movieName)
        allFileNames = dir(dirName);
        allFileNames(1:2) = [];
    else
        allFileNames = dir([dirName movieName]);
    end

    % get list of filenames
    fileNameList = {allFileNames.name}';
    numFiles = length(allFileNames);

    % find moviename. If type == 4, try until something is found
    if type == 3 || type == 4
        i=1;
        % search for filtered movie file
        while numFiles >= i && ...
                isempty(findstr(fileNameList{i},'.fim')) && ...
                isempty(findstr(fileNameList{i},'moviedat'))
            i = i + 1;
        end

        if i > numFiles
            if type == 3
                error('no filtered movie found!')
            else
                % continue search below
            end

        else % assign movieInfo etc
            movieInfo = allFileNames(i);

            % load movie header and assign correctionData
            load r3dMovieHeader
            correctionData = [];

            % reset type 4
            type = 3;
        end
    end

    % if no filtered movie found, type 4 still exists. Continue searching
    if type == 1 || type == 2 || type == 4 || type == 5
        i=1;
        % search for raw movie file
        while numFiles >= i && ...
                isempty(findstr(fileNameList{i},'.r3d')) || ...
                ~isempty(findstr(fileNameList{i},'.log'))
            i = i + 1;
        end

        if i > numFiles
            % we might not have found the movie, because it is a
            % simulation. Look for *r3c movie
            i=1;
            % search for raw movie file
            while numFiles >= i && ...
                    isempty(findstr(fileNameList{i},'.r3c')) || ...
                    ~isempty(findstr(fileNameList{i},'.log'))
                i = i + 1;
            end
            if i > numFiles
                error('no movie found!')
            else
                % set everything for synth movie
                movieInfo = allFileNames{i};
                load r3dMovieHeader
                correctionData = [];
                type = synthType;
            end
        else

            % if no problems, assign movieInfo and look for header and
            % correctionData
            movieInfo = allFileNames(i);

            if type == 1 || ((type == 4 || type == 5) && ...
                    ~exist('correctionData.mat','file'))

                % read movie header and assign empty correctionData

                r3dMovieHeader = readr3dheader(movieInfo.name);
                % no correction data
                correctionData = [];

                % reset type 4
                type = 1;

            elseif type ~= synthType

                % load movie header and correctionData
                load r3dMovieHeader
                load correctionData

                % reset type 4, 5
                type = 2;
            else
                % synthMovie is already assigned
            end
        end
    end % type == 1 || type == 2 || type == 4 || type == 5

    % handle known moviename
    if type == knownType
        % we just assume that we actually know what the movie is called
        movieInfo = dir(movieType{1}); % store all info for part-loading
        
        % load whatever necessary depending on movieData
        type = find(strcmpi([goodTypes(1:3);{'synth'}]));
        
        if isempty(type)
            error('unrecognized movie type!')
        end
        
        switch type
            case 1 % raw
                r3dMovieHeader = readr3dheader(movieInfo.name);
                correctionData = [];
            case 2 % corrected
                load r3dMovieHeader
                load correctionData
            case 3 % filtered
                load r3dMovieHeader
                correctionData = [];
            case 4 % synthetic - don't forget to adjust the type
                type = synthType;
                load r3dMovieHeader
                correctionData = [];
        end
    end % if type == knownType
        
        
    
    
    % if this is a corrected movie, and it is corrected via the first/last few
    % frames, we have to adjust frames2load (no problem with numTimepoints,
    % though - they are adjusted already)

    if ~isempty(correctionData) && ~isempty(correctionData.info.correctFrames)
        correctFrame = correctionData.info.correctFrames(1);
        movieCorrection = sum(correctionData.info.correctFrames);
    else
        correctFrame = 0;
        movieCorrection = 0;
    end


    % if jobPropertiesStruct: find size of movie, calculate minimum number of
    % divisions. Then generate loadStruct with frames2load
    % else: write startFrame/endFrame into loadStruct.frames2load
    % error if number of frames doesn't work out.
    if isstruct(loadOpt)
        % test if good structure
        if ~isfield(loadOpt,'maxSize')
            error('loadOpt-structure has no known fields!')
        end

        maxSize = loadOpt.maxSize;
        movieSize = movieInfo.bytes;

        % if it is a raw movie, the data is (mostly) uint16. The file is
        % roughly 4 times smaller than doubles - correct
        if type == 1 || type == 2
            movieSize = movieSize * 4;
        end

        movieLength = r3dMovieHeader.numTimepoints;

        % check whether there is a list of frames to load. If yes, take
        % into account the fact that the total movie will be shorter.
        if isfield(loadOpt,'frames2load')
            % I'm not testing if it's a cell here - the user will hopefully
            % be able to debug if necessary
            frameList = loadOpt.frames2load{1};
        else
            frameList = 1:movieLength;
        end

        frameListLength = length(frameList);
        sizeFactor = frameListLength/(movieLength + movieCorrection);

        % number of divisions: movieSize/maxSize
        numParts = ((sizeFactor*movieSize)/maxSize);

        if numParts > frameListLength
            warning('One single frame is bigger than maxSize. Attempting to load single frames.')
            numParts = frameListLength;
        end

        % generate frames2load. Round to make sure that we do not get too
        % big. It is well possible that there is one very long and one very
        % short part of the movie. The speed increase with having two
        % similarly sized movies is so minimal, though, that I do not try
        % to make the code more beautiful
        numFrames = floor(frameListLength/numParts);
        numParts = ceil(frameListLength/numFrames);


        % loop to build cell array with frames2load. Every entry is a list
        % of frames. Start counting at the end to prevent constant
        % reassignement of array.
        for i=numParts:-1:1
            loadStruct.frames2load{i,1} = ...
                frameList((i-1)*numFrames+1:min(i*numFrames,frameListLength))...
                + correctFrame;
        end

    elseif isnumeric(loadOpt) && ~isempty(loadOpt)

        if min(loadOpt) < 1  | max(loadOpt) > r3dMovieHeader.numTimepoints
            error('Start frame or end frame out of range!')
        end


        % add list of frames
        loadStruct.frames2load{1} = loadOpt +...
            correctFrame;

    else % no loadOpt, therefore, we want to load all

        loadStruct.frames2load{1} = [1:r3dMovieHeader.numTimepoints] + ...
            correctFrame;

    end



    % fill in loadStruct
    loadStruct.movieName = movieInfo.name;
    loadStruct.movieType = goodTypes{type};
    if ~isempty(correctionData)
        % check if we have stored an image or just a slice
        if ndims(correctionData.image)==2
            loadStruct.correctionData.image = repmat(correctionData.image,[1,1,r3dMovieHeader.numZSlices]);
        else
            loadStruct.correctionData.image = correctionData.image;
        end
    else
        loadStruct.correctionData = [];
    end

    % ok, so no it gets a little weird: For corrected movies, we need to
    % make sure that the minimum of the movie never becomes negative.
    % Therefore, the first time we load a corrected movie, we need to find
    % its global minimum, so that we can save it in the correctionData
    if ~isempty(correctionData)
        if ~isfield(correctionData,'minimumIntensity')

            % we need to load the entire movie (possibly in parts),
            % subtract the background and then find the minimum
            if isstruct(loadOpt) && isfield(loadOpt.maxSize)
                % account for the possibility that we have size limits
                testStruct.maxSize = loadOpt.maxSize;
            else
                testStruct = [1:r3dMovieHeader.numTimepoints];
            end
            
            % load movie
            [testMovie,dummy,testLoadStruct] = cdLoadMovie('raw',[],testStruct);
            
            % subtract background. The computer has to be able to handle
            % two variables of maxSize
            testMovie = testMovie - ...
                repmat(loadStruct.correctionData.image,...
                [1,1,1,size(testMovie,4),size(testMovie,5)]);

            minimumIntensity = min(testMovie(:));
            
            % loop in case we haven't loaded everything yet
            while ~isempty(testLoadStruct.frames2load)
                [testMovie,dummy,testLoadStruct] = ...
                    cdLoadMovie('raw',[],testLoadStruct);
                testMovie = testMovie - ...
                    repmat(correctionData.image,...
                    [1,1,1,size(testMovie,4),size(testMovie,5)]);
                minimumIntensity = min(min(testMovie(:)),minimumIntensity);
            end
            
            % assign minimumIntensity and save correctionData
            loadStruct.correctionData.minimumIntensity = ...
                minimumIntensity;
            correctionData.minimumIntensity = minimumIntensity;
            save('correctionData','correctionData')
            
            % try to free some space
            clear testMovie

        else
            % just assign minimumIntensity
            loadStruct.correctionData.minimumIntensity = ...
                correctionData.minimumIntensity;

        end
    end

end %if isstruct(loadOpt) && isfield(loadOpt,'movieName')




%===============================
% LOAD MOVIE
%===============================

% first: check if the loadList is continuous (for type 1,2).
loadList = loadStruct.frames2load{1};
loadListLength = length(loadList);
dll = diff(loadList);
if isempty(dll) || all(dll==1)
    isContinuous = 1;
else
    isContinuous = 0;
end

switch type
    case 1
        if isContinuous
            % load all at once
            % r3dread: filename, start, numFrames
            movie = r3dread(loadStruct.movieName,loadList(1), ...
                loadListLength);
        else
            % load movie timepoint by timepoint
            for t = loadListLength:-1:1
                movie(:,:,:,:,t) = r3dread(loadStruct.movieName,...
                    loadList(t),1);
            end
        end

    case 2
        if isContinuous
            % load all at once
            % r3dread: filename, start, numFrames
            movie = r3dread(loadStruct.movieName,loadList(1), ...
                loadListLength);
        else
            % load movie timepoint by timepoint
            for t = loadListLength:-1:1
                movie(:,:,:,:,t) = r3dread(loadStruct.movieName,...
                    loadList(t),1);
            end
        end

        % subtract background. If the computer cannot cope with two
        % matrices of maxSize, all is lost, anyway.
        movie = movie - ...
            repmat(loadStruct.correctionData.image,...
            [1,1,1,size(movie,4),size(movie,5)]) - ...
            loadStruct.correctionData.minimumIntensity;



    case 3
        % readmat
        movie = readmat(loadStruct.movieName,loadList);

    case synthType
        % readmat
        movie = readmat(loadStruct.movieName,loadList);
end


%============================
% ASSIGN OUTPUT
%============================

% remove the first entry of files2load
loadStruct.frames2load(1) = [];

% rename header
if emptyHeader
    movieHeader = [];
else
    movieHeader = r3dMovieHeader;
end

% go back to old dir
cd(oldDir);