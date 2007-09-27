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
%                           "synth". MovieName should contain directory!
%
%           dirName (opt): directory name in which the movie can be found.
%                          If empty, current directory will be used
%           loadOpt (opt): Load options
%                            - array with a vector of timepoints to be
%                              loaded (e.g. [startFrame:endFrame])
%                            - -1: empty movie, will be returned, but
%                              complete movieHeader and loadStruct
%                            - structure with field .maxSize
%                              indicating maximum array size in bytes. If
%                              the movie is bigger than maxSize, only the
%                              first few frames will be loaded. With the
%                              information in loadStruct, cdLoadMovie can
%                              be called in a loop so that eventually, all
%                              frames of the movies will have been
%                              analyzed.
%                              * If .maxSize = 'check', cdLoadMovie will
%                              attempt to load dataProperties and find a
%                              field maxSize there to find maxSize. If
%                              there is no variable dataProperties (or
%                              tmpDataproperties), default dataProperties
%                              will be used; if there's no field .maxSize,
%                              the entire movie will be loaded.
%                              * If the structure contains an additonal field
%                              .frames2load, the program will load the
%                              appropriate number of frames from the list
%                              in .frames2load{1}.
%                              * If there is an additional field .noMovie
%                              that is 1, no movie will be loaded, but
%                              everything else will be returned (if this is
%                              the only option you want to set, supply -1
%                              for loadOpt instead; see above!)
%                              * .waveOrder, .waveIdx: If there are
%                              multiple wavelenghts, .waveOrder indicates
%                              how wavelenght is ordered compared to
%                              timepoints (see help r3dread). .waveIdx are
%                              the indices into the list of wavelenghts
%                              that should be loaded.
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
%                              .loadedFrames (output only): List of the
%                               frames that have just been loaded
%                              .waveOrder: [] if not otherwise specified
%                              .waveIdx: [] if not otherwise specified
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
goodTypes = {'raw';'corrected';'filtered';'latest';'corr/raw';'ask';'synth'};
knownType = length(goodTypes) + 1; % known movieName

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
    % read dirName
    type = knownType; % assign everything later
end

if nargin < 2 || isempty(dirName)
    dirName = pwd;
end

if nargin < 3
    loadOpt = [];
elseif isstruct(loadOpt) && isfield(loadOpt,'moviePath')
    dirName = loadOpt.moviePath;
end

% take care of noMovie
if isequal(loadOpt,-1)
    noMovie = 1;
    loadOpt = [];
else
    noMovie = 0;
end

if isstruct(loadOpt) && isfield(loadOpt, 'noMovie')
    noMovie = loadOpt.noMovie;
end

if isstruct(loadOpt) && ~isfield(loadOpt, 'waveIdx')
    loadOpt.waveIdx = [];
end

if isstruct(loadOpt) && ~isfield(loadOpt, 'waveOrder')
    loadOpt.waveOrder = [];
end

if isstruct(loadOpt) && ~isfield(loadOpt, 'crop')
    loadOpt.crop = [];
end
if isstruct(loadOpt) 
    if isempty(loadOpt.crop)
        doCrop = false;
    else
        doCrop = true;
    end
else
    doCrop = false;
end



%==========================


%==========================
% PREPARE LOADING
%==========================

% find movie file if necessary
if type == 6
    oldDir = cd(dirName);
    [movieName, dirName, chooseIdx] = ...
        uigetfile({'*.fim;moviedat*','filtered movie';...
        '*.r3d;*.r3c;*.dv','corrected movie';...
        '*.r3d;*.r3c;*.dv','raw movie'},...
        'Please choose movie type and location!');
    cd(oldDir);

    if movieName == 0
        % user aborted
        warning('CDLOADMOVIE:UserAbort','No movie loaded');
        if nargout > 0
            movie = 0;
            movieHeader = 0;
            loadStruct = 0;
        end
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
emptyHeader = false;

% if there is a loadStruct, we know the movie already. Otherwise, use type
% to decide on filename (*.r3d or *.fim)
if isstruct(loadOpt) && isfield(loadOpt,'movieName') ...
        && isfield(loadOpt,'frames2load')
    loadStruct = loadOpt;
    emptyHeader = true;
else
    % generate loadStruct
    loadStruct(1) = struct('movieName',[],...
        'frames2load',[],'correctionData',[], 'movieType',[],...
        'waveIdx',[],'waveOrder',[]);
    if type < knownType
        % get all files in the directory. If we know the file already, this
        % is a lot faster! (of course, it's even better if we don't have to
        % search at all)
        if isempty(movieName)
            allFileNames = dir(dirName);
            allFileNames(1:2) = [];
        else
            allFileNames = dir(fullfile(dirName, movieName));
        end

        % get list of filenames
        fileNameList = {allFileNames.name}';
        numFiles = length(allFileNames);
    end
    % find moviename. If type == 4, try until something is found
    if type == 3 || type == 4

        % find via regexpr
        regCell = regexp(fileNameList,'fim$|moviedat');
        % find index where there is something
        idx = find(~cellfun('isempty',regCell));

        % check whether there is a problem
        if isempty(idx)
            if type == 3
                error('no filtered movie found!')
            else
                % continue search below
            end

        else
            movieInfo = allFileNames(idx);

            if length(movieInfo) > 1
                warning('more than one filtered movie found. Loading the first')
                movieInfo = movieInfo(1);
            end

            % load movie header and assign correctionData
            r3dMovieHeader = loadMovieHeader(movieInfo.name);
            correctionData = [];

            % reset type 4
            type = 3;
        end
    end

    % if no filtered movie found, type 4 still exists. Continue searching
    if type == 1 || type == 2 || type == 4 || type == 5

        % find via regexpr - $ searches at end of name only - no problem
        % with .log
        regCell = regexp(fileNameList,'r3d$|3D.dv$');
        % find index where there is something
        idx = find(~cellfun('isempty',regCell));

        if isempty(idx)
            % we might not have found the movie, because it is a
            % simulation. Look for *r3c movie

            % search for raw movie file
            regCell = regexp(fileNameList,'r3c$');
            % find index where there is something
            idx = find(~cellfun('isempty',regCell));

            if isempty(idx)
                error('no movie found!')
            else
                % set everything for synth movie
                movieInfo = allFileNames(idx);
                load r3dMovieHeader
                correctionData = [];
                type = 7;
            end
        else
            % check whether there are multiple files. If yes, and one of
            % them is a DIC, then we take the other. If no DIC, we error.
            if length(idx) > 1
                regCell = regexpi(fileNameList(idx),'^DIC_|^pre_|^post_');
                % idxIdx is an index into idx
                idxIdx = find(~cellfun('isempty',regCell));
                if ~isempty(idxIdx) && (length(idxIdx) == length(idx) - 1)
                    % remove DIC-indices
                    idx(idxIdx) = [];
                else
                    error('multiple movies found in directory!')
                end
            end

            % if no problems, assign movieInfo and look for header and
            % correctionData
            movieInfo = allFileNames(idx);

            if type == 1 || ((type == 4 || type == 5) && ...
                    ~exist(fullfile(dirName,'correctionData.mat'),'file'))

                % read movie header and assign empty correctionData

                r3dMovieHeader = readr3dheader(movieInfo.name);
                % no correction data
                correctionData = [];

                % reset type 4
                type = 1;

            elseif type ~= 7

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
        movieName = movieType{1};
        movieInfo = dir(movieName); % store all info for part-loading

        % get correct directory
        tmpDirName  = fileparts(movieName);
        if ~isempty(tmpDirName)
            dirName = tmpDirName;
        end
        oldDir = cd(dirName);

        % load whatever necessary depending on movieData
        type = find(strcmpi(goodTypes,movieType{2}));

        if isempty(type)
            error('unrecognized movie type!')
        end

        switch type
            case 1 % raw - read header from movie itself!
                r3dMovieHeader = readr3dheader(movieInfo.name);
                correctionData = [];
            case 2 % corrected
                r3dMovieHeader = loadMovieHeader(movieInfo.name);
                load correctionData
            case 3 % filtered
                r3dMovieHeader = loadMovieHeader(movieInfo.name);
                correctionData = [];
            case 7 % synthetic - don't forget to adjust the type
                load r3dMovieHeader
                correctionData = [];

            case 5
                % if we weren't quite sure initially whether the movie was
                % corrected, case 5 occurs
                if exist(fullfile(dirName,'correctionData.mat'),'file')
                    % corrected movie
                    type = 2;
                    r3dMovieHeader = loadMovieHeader(movieInfo.name);
                    load correctionData
                else
                    % raw movie
                    type = 1;
                    r3dMovieHeader = readr3dheader(movieInfo.name);
                    correctionData = [];
                end
            otherwise
                error('current type cannot be handled!')
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


    % take care of wave-options. Remove fields from loadOpt to make sure
    % that we don't accidently get in trouble below with a missing field
    % "maxSize"
    if isstruct(loadOpt)
        loadStruct.waveIdx = loadOpt.waveIdx;
        loadStruct.waveOrder = loadOpt.waveOrder;
        loadOpt = rmfield(loadOpt,'waveIdx');
        loadOpt = rmfield(loadOpt,'waveOrder');
        loadStruct.crop = loadOpt.crop;
        loadOpt = rmfield(loadOpt,'crop');

        % check whether there are any fields left for loadOpt
        if isempty(fieldnames(loadOpt))
            loadOpt = [];
        end
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

        if ischar(loadOpt.maxSize) && strcmp(loadOpt.maxSize,'check')
            % go and look for dataProperties.maxSize
            if ~isempty(dir('dataProperties.mat'))
                load dataProperties
                if isfield(dataProperties,'maxSize')
                    loadOpt.maxSize = dataProperties.maxSize;
                end
            elseif ~isempty(dir('tmpDataProperties.mat'))
                load tmpDataProperties
                if isfield(dataProperties,'maxSize')
                    loadOpt.maxSize = dataProperties.maxSize;
                end
            else
                dpName = dir('*dataProperties*.mat');
                if ~isempty(dpName)
                    load(dpName)
                    if isfield(dataProperties,'maxSize')
                        loadOpt.maxSize = dataProperties.maxSize;
                    end
                end

            end
            % if we still don't know the size, go for defaultDataProperties
            if strcmp(loadOpt.maxSize,'check')
                dataProperties = defaultDataProperties;
                loadOpt.maxSize = dataProperties.maxSize;
            end
        end
        maxSize = loadOpt.maxSize;
        movieSize = movieInfo.bytes;

        % if it is a raw movie, the data is (mostly) uint16. The file is
        % roughly 4 times smaller than doubles - correct
        if type == 1 || type == 2
            movieSize = movieSize * 4;
        end
        % if we're not loading all the wavelengths, the movieSize gets
        % effectively smaller
        if ~isempty(loadStruct.waveIdx)
            movieSize = movieSize * ...
                (length(loadStruct.waveIdx)/r3dMovieHeader.numWvs);
        end
        % if we're cropping, movieSize gets even smaller
        if doCrop
            isCrop = any(loadStruct.crop(:,1:3),1);
            if any(isCrop)
                croppedSize = diff(loadStruct.crop(:,1:3),1,1)+1;
                frameSize = [r3dMovieHeader.numRows,...
                    r3dMovieHeader.numCols, r3dMovieHeader.numZSlices];
                croppedSize(1,~isCrop) = frameSize(1,~isCrop);
                movieSize = movieSize * ...
                    prod(croppedSize)/prod(frameSize);
                % if we don't crop along one of the dimensions, we say we take
                % from 1:end instead 
                loadStruct.crop(1,~isCrop) = 1;
                loadStruct.crop(2,~isCrop) = frameSize(1,~isCrop); 
                loadStruct.crop = loadStruct.crop(:,1:3);
            else
                % don't crop if not necessary
                doCrop = false;
            end
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

        if min(loadOpt) < 1  || max(loadOpt) > r3dMovieHeader.numTimepoints
            error('Start frame or end frame out of range!')
        end


        % add list of frames
        loadStruct.frames2load{1} = loadOpt +...
            correctFrame;

    else % no loadOpt, therefore, we want to load all

        loadStruct.frames2load{1} = (1:r3dMovieHeader.numTimepoints) + ...
            correctFrame;

    end



    % fill in loadStruct
    loadStruct.movieName = movieInfo.name;
    loadStruct.movieType = goodTypes{type};
    loadStruct.moviePath = dirName;
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
            if isstruct(loadOpt) && isfield(loadOpt,'maxSize')
                % account for the possibility that we have size limits
                testStruct.maxSize = loadOpt.maxSize;
            else
                testStruct = 1:r3dMovieHeader.numTimepoints;
            end

            % add wave options
            testStruct.waveIdx = loadStruct.waveIdx;
            testStruct.waveOrder = loadStruct.waveOrder;

            % check whether this works at all (at the moment there is no
            % possibility of using dark frame subtraction with multiple
            % wavelengths
            if length(testStruct.waveIdx) > 1
                error('multicolor darkframes not implemented yet')
            end

            % load movie
            [testMovie,dummy,testLoadStruct] = cdLoadMovie('raw',[],testStruct);

            % subtract background. The computer has to be able to handle
            % two variables of maxSize
            testMovie = testMovie - ...
                repmat(loadStruct.correctionData.image,...
                [1,1,1,1,size(testMovie,5)]);

            minimumIntensity = min(testMovie(:));

            % loop in case we haven't loaded everything yet
            while ~isempty(testLoadStruct.frames2load)
                [testMovie,dummy,testLoadStruct] = ...
                    cdLoadMovie('raw',[],testLoadStruct);
                testMovie = testMovie - ...
                    repmat(correctionData.image,...
                    [1,1,...
                    size(testMovie,3)/size(correctionData.image,3),...
                    size(testMovie,4),size(testMovie,5)]);
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

% load only if necessary
if noMovie
    % don't load
    movie = [];
else

    % first: check if the loadList is continuous (for type 1,2).
    loadList = loadStruct.frames2load{1};
    loadListLength = length(loadList);
    dll = diff(loadList);
    if (isempty(dll) || all(dll==1)) && ~doCrop
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
                    loadListLength,[],loadStruct.waveIdx,loadStruct.waveOrder);
            else
                % load movie timepoint by timepoint
                for t = loadListLength:-1:1
                    tmp = r3dread(loadStruct.movieName,...
                        loadList(t),1,[],loadStruct.waveIdx,loadStruct.waveOrder);
                    if doCrop
                        movie(:,:,:,:,t) = tmp(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                            loadStruct.crop(1,2):loadStruct.crop(2,2),...
                            loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                    else
                        movie(:,:,:,:,t) = tmp;
                    end

                end
            end

        case 2
            if isContinuous
                % load all at once
                % r3dread: filename, start, numFrames
                movie = r3dread(loadStruct.movieName,loadList(1), ...
                    loadListLength,[],loadStruct.waveIdx,loadStruct.waveOrder);
            else
                % load movie timepoint by timepoint
                for t = loadListLength:-1:1
                    tmp = r3dread(loadStruct.movieName,...
                        loadList(t),1,[],loadStruct.waveIdx,loadStruct.waveOrder);
                    if doCrop
                        movie(:,:,:,:,t) = tmp(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                            loadStruct.crop(1,2):loadStruct.crop(2,2),...
                            loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                    else
                        movie(:,:,:,:,t) = tmp;
                    end
                end
            end

            % subtract background. If the computer cannot cope with two
            % matrices of maxSize, all is lost, anyway.
            if doCrop
                movie = movie - ...
                    repmat(loadStruct.correctionData.image(...
                    loadStruct.crop(1,1):loadStruct.crop(2,1),...
                    loadStruct.crop(1,2):loadStruct.crop(2,2),...
                    loadStruct.crop(1,3):loadStruct.crop(2,3)),...
                    [1,1,1,1,size(movie,5)]) - ...
                    loadStruct.correctionData.minimumIntensity;

                movie = movie - ...
                    repmat(loadStruct.correctionData.image,...
                    [1,1,1,1,size(movie,5)]) - ...
                    loadStruct.correctionData.minimumIntensity;
            end



        case 3
            % readmat - check for single-frame movie
            if ~emptyHeader && r3dMovieHeader.numTimepoints == 1
                movie = readmat(loadStruct.movieName);
                if doCrop
                    movie = movie(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                        loadStruct.crop(1,2):loadStruct.crop(2,2),...
                        loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                end
            elseif doCrop
                for t = loadListLength:-1:1
                    tmp = readmat(loadStruct.movieName,loadList(t));

                    if all(diff(loadStruct.crop,1,1)+1 == size(tmp,1:3))
                        % already cropped
                        movie(:,:,:,:,t) = tmp;
                    else

                        movie(:,:,:,:,t) = tmp(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                            loadStruct.crop(1,2):loadStruct.crop(2,2),...
                            loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                    end

                end
            else
                movie = readmat(loadStruct.movieName,loadList);
            end

        case 7
            % readmat - check for single-frame movie
            if r3dMovieHeader.numTimepoints == 1
                movie = readmat(loadStruct.movieName);
                if doCrop
                    movie = movie(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                        loadStruct.crop(1,2):loadStruct.crop(2,2),...
                        loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                end
            elseif doCrop
                for t = loadListLength:-1:1
                    tmp = readmat(loadStruct.movieName,loadList(t));
                    if all(diff(loadStruct.crop,1,1)+1 == size(tmp,1:3))
                        % already cropped
                        movie(:,:,:,:,t) = tmp;
                    else

                        movie(:,:,:,:,t) = tmp(loadStruct.crop(1,1):loadStruct.crop(2,1),...
                            loadStruct.crop(1,2):loadStruct.crop(2,2),...
                            loadStruct.crop(1,3):loadStruct.crop(2,3),:,:);
                    end
                end
            else
                movie = readmat(loadStruct.movieName,loadList);
            end

        otherwise
            error('type %i not handled!',type)
    end

end


%============================
% ASSIGN OUTPUT
%============================

% move the first entry of files2load
% into the list of loaded frames
if noMovie
    loadStruct.loadedFrames = [];
else
    loadStruct.loadedFrames = loadStruct.frames2load{1};
    loadStruct.frames2load(1) = [];
end

% rename header
if emptyHeader
    movieHeader = [];
else
    movieHeader = r3dMovieHeader;
end

% go back to old dir
cd(oldDir);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   SUBFUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loadMovieHeader
function movieHeader = loadMovieHeader(movieName)

% this subfunction was written to be able to handle "misteli"-movies that
% are multi-color metamorph stacks. Since there are two movieHeaders in a
% directory, we need to be able to load the correct one.

% movieHeader is the movieHeader corresponding to movieName

% find movieHeaders in the current directory. Be careful for
% capitalizations!
movieHeaderNames = dir('*ovieHeader*');

% switch on the number of movieHeaders we find: if there is none, we go and
% look for an r3d-file (if unsucessful, too, we die).
% If there's one header: load it.
% If there are multiple headers: There will be an identifier between
% "Header" and ".mat". Find this part for every header and doublecheck with
% the movieName.
nMovieHeaders = length(movieHeaderNames);
switch nMovieHeaders
    case 0 % no r3dMovieHeader. Read from original file
        try
            rawMovie = dir('*.r3d');
            movieHeader = readr3dheader(rawMovie.name);
        catch
            error('No movie header found for %s in %s',movieName,pwd)
        end
    case 1 % only one header. If only if was always that easy
        % just to be sure: allow movieHeaders named any way
        file = load(movieHeaderNames.name);
        fn = fieldnames(file);
        movieHeader = file.(fn{1});
    otherwise % multiple headers.
        % find identifier in filenames
        [movieHeaderNameList{1:nMovieHeaders}] = deal(movieHeaderNames.name);
        identifiers = ...
            regexpi(movieHeaderNameList,'\w*header(.*).mat','tokens');
        % since the output is very inconveniently packed into cells, we
        % loop here until we find the correct header
        found = 0;
        headerIdx = 0;
        while ~found && headerIdx < nMovieHeaders
            headerIdx = headerIdx + 1;
            currentID = identifiers{headerIdx};
            currentID = currentID{1};

            % find currentID in movieName
            if isempty(regexpi(movieName,currentID))
                % continue loop
            else
                found = 1;
            end
        end

        % check if we found anything
        if ~found
            error('No movie header found for %s in %s',movieName,pwd)
        else
            file = load(movieHeaderNames(headerIdx).name);
            fn = fieldnames(file);
            movieHeader = file.(fn{1});
        end

end % switch on number of headers found