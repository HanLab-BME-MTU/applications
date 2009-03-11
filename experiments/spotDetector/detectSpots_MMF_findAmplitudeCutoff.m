function [dataProperties, testRatios, debugData] = detectSpots_MMF_findAmplitudeCutoff(rawMovieName, coordinates, dataProperties,movieLoader,verbose,debug)
%DETECTSPOTS_MMF_FINDAMPLITUDECUTOFF finds the cutoff that distinguishes the amplitude of true spots from noise
%
% SYNOPSIS: [dataProperties, testRatios] = detectSpots_MMF_findAmplitudeCutoff(rawMovieName, coordinates, dataProperties,movieLoader,verbose,debug)
%
% INPUT rawMovieName: full name of raw movie file or imaris handle or
%           imageDataObj
%		coordinates: coordinates of local maxima from spotfind
%		dataProperties: dataProperties structure
%		movieLoader: (opt) function for loading the movie
%               ['imaris'/{'cdLoadMovie'}/'sedat'/'obj']
%       verbose: (opt)  0: nothing at all.
%                      {1} show waitbar.
%                       2: show waitbar and cutoff-plot
%
% OUTPUT dataProperties: dataProperties with filled field "amplitudeCutoff"
%        testRatios:     nTimepoints-by-1 cell array with
%                        {timepoint, testValue} for the amplitude fits
%                        corresponding to the coordinates
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: Jonas Dorn
% DATE: 07-Feb-2006
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input - not really done, b/c this is merely a helper function

% assign default movieLoader
if nargin < 4 || isempty(movieLoader)
    movieLoader = 'cdLoadMovie';
end

if nargin < 5 || isempty(verbose)
    verbose = 1;
end

if nargin < 6 || isempty(debug)
    debug = 0;
end
% assign movieLoader in case it's a movieDataObj
if ~isempty(rawMovieName) && isa(rawMovieName,'imageDataObj')
    movieLoader = 'obj';
end

% find maximum allowed data size
if isfield(dataProperties,'maxSize')
    loadOptions.maxSize = dataProperties.maxSize;
else
    % set maximum size to 100 Mb
    loadOptions.maxSize = 100000000;
end

% turn off singularMatrix warning
warningState = warning;
warning off MATLAB:nearlySingularMatrix

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
            imarisImread(rawMovieName,[],dataProperties.crop,loadOptions.maxSize);
        deltaFrames = 0;
    case 'sedat'
        % load one frame
        rawMovie = sedatLoadRaw(1,rawMovieName,dataProperties);
        loadStruct.loadedFrames = 1;
        loadStruct.frames2load = 2:dataProperties.movieSize(4);
        deltaFrames = 0;
    case 'none'
        % entire movie is being passed down as 'rawMovieName'
        rawMovie = rawMovieName;
        loadStruct.loadedFrames = 1:size(rawMovie,5);
        loadStruct.frames2load = [];
        deltaFrames = 0;
    case 'obj'
        rawMovie = rawMovieName;
        % 'loadedFrames' are from minT to maxT
        loadStruct.loadedFrames = rawMovie.cropInfo.crop(1,4):...
            rawMovie.cropInfo.crop(2,4);
        % frames2load is empty, because we are going to go through the
        % entire movie in one loop
        loadStruct.frames2load = [];
        deltaFrames = 0;
    otherwise
        error('movieLoader %s not recognized',movieLoader)
end

debugData = [];

%===================
% FIND CUTOFF
%===================

% find cutoff for amplitude. Make it the easy way first (we can try to
% program it faster later): For all frames, fit the local maxima (just N
% spots) and store amp/sigmaAmp. Make a histogram of that to find the best
% cut-off for good spots vs. noise. Then do the real fitting.


% define store variable for testRatios
nTimepoints = dataProperties.movieSize(4);
testRatios = cell(nTimepoints,1);
if debug
    debugData.residualImages = cell(nTimepoints*(dataProperties.MAXSPOTS+3),1);
    debugData.resAndGauss = cell(nTimepoints*(dataProperties.MAXSPOTS+3),1);
    imCt = 1;
end



if verbose
    h= mywaitbar(0,[],nTimepoints,'Finding amplitude cutoff...');
end


% loop and load movieChunks
done = 0;
maxNumSpots = 9999;
maxNumCoords = 9999;
while ~done
    
    % loop t to fit local maxima
    for t = 1:length(loadStruct.loadedFrames)
        currentTime = loadStruct.loadedFrames(t) - deltaFrames;
        % find spots to fit
        if isfield(coordinates(currentTime),'sp') &&...
                ~isempty(coordinates(currentTime).sp)
            %
            % have try-catch loop to watch out for memory issues
            try
                
                cordList=cat(1,coordinates(currentTime).sp.cord);
                %change to matlab coords
                tc2=cordList(:,2);
                cordList(:,2)=cordList(:,1);
                cordList(:,1)=tc2;
                
                % if we have an imageDataObj, we get the frame via getFrame
                % directly so that there won't be memory issues
                if strcmp(movieLoader,'obj')
                    imgStk = rawMovie.getFrame(currentTime);
                else
                    imgStk=rawMovie(:,:,:,1,t);
                end
                
                % preassign current testRatios
                nSpots = size(cordList,1);
                testRatios{currentTime} = zeros(nSpots,2);
                spotIdxList = 1:nSpots;
                
                while (~isempty(cordList))
                    discerned = 0;
                    if length(cordList) >= maxNumCoords
                        error('MATLAB:nomem',...
                            'too many coords (%i) - Matlab will run out of memory!',...
                            length(cordList));
                    end
                    [spotsidx, mask] = discernspots(cordList,size(imgStk),dataProperties);
                    % remember which spots we're currently fitting!
                    currentSpots = spotIdxList(spotsidx);
                    discerned = 1;
                    if length(spotsidx) >= maxNumSpots
                        error('MATLAB:nomem',...
                            'too many spots (%i) - Matlab will run out of memory!',...
                            length(spotsidx));
                    end
                    idxList=find(mask);
                    mskData=imgStk(idxList);
                    
                    % if the spot is at the very edge of a cropped region,
                    % there may be NaNs in the mask. As edge spots are not
                    % interesting, anyway, discard the spot
                    
                    if any(isnan(mskData))
                        
                        % do nothing
                    else
                        
                        try
                            
                            % do the maximum fitting
                            if debug
                                [testValue debugData.residualImages{imCt},debugData.resAndGauss{imCt}]=...
                                    fitTestFirstRound(mskData,cordList(spotsidx,:),idxList,...
                                    size(imgStk),dataProperties);
                                imCt = imCt+1;
                            else
                                testValue =...
                                    fitTestFirstRound(mskData,cordList(spotsidx,:),idxList,...
                                    size(imgStk),dataProperties);
                            end
                            
                            % collect testValues. Be careful:
                            % - discernspots picks groups that can be any combination
                            % - fitTest reverses the order of the spots
                            currentRatios = testRatios{currentTime};
                            currentRatios(currentSpots(end:-1:1),:) = ...
                                [ones(length(testValue),1)*currentTime,testValue(:)];
                            testRatios{currentTime} = currentRatios;
                            
                        catch err
                            if strcmp(err.identifier,'optim:snls:InvalidUserFunction')
                                % do nothing. Simply move on and remove
                                % spots
                            else
                                % probably an out of memory error. Rethrow
                                % and catch later
                                rethrow(err)
                            end
                        end % try and fail gently in case of nan in msk
                    end % check for nan in msk
                    
                    % remove spots from cordList - whether we fitted or
                    % discarded
                    cordList(spotsidx,:)=[];
                    spotIdxList(spotsidx) = [];
                    
                end
            catch
                % check the type of error
                err = lasterror;
                if strmatch(err.identifier,'MATLAB:nomem')
                    % we have a memory problem
                    % remember length of spotsIdx
                    if discerned
                        maxNumSpots = min(maxNumSpots,length(spotsidx));
                    else
                        maxNumCoords = min(maxNumCoords,length(cordList));
                    end
                    
                    % remove the spots that have already been analyzed
                    testRatios{t} = [];
                    
                    disp(sprintf(...
                        'MMF_findAmplitude aborted in frame %i:\n%s',...
                        currentTime,err.message))
                else
                    % this is not a memory problem. rethrow the error
                    rethrow(err)
                end
            end
        end % if ~isempty
        if verbose
            mywaitbar(currentTime/nTimepoints,h,nTimepoints);
        end
    end % loop time
    
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
end % while loop

%==================================
% GET CUTOFF
%==================================


% get cutoff. Transform cell to matrix and read out ratios only
tmp = cat(1,testRatios{:});
ratios = tmp(:,2);
times = tmp(:,1);
clear tmp

% if there's only one timepoint, there'll be trouble. Just take maxspots
% points
if nTimepoints < 2
    % get testRatios
    ratios = -sort(-ratios);
    nRatios = length(ratios);
    % cut at MAXSPOTS
    cutIdx = min(nRatios,dataProperties.MAXSPOTS);
    cutValue = 0.99999 * ratios(cutIdx);
    
    
else
    
    % remove very small values from ratios. If there are multiple almost zeros,
    % for example, there will be a peak at 0 that will screw up
    % cutFirstHistMode
    badIdx = ratios < 1e-3;
    if sum(badIdx) > 0
        %disp(sprintf('findAmpCutoff: %i low int spots',sum(badIdx)))
        times(badIdx) = [];
        ratios(badIdx) = [];
    end
    
    % find cutoff
    
    % first guess via cutFirstHistMode
    [cutIdx, cutVal,sp] = cutFirstHistMode(ratios,0);
    
    % now check the local minima in the vicinity of the cutoff
    spder = fnder(sp);
    zeroList = fnzeros(spder);
    zeroList = zeroList(1,:);
    % evaluate
    zeroVals = fnval(sp,zeroList);
    
    % look in zeroList. Find one value before cutVal, three after. Go into
    % zeroVals and find lowest minimum
    [dummy,closestIdx] = min(abs(zeroList - cutVal));
    
    % check only the minimas that are close by; two to the right and
    % one to the left (don't forget that between two minima there will
    % always be a maximum!)
    indexList = (closestIdx-2):(closestIdx + 4);
    indexList(indexList < 1 | indexList > length(zeroVals)) = [];
    
    % also, never allow going left of the highest maximum!
    [dummy,maxValIdx] = max(zeroVals);
    indexList(indexList < maxValIdx) = [];
    
    % find lowest
    [dummy, cutIdx] = min(zeroVals(indexList));
    % and determine break value
    cutValue = zeroList(indexList(cutIdx));
    
    % if there is no tag that could possibly be accepted: use cutVal.
    if nnz(ratios > cutValue) < 3
        cutValue = cutVal;
    end
    
    % plot cutoff
    if verbose > 1
        figure('Name',sprintf('cutoff for %s',dataProperties.name)),
        axesH(1)=subplot(2,1,1);
        t1idx = ismember(times,1:3:nTimepoints);
        t2idx = ismember(times,2:3:nTimepoints);
        t3idx = ismember(times,3:3:nTimepoints);
        plot(times(t1idx),ratios(t1idx),'+r',...
            times(t2idx),ratios(t2idx),'+g',...
            times(t3idx),ratios(t3idx),'+b'),
        axesH(2)=subplot(2,1,2);
        cutFirstHistMode(axesH(2),ratios);
        set(axesH(1),'NextPlot','add');
        plot(axesH(1),[1,nTimepoints],[cutVal,cutVal],'g')
        plot(axesH(1),[1,nTimepoints],[cutValue,cutValue],'r')
        set(axesH(2),'NextPlot','add');
        plot(axesH(2),[cutVal,cutVal],[0,100],'g',...
            [cutValue,cutValue],[0,100],'r')
    end
    
end

% clean up
if verbose
    close(h);
end

warning(warningState)

% store cutValue in dataProperties
dataProperties.amplitudeCutoff = cutValue;

% write cutoffs into debugData
debugData.cutoff = [cutVal,cutValue];