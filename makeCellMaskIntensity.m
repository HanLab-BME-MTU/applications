function [movieArray] = makeCellMaskIntensity(experiment,forceRun,closureRadius,fillHoles,smoothSigma,selObject)

% makeCellAreaMask makes mask of cell area based on pixel intensities
%
% Input:
%       experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%                          
%
%       closureRadius  -  Optional. Positive integer. If input, this is the radius of the
%                                   structuring element used to close the mask. Basically it fills in
%                                   any gaps smaller than a disk of this radius. If not input, no
%                                   closure is done.
%
%       fillHoles -         Optional. Logical. If 1, holes in the mask are filled
%                               in. If 0, they aren't.
%
%       smoothSigma  - Optional. Positive scalar real number. This is the
%                                 sigma in pixels of the gaussian spatial filter used to smooth the
%                                 image before masking. Use 1 for most images. Use larger values for
%                                 speckle images, other spotty stuff. If not input, no filtering is
%                                 done.
%
%       selObject -        Optional. Binary. If true, only the largest
%                                 topologically distinct object in the mas is kept
% 
% 
%     forceRun -        If true, masks are created even if the movie has
%                               been previously segmented.
% OUTPUT
%           movieArray
% Uses:
%       batchSegmentMovies
%
% Daniel Nunez, updated May 06, 2009

%save current working directory
od = pwd;

%if run not input for forceRun default to not forcing the run
if nargin < 2 || isempty(forceRun)
    forceRun = 0;
end

%DEFAULTS
if nargin < 3 || isempty(closureRadius)
    closureRadius = 3;
end
if nargin < 4 || isempty(fillHoles)
    fillHoles = 0;
end
if nargin < 5 || isempty(smoothSigma)
    smoothSigma = 2;
end
if nargin < 6 || isempty(selObject)
    selObject = true;
end

%SETUP MOVIEDATA FOR EACH MOVIE
%start variables
movieDataEmpty = struct('analysisDirectory',[],'imageDirectory',[],'nImages',[],'pixelSize_nm',67,'timeInterval_s',[]);
%prep data
for iexp = 1:length(experiment)
    directory = experiment(iexp).source;
    movieDataExists = 0;
        if exist([directory filesep 'movieData.mat'],'file') == 2
            movieDataExists = 1;
        end
    movieData = movieDataEmpty;
    
    if movieDataExists == 0
        movieData.analysisDirectory = directory;
        movieData.imageDirectory = directory;
        movieData.nImages = experiment(iexp).movieLength;
        movieData.timeInterval_s = experiment(iexp).framerate;
        %SAVE MOVIE DATA
        updateMovieData(movieData)
    else
        load([directory filesep 'movieData.mat'])
    end
    movieArray{iexp} = movieData;

end %for each experiment

%SEGMENT MOVIES IN ARRAY
movieArray = batchSegmentMovies(movieArray,closureRadius,fillHoles,smoothSigma,[],selObject,[],forceRun);

%go back to old working directory
cd(od)
end %of function