function [autoCorrDist,autoCorrVel] = makiAutocorrRandPairs(jobType,analysisStruct,verbose)
%MAKIAUTOCORRRANDPAIRS construct random sister pairs and calculates the autocorrelation of sister separation and rate of chance of sister separation
%
%SYNOPSIS [autoCorrDist,autoCorrVel] = makiAutocorrRandPairs(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%               'MCAINSH'
%       analysisStruct: Structure with field movies indicating the movies
%                       to be analyzed. Optional. If not input, GUI to load
%                       movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT autoCorrDist: Autocorrelation of sister separation.
%       autoCorrVel : Autocorrelation of rate of change in sister
%                     separation.
%
%Khuloud Jaqaman, Februaru 2008

%% input

if nargin < 1 || isempty(jobType)
    jobType = 'DANUSER';
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStructs belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%% make random sisters and collect their information

%initialize variables
distanceAll = [];
changeDistAll = [];

%go over all movies
for iMovie = 1 : numMovies
    
    %get tracks in movie
    tracks = dataStruct(iMovie).tracks;
    tracksMat = convStruct2MatNoMS(tracks);

    %get number of tracks and number of pairs to be contructed
    numTracks = length(tracks);
    numPairs = floor(numTracks/2);

    %generate random index of tracks to make pairs
    %     randIndx = randperm(numPairs*2);
    randIndx = (1:numPairs*2);

    %reserve memory for variables
    distance = repmat(struct('observations',[]),numPairs,1);
    changeDist = distance;
    
    %go over all pairs
    for iPair = 1 : numPairs

        %get coordinates of sister 1
        coordX1 = tracksMat(randIndx(2*iPair-1),1:8:end)';
        coordY1 = tracksMat(randIndx(2*iPair-1),2:8:end)';
        coordZ1 = tracksMat(randIndx(2*iPair-1),3:8:end)';

        %get coordinates of sister 2
        coordX2 = tracksMat(randIndx(2*iPair),1:8:end)';
        coordY2 = tracksMat(randIndx(2*iPair),2:8:end)';
        coordZ2 = tracksMat(randIndx(2*iPair),3:8:end)';

        %calculate distance between sisters
        distValues = sqrt( (coordX1-coordX2).^2 + (coordY1-coordY2).^2 + ...
            (coordZ1-coordZ2).^2 );

        %calculate rate of change of distance
        velValues = diff(distValues);

        %store values in structure
        distance(iPair).observations = distValues;
        changeDist(iPair).observations = velValues;

    end

    %add information from current movie to ensemble information
    distanceAll = [distanceAll; distance];
    changeDistAll = [changeDistAll; changeDist];
    
end
   
%% calculate autocorrelation

maxLag = 10;
autoCorrDist = autoCorr(distanceAll,maxLag);
autoCorrVel = autoCorr(changeDistAll,maxLag);

%% plot

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);

    %open figure and hold on
    figure
    hold on

    %plot distance and velocity autocorrelation
    plot((0:maxLag)*timeLapse,autoCorrDist(:,1),'k','marker','.');
    plot((0:maxLag)*timeLapse,autoCorrVel(:,1),'r','marker','.');

    %set axes limit
    axis([0 maxLag*timeLapse min(0,1.1*min([autoCorrDist(:,1);autoCorrVel(:,1)])) 1.1]);

    %write axes labels
    xlabel('Lag (s)');
    ylabel('Autocorrelation');

    %write legend
    text(1*timeLapse,0.9,sprintf([' Black: Sister separation \n Red: ' ...
        'Rate of change of sister separation']));

    %hold off figure
    hold off

end

%% ~~~ the end ~~~
