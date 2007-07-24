function [trackStats,overallStats] = makiEvalTracks

%allow user to choose directory
basePath = uigetdir('O:','Please select data-dir');

%find all makiData files in chosen directory
fileList = searchFiles('makiData',[],basePath,1);

%allow user to choose movies
selectIdx = listSelectGUI(fileList(:,1),[],'move');

%get number of chosen movies
numMovies = length(selectIdx);

%keep only the chosen movies
movieList = fileList(selectIdx,:);

%load the dataStructs of the chosen movies
for iMovie = 1 : numMovies
    dataFileName = [movieList{iMovie,2} fileSep movieList{iMovie,1}];
    dataStruct(iMovie) = makiLoadDataFile(dataFileName);
end

trackStats(1:numMovies) = struct('numFrames',[],'numSpots',[],...
    'timeWindow',[],'lifeTime',[],'gapSize',[],'aveNumSpots',[],...
    'aveLifeTime',[],'aveGapSizeAll',[],'aveGapSizeInd',[]);

%for each movie ...
for iMovie = 1 : numMovies
    
    %get number of frames in movie
    trackStats(iMovie).numFrames = dataStruct(iMovie).dataProperties.movieSize(end);
    
    %get number of spots in each frame
    initCoord = dataStruct(iMovie).initCoord;
    trackStats(iMovie).numSpots = vertcat(initCoord.nSpots);
    
    %get time window used in gap closing
    trackStats(iMovie).timeWindow = dataStruct(iMovie).dataProperties.tracksParam.gapCloseParam.timeWindow;
    
    %get track lifetime and normalize by number of frames
    trackSEL = getTrackSEL(dataStruct(iMovie).tracks);
    trackStats(iMovie).lifeTime = [trackSEL(:,3) trackSEL(:,3)/trackStats(iMovie).numFrames];
    
    %get gaps in track and normalize by time window
    gapsInTrack = findTrackGaps(dataStruct(iMovie).tracks);
    trackStats(iMovie).gapSize = [gapsInTrack(:,4) gapsInTrack(:,4)/trackStats(iMovie).timeWindow];
    
    %calculate average number of spots in a frame
    trackStats(iMovie).aveNumSpots = mean(trackStats(iMovie).numSpots);
    
    %calculate average lifetime
    trackStats(iMovie).aveLifeTime = mean(trackStats(iMovie).lifeTime);
    
    %calculate average gap size
    trackStats(iMovie).aveGapSizeAll = mean(trackStats(iMovie).gapSize);
    
    %calculate average gap size in each track
    aveGapSizeInd = zeros(size(trackSEL,1),1);
    for iTrack = 1 : size(trackSEL,1)
        aveGapSizeInd(iTrack,1) = mean(gapsInTrack(gapsInTrack(:,1)==iTrack,4));
        trackStats(iMovie).aveGapSizeInd = aveGapSizeInd/trackStats(iMovie).timeWindow;
    end
    
end %(for iMovie = 1 : numMovies)

%plot average number of spots per frame in each movie
figure, hold on
numSpotsAll = vertcat(trackStats.aveNumSpots);
plot(numSpotsAll,'marker','.');
xlabel('movie number');
ylabel('average number of spots per frame');

%plot average normalized lifetime for each movie and normalized lifetime
%histogram for all movies together
figure, hold on
subplot(2,1,1)
lifeTimeAll = vertcat(trackStats.aveLifeTime);
plot(lifeTimeAll(:,2),'marker','.')
xlabel('movie number');
ylabel('average track lifetime');
subplot(2,1,2)
lifeTimeAll = vertcat(trackStats.lifeTime);
hist(lifeTimeAll(:,2),(0:0.05:1))
xlabel('track lifetime');
ylabel('number of tracks');

%plot average normalized gap size for each movie and normalized gap size
%histogram for all movies together
figure, hold on
subplot(2,1,1)
gapSizeAll = vertcat(trackStats.aveGapSizeAll);
plot(gapSizeAll(:,2),'marker','.')
xlabel('movie number');
ylabel('average gap size');
subplot(2,1,2)
gapSizeAll = vertcat(trackStats.gapSize);
hist(gapSizeAll(:,2),(0:0.05:1))
xlabel('gap size');
ylabel('number of tracks');

%plot average size of gaps in each track vs. its lifetime
figure, hold on
aveGapSizeInd = vertcat(trackStats.aveGapSizeInd);
plot(lifeTimeAll(:,2),aveGapSizeInd,'.');

overallStats = [];

