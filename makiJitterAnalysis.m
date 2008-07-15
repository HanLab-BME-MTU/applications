function results = makiJitterAnalysis(jobType,analysisStruct,verbose)
%MAKIJITTERANALYSIS calculate the jitter of fixed kinetochores relative
%
%
%SYNOPSIS results = makiJitterAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%               'MCAINSH'
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If no input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 1.
%
%OUTPUT results: 
%           .meanDrift (array with mean drift per movie)
%           .meanJitter (array with mean jitter per movie)
%           .driftMag (cell array {nMovies}(nTimePoints) with drift magnitude [um/frame])
%           .driftDirXY (cell array {nMovies}(nTimePoints) with drift direction [rad])
%           .driftDirZ (cell array {nMovies}(nTimePoints) with drift direction [rad])
%           .jitter (cell array {nMovies}(nTimePoints) with jitter values)
%
%Gaudenz Danuser, July 2008

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'DANUSER';
end

if nargin < 3 || isempty(verbose)
    verbose = 1;
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

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%convert tracks into track matrix
numTracks = zeros(numMovies,1);
for iMovie = 1 : numMovies
    trackMats{iMovie} = convStruct2MatNoMS(dataStruct(iMovie).tracks);
    numTracks(iMovie) = length(dataStruct(iMovie).tracks);
end

%find number of frames in each movie
numFrames = zeros(numMovies,1);
for iMovie = 1 : numMovies
    numFrames(iMovie) = dataStruct(iMovie).dataProperties.movieSize(end);
end

%get time between frames
timeLapse = round(2*dataStruct(1).dataProperties.timeLapse)/2;

% calculate jitter and drift
for iMovie = 1 : numMovies
    % extract X, Y, Z coordinates
    trackX = trackMats{iMovie}(:,[1:8:(numFrames(iMovie)-1)*8+1]);
    trackY = trackMats{iMovie}(:,[2:8:(numFrames(iMovie)-1)*8+2]);
    trackZ = trackMats{iMovie}(:,[3:8:(numFrames(iMovie)-1)*8+3]);
    
    % calculate displacement matrix
    dX = diff(trackX,1,2);
    dY = diff(trackY,1,2);
    dZ = diff(trackZ,1,2);
    
    % calculate drift vector
    driftX = nanmean(dX);
    driftY = nanmean(dY);
    driftZ = nanmean(dZ);
    
    results.driftMag{iMovie} = sqrt(driftX.^2+driftY.^2+driftZ.^2);
    results.driftOriZ{iMovie} = acos(driftZ./results.driftMag{iMovie});
    results.driftOriXY{iMovie} = acos(driftX./(sqrt(driftX.^2+driftY.^2)));
    
    % subtract the drift from the displacement
    jitX = dX - ones(numTracks(iMovie),1)*driftX;
    jitY = dY - ones(numTracks(iMovie),1)*driftY;
    jitZ = dZ - ones(numTracks(iMovie),1)*driftZ;
    
    results.jitter{iMovie} = sqrt(nanmean(jitX.^2)+nanmean(jitY.^2)+nanmean(jitZ.^2));
    
    % get overall statistics over one movie 
    results.meanDrift(iMovie) = mean(results.driftMag{iMovie});
    results.meanJitter(iMovie) = mean(results.jitter{iMovie});
end

if verbose
    
    % prepare graphics
    colVec = ['b','g','r','c','m','k'];
    driftPlot = figure;
    title('Drift magnitude');
    hold on;
    
    jitPlot = figure;
    title('Jitter magnitude');
    hold on;
    
    for iMovie = 1 : numMovies
        
        figure(driftPlot);
        plot(1:(numFrames(iMovie)-1),results.driftMag{iMovie},colVec(mod(iMovie-1,length(colVec))+1));
        
        figure(jitPlot);
        plot(1:(numFrames(iMovie)-1),results.jitter{iMovie},colVec(mod(iMovie-1,length(colVec))+1));
    end
end
%% ~~~ the end ~~~ %%

