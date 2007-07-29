function analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%MAKISISTERCONNECTIONANALYSIS analyzes in multiple ways the connection between sisters
%
%SYNOPSIS analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: 1: test job (default)
%                2: hercules run
%       analysisStruct: Structure with field movies indicating the movies
%                       to be analyzed. Optional. If not input, GUI to load
%                       movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with addition field
%           .sisterConnection: 
%
%Khuloud Jaqaman, July 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 1;
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%% means and stds
%collect all sisterList's into one array
sisterList = vertcat(dataStruct.sisterList);
numSisterList = size(sisterList,1);

%collect all inter-sister distances into one array
sisterDist = vertcat(sisterList.distances);
sisterDist = sisterDist(~isnan(sisterDist(:,1)),1);

%calculate their mean and standard deviation    
sisterDistMeanStd = [mean(sisterDist) std(sisterDist)];

%calculate change per frame in inter-sister distance
for iList = 1 : numSisterList
    sisterList(iList).vel = diff(sisterList(iList).distances(:,1));
end
sisterVel = vertcat(sisterList.vel);
sisterVel = sisterVel(~isnan(sisterVel));

%convert to nm/s
sisterVel = sisterVel*1000/15;

%calculate the mean and standard deviation of positive velocities
sisterVelPos = sisterVel(sisterVel>0);
sisterVelPosMeanStd = [mean(sisterVelPos) std(sisterVelPos)];

%calculate the mean and standard deviation of negative velocities
sisterVelNeg = abs(sisterVel(sisterVel<0));
sisterVelNegMeanStd = [mean(sisterVelNeg) std(sisterVelNeg)];

%delete some variables
clear sisterDist sisterVel sisterVelPos sisterVelNeg

%% autocorrelation

%define maximum lag
maxLag = 10;

%put the distances and velocities in the correct structure to calculate
%their autocorrelation
for iList = 1 : numSisterList
    sisterDist(iList).observations = sisterList(iList).distances;
    sisterVel(iList).observations = sisterList(iList).vel;
end

%calculate distance autocorrelation
autoCorrDist = autoCorr(sisterDist,maxLag);

%calculate velocity autocorrelation
autoCorrVel = autoCorr(sisterVel,maxLag);

%% ARMA

%% output

%store results in one structure
sisterConnection = struct('separation',sisterDistMeanStd,'rateChangeSepPos',...
    sisterVelPosMeanStd,'rateChangeSepNeg',sisterVelNegMeanStd,...
    'autoCorrDist',autoCorrDist,'autoCorrVel',autoCorrVel);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'sisterConnection');

%store results in analysisStruct
analysisStruct.sisterConnection = sisterConnection;

%if sisterConnection field already existed, add 1 to the version number in
%the file name where analysisStruct will be stored stored
if fieldExists
    [versionNum,fileBody] = makiGetVersion(fileName);
    fileName = [fileBody '_' num2str(versionNum+1) '.mat'];
    analysisStruct.fileName = fileName;
end

%save analysisStruct
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);
    
    %get number of frames in each movie
    numFrames = dataStruct(1).dataProperties.movieSize(end);

    %open figure and write title
    figure('Name',fileName,'NumberTitle','off');
    
    %plot a sample of trajectories
    
    %create subplot 1
    subplot(2,1,1);
    hold on;
    
    %put all distances together in one matrix
    distanceMat = [sisterDist.observations];
    distanceMat = distanceMat(:,1:2:end);
    
    %plot distance over time for all sisters
    plot((1:numFrames)*timeLapse,distanceMat);

    %set axes limit
    axis([timeLapse numFrames*timeLapse 0 max(distanceMat(:))+1]);
    
    %write axes labels
    xlabel('Time (s)');
    ylabel('Sister separation (um)');
    
    %write averaging information
    text(2*timeLapse,max(distanceMat(:))+0.6,sprintf([' Sister separation' ...
        ' (um): %4.2f +- %4.2f \n Rate of change of sister separation ' ...
        '(nm/s): \n +ve: %4.2f +- %4.2f, -ve: %4.2f +- %4.2f'],...
        sisterDistMeanStd(1),sisterDistMeanStd(2),sisterVelPosMeanStd(1),...
        sisterVelPosMeanStd(2),sisterVelNegMeanStd(1),sisterVelNegMeanStd(2)));
    
    %hold off figure
    hold off
    
    %plot autocorrelation functions

    %create subplot 2
    subplot(2,1,2); 
    hold on;

    %plot the distance and velocity autocorrelations
    plot((0:maxLag)*timeLapse,autoCorrDist,'k','marker','.');
    plot((0:maxLag)*timeLapse,autoCorrVel,'r','marker','.');

    %set axes limit
    axis([0 maxLag*timeLapse min(0,1.1*min([autoCorrDist;autoCorrVel])) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Autocorrelation');

    %write legend
    text(1*timeLapse,0.9,sprintf([' Black: Sister separation \n Red: ' ...
        'Rate of change of sister separation']));

    %hold off figure
    hold off

end

%%