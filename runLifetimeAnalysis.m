function [] = runLifetimeAnalysis(restrict,shape,expData)

% runLifetimeAnalysis takes all data folders under a given condition and fills in
% missing tracking and lifetime data
%
% SYNOPSIS [] = runLifetimeAnalysis()
%
% INPUT     restrict    = (optional) time restriction in seconds, e.g. 300
%           shape       = (optional) shape for populations, e.g. [2 2 1], where 1
%                   indicates an exponential distribution, and 2 a Rayleigh
%                   distribution
%           expData  = (optional) if specified, data is taken from this
%                   structure instead of loaded via user interface
%
% OUTPUT
%
% REMARKS   When specifying either the restrict or the shape parameters and
%           not specifying the other parameter leave the non specified
%           parameter as an empty field (e.g. runLiftimeAnalysis(300,[]) or
%           runLifetimeAnalysis([],[2 2 1])
%
% Daniel Nunez, March 5, 2008
% Dinah Loerke, March 15, 2009
% Francois Aguet, Feb 2010

%Default Values
if nargin < 2 || isempty(shape)
    shape = [2 2 1]
end
if nargin < 1 || isempty(restrict)
    restrict = 300;
end
if nargin < 3 || isempty(expData)
    %Load Data, ask user to specify movies to analyse
    [experiment] = loadIndividualMovies();
else
    experiment = expData;
end

currentDate = datestr(now, 'yyyymmdd');

% get condition directory from source path (two levels up)
condDir = strfind(experiment(1).source, filesep);
condDir = experiment(1).source(1:condDir(end-2));

%get condition directory name to use as default for analysis results
dirName = findstr(condDir,filesep);
dirName = condDir(dirName(end-1)+1:dirName(end)-1);

directory = uigetdir(condDir,['Specify folder to store lifetime analysis result.\n If empty default (' condDir ') will be used']);
if (directory==0)
    directory = condDir;
end

%Ask user to name file
fileName = input('Specify name for lifetime analysis result files (date will be included automatically).','s');
if isempty(fileName)
    fileName = [dirName 'LifetimeAnalysisResults'];
end


%Fill in movie length and image size
[experiment] = determinePitDensities(experiment);
%Create lftInfo file if not already there
[experiment] = determineLifetimeInfo(experiment);
%load lftInfo into experiment structures
for k = 1:length(experiment)
    load([experiment(k).source 'LifetimeInfo' filesep 'lftInfo.mat']);
    experiment(k).lftInfo = lftInfo;
    experiment(k).meanPitDensity = mean(experiment(k).pitDensity);
end

%Analyse Data
compactRes = lifetimeCompactFitData(experiment, restrict, shape);

%Save Data
filePath = [directory filesep fileName currentDate];
secureSave(filePath,'compactRes');

if exist([filePath '.txt'],'file')
    directoryFiles = ls(directory);
    fileNumber = length(findstr([fileName currentDate '.txt'], directoryFiles));
    filePath = [directory filesep fileName currentDate '_' num2str(fileNumber)];
    while exist([filePath '.txt'],'file')
        fileNumber = fileNumber + 1;
        filePath = [directory filesep fileName currentDate '_' num2str(fileNumber)];
    end
end

%calculate mean and error of movie densities for condition
meanConditionDensity = nanmean([experiment.meanPitDensity]);
errorConditionDensity = nanstd([experiment.meanPitDensity]);


fid = fopen([filePath '.txt'],'w');
fprintf(fid, '%s',[fileName ',shape,C,deltaC,Tau,deltaTau,Tau50']);
fprintf(fid, '\n');
fprintf(fid, '%s\n',['PP,' num2str(0) ',' num2str(compactRes.contr(1)) ',' num2str(compactRes.contrError(1)) ',' num2str(compactRes.tau(1)) ',' num2str(compactRes.tauError(1)) ',' num2str(compactRes.tau50(1)) ',' num2str(compactRes.numcells) ',cell number']);
fprintf(fid, '%s\n',['P1,' num2str(shape(1)) ',' num2str(compactRes.contr(2)) ',' num2str(compactRes.contrError(2)) ',' num2str(compactRes.tau(2)) ',' num2str(compactRes.tauError(2)) ',' num2str(compactRes.tau50(2)) ',' num2str(compactRes.numtraj) ',traj number']);
fprintf(fid, '%s\n',['P2,' num2str(shape(2)) ',' num2str(compactRes.contr(3)) ',' num2str(compactRes.contrError(3)) ',' num2str(compactRes.tau(3)) ',' num2str(compactRes.tauError(3)) ',' num2str(compactRes.tau50(3)) ',' num2str( 1000*meanConditionDensity) ',density']);
fprintf(fid, '%s\n',['P3,' num2str(shape(3)) ',' num2str(compactRes.contr(4)) ',' num2str(compactRes.contrError(4)) ',' num2str(compactRes.tau(4)) ',' num2str(compactRes.tauError(4)) ',' num2str(compactRes.tau50(4)) ',' num2str(1000*errorConditionDensity) ',delta density']);
fprintf(fid, '%s\n',['Max lifetime for analysis (restrict):' num2str(restrict)]);
fprintf(fid, 'Movies Used:\n');
for iexperiment = 1:length(experiment)
    fprintf(fid,'%s\n',[experiment(iexperiment).source]);
end
fclose(fid);