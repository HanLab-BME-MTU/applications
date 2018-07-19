% ind - index of the experiment that is linked to (node,task) through 
% Assaf Jan. 2017
function [] = runQLCHOneAtATimeThirdGenMD(ind)

close all;

addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/extern'));

% addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/algs'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/timeLapseAnalysis'));
addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));

% warning('off','all');

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

% pixelSize = 0.325;
timePerFrame = 1;
% mainDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/PrimaryMelanoma/';
% analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
metaDataFname = [analysisDirname 'MetaData/ThirdGen20170103.mat'];

load(metaDataFname);

if ind > metaData.tasks.N
    return;
end

curExp = metaData.tasks.exps(ind);
curTask = metaData.tasks.tasks(ind);
curFname = metaData.experiments.fnames{curExp};
% curSource = metaData.experiments.source1{curExp};

fprintf(sprintf('\nStarting to process %s_s%02d\n',curFname,curTask));

%% make sure pcInitiateData was run earlier
mdFname = [analysisDirname 'ThirdGenData/' curFname filesep...
    curFname '_s' sprintf('%02d',curTask) filesep...
    curFname '_s' sprintf('%02d',curTask) '.mat'];

if ~exist(mdFname,'file')
    mdFname = [analysisDirname 'ThirdGenData/' curFname filesep...
    curFname '_s' sprintf('%d',curTask) filesep...
    curFname '_s' sprintf('%d',curTask) '.mat'];
end

if ~exist(mdFname,'file')
    error('No MovieData file named %s\n',mdFname);
end

MD =  MovieData.load(mdFname);

flags = -1;
pcProcessTimeLapseMD(MD, curFname, curTask, timePerFrame, flags);

fprintf(sprintf('%s',strcat('program with input (', int2str(ind), '): job = ', int2str(curExp), ' task = ', int2str(curTask), ' ends here.\n')));

end
