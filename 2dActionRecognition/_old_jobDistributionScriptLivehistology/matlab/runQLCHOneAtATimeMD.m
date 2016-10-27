% ind - index of the experiment that is linked to (node,task) through 
function [] = runQLCHOneAtATimeMD(ind)

addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/extern'));

addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/algs'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/timeLapseAnalysis'));
addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));

warning('off','all');

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

% pixelSize = 0.325;
timePerFrame = 1;
% mainDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/PrimaryMelanoma/';
analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
metaDataFname = [analysisDirname 'MetaData/Experiments20150521.mat'];

load(metaDataFname);

curExp = metaData.tasks.exps(ind);
curTask = metaData.tasks.tasks(ind);
curFname = metaData.experiments.fnames{curExp};
if curTask <= metaData.experiments.n1{curExp}
curSource = metaData.experiments.source1{ind};
else
    curSource = metaData.experiments.source2{ind};
end

fprintf(sprintf('\nStarting to process %s_s%02d\n',curFname,curTask));

%% make sure pcInitiateData was run earlier
MD =  MovieData.load([analysisDirname 'Data/' curSource filesep curFname filesep...
    curFname '_s' sprintf('%02d',curTask) filesep...
    curFname '_s' sprintf('%02d',curTask) '.mat']);

flags = -1;
pcProcessTimeLapseMD(MD, timePerFrame, flags);

fprintf(sprintf('%s',strcat('program with input (', int2str(ind), '): job = ', int2str(curExp), ' task = ', int2str(curTask), ' ends here.\n')));

end
