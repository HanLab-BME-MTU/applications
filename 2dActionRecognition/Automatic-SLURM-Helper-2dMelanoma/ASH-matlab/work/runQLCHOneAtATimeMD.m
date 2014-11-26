% job - the experiment (one per node), task - index with in experiment (thread)
function [] = runQLCHOneAtATimeMD(job,task)
% Assumes the current directory is at 
workdir = '/home2/azaritsky/code/applications/2dActionRecognition/Automatic-SLURM-Helper-2dMelanoma/ASH-matlab/work';
cd(workdir);

addpath(genpath('../../../../../common'));
addpath(genpath('../../../../../extern'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/MET/woundHealing/utils'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/algs/'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/utils/'));
addpath(genpath('../../../../monolayer/utils'));
addpath(genpath('../../../../monolayer/algs'));
addpath(genpath('../../../../monolayer/timeLapseAnalysis'));

% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs'));
addpath(genpath('../../../../2dActionRecognition'));% probably same as just '../../'

warning('off','all');

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

% pixelSize = 0.325;
timePerFrame = 1;
% mainDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/PrimaryMelanoma/';
analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
metaDataFname = [analysisDirname 'Experiments20141119.mat'];

load(metaDataFname);

if ismember(task,metaData.experiments.exclude{job})
    fprintf(sprintf('%s',strcat('program with input : job = ', int2str(job), ' task = ', int2str(task), ' ends here (exclusion list).\n')));
    return;
end

jobs = metaData.experiments;
fname = jobs.fnames{job};

fprintf(sprintf('\nStarting to process %s_s%s\n',fname,pad(task,2)));

%% make sure pcInitiateData was run earlier
MD =  MovieData.load([analysisDirname jobs.source{job} filesep fname filesep fname '_s' pad(task,2) filesep fname '_s' pad(task,2) '.mat']);

flags = -1;
pcProcessTimeLapseMD(MD, timePerFrame, flags);

fprintf(sprintf('%s',strcat('program with input : job = ', int2str(job), ' task = ', int2str(task), ' ends here.\n')));

end
