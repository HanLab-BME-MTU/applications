% job - the experiment (one per node), task - index with in experiment (thread)
function [] = runQLCHOneAtATime(job,task)
addpath(genpath('/work/gdanuser/azaritsky/UTSW/HarvardCode/common'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code'));
addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code'));
warning('off','all');

pixelSize = 0.325;
timePerFrame = 1;
mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/Shanghua/';%'/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/';
metaDataFname = [mainDirname 'PhaseContrast2DExperiments20140410.mat'];

load(metaDataFname);

jobs = metaData.experiments;

% TODO: add defensive programming

fname = jobs.fnames{job};
mainDirname = [mainDirname fname(1:6) '/']; % assuming this directory exists
dirname = [mainDirname fname pad(task,2)];


if task > 20
    fprintf(sprintf('Number of tasks 20 < index %d\n',task));
    return;    
end
fprintf(sprintf('processing %s\n',jobs.fnames{job}));


%% make sure whArrangeData(mainDirname) was run earlier
flags = -1;
pcProcessTimeLapse_stable(pixelSize, timePerFrame, mainDirname, [fname pad(task,2)], flags);

msg = strcat('program with input : job = ', int2str(job), ' task = ', int2str(task), ' ends here.')

end
