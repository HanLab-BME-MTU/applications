function [] = runHallOneAtATime(ind)
% Temporarily!
% addpath(genpath('../../../../../common'));
% % addpath(genpath('../../../../../extern'));
% addpath(genpath('../../../../monolayer'));

% Temporarily! -- to remove!!
addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));

% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/MET/woundHealing/utils'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/algs/'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/utils/'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/Hall'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs'));
addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/MET/woundHealing/motion/clusterMF'));
addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/algs/'));

warning('off','all');

timePerFrame = 5;
initParams.nRois = 1;
% initParams.always = true; % TEMPORARY!
% mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/Data20140224/';
% mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/POC_Nov21/';
% mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/Data20140331/';
% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data20140603/';
% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data20140627/';
% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data20140701/';
% mainDirname ='/project/cellbiology/gdanuser/collab/hall/Data20141013/';
% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data20141201/';

mainDirname = '/project/cellbiology/gdanuser/collab/hall/allData/';

% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data201504/';
% mainDirname = '/project/cellbiology/gdanuser/collab/hall/Data20150424/';


% load([mainDirname 'GTPasesScreenMetaData20140121_new.mat']); % exps
% load([mainDirname 'GTPasesScreenMetaData20140401_new.mat']); % exps
% load([mainDirname 'GTPasesScreenMetaData20140602_new.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20140627_Angeles.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20140703_Angeles_new.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20140703_all.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20140829_all_exclude112.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20141013_new_kd50.mat']); % metaData
% load([mainDirname 'GEFsScreenMetaData20141013_kd50.mat']);

% load([mainDirname 'GEFsScreenMetaData20150415_new_kd50.mat']);
% load([mainDirname 'GEFsScreenMetaData20150424_Angeles_new_kd0.mat']);
load([mainDirname 'GEFsScreenMetaData20150424_kd0.mat']);

if ind > length(metaData.fnames)
    fprintf(sprintf('Number of experiments %d < index %d\n',length(metaData.fnames), ind));
    return;
end
fprintf(sprintf('processing %s\n',metaData.fnames{ind}));
whProcessTimeLapse(metaData.pixelSize{ind}, timePerFrame, mainDirname, metaData.fnames{ind}, initParams);
close all;
end