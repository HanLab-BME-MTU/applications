function [] = pcInitiateData(dataDirname,analysisDirname,metaDataFname)
% Reads all nds2 files from the dataDirname that exist in the metaDataFname, initiates corresponding MovieData objects, 
% creates corresponding folders in the analysisDirname
% Assaf Zaritsky, May 2014

if nargin  == 0
    dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/';
    analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
    metaDataFname = [analysisDirname 'Experiments20141119.mat'];
end

addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/common/toolfun'));
addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/extern/bioformats'));

% loci.common.DebugTools.enableLogging('INFO');

load(metaDataFname);

jobs = metaData.experiments;
njobs = jobs.N;

for i = 1 : njobs
    curFname = jobs.fnames{i};
    source = metaData.tumors.source{ismember(metaData.tumors.ids,metaData.experiments.tumor1ID{i})}; % assumes tumor1 is present && both tumor's sources are the same
    name = [source filesep curFname];
    if exist([analysisDirname name],'dir')
        continue; % the MovieData file was created in corresponding analysis folder
    end
    MDs = MovieData([dataDirname name '.nd2'],'outputDirectory', [analysisDirname name]); % job1_task1, job1_task2... (for job1/task1, task2... use [analysisDirname name filesep curFname]
    %     MDs = MovieData([dataDirname name '.nd2'],'outputDirectory', [analysisDirname name],'reuseReader',true);
    %     MDs = MovieData([dataDirname curFname '.nd2'],'outputDirectory', [analysisDirname curFname '/' curFname]);
end

end
