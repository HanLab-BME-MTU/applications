function [] = pcInitiateData(dataDirname,analysisDirname,metaDataFname)
% Reads all nds2 files from the dataDirname that exist in the metaDataFname, initiates corresponding MovieData objects, 
% creates corresponding folders in the analysisDirname

% dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/PrimaryMelanoma/';
% analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/';
% metaDataFname = [analysisDirname 'PhaseContrast2DExperimentsMD_debug.mat'];
% metaDataFname = [analysisDirname 'PhaseContrast2DExperiments_pilot.mat'];
% metaDataFname = [analysisDirname 'PhaseContrast2DExperiments20140430.mat'];

% Assaf Zaritsky, May 2014

% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/LCCBcommon/toolfun/movieManagement'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/LCCBcommon/toolfun'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/LCCBextern/bioformats'));

% loci.common.DebugTools.enableLogging('INFO');

load(metaDataFname);

jobs = metaData.experiments;
njobs = jobs.N;

for i = 1 : njobs
    curFname = jobs.fnames{i};
    if exist([curFname '.mat'],'file')
        continue; % the MovieData file was created
    end
    MDs = MovieData([dataDirname curFname '.nd2'],'outputDirectory', [analysisDirname curFname '/' curFname],'reuseReader',true);
    %     MDs = MovieData([dataDirname curFname '.nd2'],'outputDirectory', [analysisDirname curFname '/' curFname]);
end

end