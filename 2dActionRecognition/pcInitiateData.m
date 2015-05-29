% dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/Data/';
% analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/Data/';
function [] = pcInitiateData(dataDirname,analysisDirname)
% Reads all nds2 files from the dataDirname, initiates corresponding MovieData objects, 
% creates corresponding folders in the analysisDirname
% Assaf Zaritsky, May 2015

if nargin  == 0
    dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/Data/';
    analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/Data/';
end

addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/common/toolfun'));
addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/extern/bioformats'));

% loci.common.DebugTools.enableLogging('INFO');

cellLinesDataDir = [dataDirname filesep 'CellLines' filesep];
cellLinesAnalysisDir = [analysisDirname filesep 'CellLines' filesep];
prepareMovieData(cellLinesDataDir,cellLinesAnalysisDir);

melanocytesDataDir = [dataDirname filesep 'Melanocytes' filesep];
melanocytesAnalysisDir = [analysisDirname filesep 'Melanocytes' filesep];
prepareMovieData(melanocytesDataDir,melanocytesAnalysisDir);

TumorsDataDir = [dataDirname filesep 'Tumors' filesep];
TumorsAnalysisDir = [analysisDirname filesep 'Tumors' filesep];
prepareMovieData(TumorsDataDir,TumorsAnalysisDir);

end

%%
function [] = prepareMovieData(dataDir,analysisDir)
filenames = dir(dataDir);
nfiles = length(filenames);

for i = 1 : nfiles
    filename = filenames(i).name;
    
    [pathstr, name, ext] = fileparts(filename);
    
    if (strcmp(ext, '.nd2'))
        if exist([analysisDir name],'dir')
            continue; % the MovieData file was created in corresponding analysis folder
        end
        MDs = MovieData([dataDir name '.nd2'],'outputDirectory', [analysisDir name filesep name]);
    end
end
end
