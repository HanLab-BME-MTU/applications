% dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/Data/';
% analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/Data/';
function [] = pcInitiateData(dataDirname,analysisDirname,always)
% Reads all nds2 files from the dataDirname, initiates corresponding MovieData objects, 
% creates corresponding folders in the analysisDirname
% Assaf Zaritsky, May 2015

addpath(genpath('/home2/azaritsky/code/common/toolfun'));
addpath(genpath('/home2/azaritsky/code/extern/bioformats'));

if nargin  == 0
    dataDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/raw/';
    %     analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/Data/';
    analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Data/';    
    always = true;
end

% addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/common/toolfun'));
% addpath(genpath('/project/cellbiology/gdanuser/melanomaModel/Analysis/code/extern/bioformats'));

% loci.common.DebugTools.enableLogging('INFO');

cellLinesDataDir = [dataDirname filesep 'CellLines' filesep];
cellLinesAnalysisDir = [analysisDirname filesep 'CellLines' filesep];
if ~exist(cellLinesAnalysisDir,'dir')
    mkdir(cellLinesAnalysisDir);
end
prepareMovieData(cellLinesDataDir,cellLinesAnalysisDir,always);

melanocytesDataDir = [dataDirname filesep 'Melanocytes' filesep];
melanocytesAnalysisDir = [analysisDirname filesep 'Melanocytes' filesep];
if ~exist(melanocytesAnalysisDir,'dir')
    mkdir(melanocytesAnalysisDir);
end
prepareMovieData(melanocytesDataDir,melanocytesAnalysisDir,always);

TumorsDataDir = [dataDirname filesep 'Tumors' filesep];
TumorsAnalysisDir = [analysisDirname filesep 'Tumors' filesep];
if ~exist(TumorsAnalysisDir,'dir')
    mkdir(TumorsAnalysisDir);
end
prepareMovieData(TumorsDataDir,TumorsAnalysisDir,always);

end

%%
function [] = prepareMovieData(dataDir,analysisDir,always)
if nargin  < 3
    always = false;
end
filenames = dir(dataDir);
nfiles = length(filenames);

for i = 1 : nfiles
    filename = filenames(i).name;
    
    [pathstr, name, ext] = fileparts(filename);
    
    if (strcmp(ext, '.nd2'))
        if exist([analysisDir name],'dir') && ~always
            continue; % the MovieData file was created in corresponding analysis folder
        end
        try
            MDs = MovieData(fullfile([dataDir name '.nd2']),'outputDirectory', fullfile([analysisDir name filesep name]));
        catch ee
            warning(['file ' [dataDir name '.nd2'] 'can not be opend by MovieData']);
            continue;
        end
    end
end
end
