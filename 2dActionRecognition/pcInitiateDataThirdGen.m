% dataDirname = '/project/cellbiology/gdanuser/melanomaModel/RawData/2DMorphodynamicsNikon/Data/';
% analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/Data/';
function [] = pcInitiateDataThirdGen(dataDirname,analysisDirname,always)
% Reads all nds2 files from the dataDirname, initiates corresponding MovieData objects, 
% creates corresponding folders in the analysisDirname
% Assaf Zaritsky, January 2017

addpath(genpath('/home2/azaritsky/code/common/toolfun'));
addpath(genpath('/home2/azaritsky/code/extern/bioformats'));
addpath('/home2/azaritsky/code/applications/2dActionRecognition');

if nargin  == 0
    dataDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/raw/ThirdGenData/';    
    analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/ThirdGenData/';    
    always = true;
end

if ~exist(analysisDirname,'dir')
    mkdir(analysisDirname);
end

prepareMovieData(dataDirname,analysisDirname,always);

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
