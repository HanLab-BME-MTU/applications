% List input files
mainFolder = '/Users/sebastien/Documents/Julie/kMTProject';
tiffFiles = dir(fullfile(mainFolder,'*.tif'));
nMovies = numel(tiffFiles);

% Create array of movies
MD(nMovies,1) = MovieData();
for i=1:nMovies
   
   outputDir = fullfile(mainFolder,tiffFiles(i).name(1:end-4));
   if ~isdir(outputDir),mkdir(outputDir); end
   MD(i) = MovieData.load(fullfile(mainFolder,tiffFiles(i).name),false,...
       'outputDirectory',outputDir,'movieDataPath_',outputDir,...
       'movieDataFileName_',[tiffFiles(i).name(1:end-4) '.mat']);
end

% Create list of movies
listFolder = fullfile(mainFolder,'kMTAnalysis');
if ~isdir(listFolder),  mkdir(listFolder); end
ML=MovieList(MD,listFolder);
ML.setPath(listFolder);
ML.setFilename('kMTList.mat');
ML.save


%% Create processes
arrayfun(@reset,MD);
for i=1:nMovies
    MD(i).addPackage(UTrackPackage(MD(i)));
    MD(i).packages_{1}.createDefaultProcess(1);
    MD(i).packages_{1}.createDefaultProcess(2);
    MD(i).addProcess(SisterGroupingProcess(MD(i)));
    MD(i).addProcess(KMTDetectionProcess(MD(i)));
end

%% Set parameters
for i=1:nMovies
    iProc=MD(i).getProcessIndex('SubResolutionProcess',1,false);
    funParams = MD(i).processes_{1}.funParams_;
    funParams.ChannelIndex=2;
    funParams.detectionParam.psfSigma=1.5;
    parseProcessParams(MD(i).processes_{1},funParams);
    
    funParams = MD(i).processes_{2}.funParams_;
    funParams.ChannelIndex=2;
    funParams.costMatrices(1).parameters.diagnostics=[];
    parseProcessParams(MD(i).processes_{2},funParams);
    
    
    parseProcessParams(MD(i).processes_{3},struct('ChannelIndex',2));
    
    parseProcessParams(MD(i).processes_{4},struct('ChannelIndex',2));

end

%%
for i=1:nMovies
    cellfun(@run,MD(i).processes_);
    close all
end


