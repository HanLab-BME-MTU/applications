
function [] = pcCalcDetectionScaleMD()

% Usese the motion-estimation cross-correlation scores to calculate scale
% for cell detection

% Assaf Zaritsky, June 2015

addpath(genpath('/home2/azaritsky/code/common/mathfun/psfModels'));
addpath(genpath('/home2/azaritsky/code/common/detectionAlgorithms'));

% addpath(genpath('/home2/azaritsky/code/common'));
% addpath(genpath('/home2/azaritsky/code/extern'));
% 
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/algs'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/timeLapseAnalysis'));
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));

warning('off','all');

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

pixelSize = 0.325;
patchSize = ceil(10.0/pixelSize);
timePerFrame = 1;

analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
metaDataFname = [analysisDirname 'MetaData/Experiments20150616_AN.mat'];

load(metaDataFname);

imagesList = [];
for i = 1 : 1 : metaData.tasks.N
    if i > metaData.tasks.N
        return;
    end
    curExp = metaData.tasks.exps(i);
    curTask = metaData.tasks.tasks(i);
    curFname = metaData.experiments.fnames{curExp};
    if curTask <= metaData.experiments.n1{curExp}
        curSource = metaData.experiments.source1{curExp};
    else
        curSource = metaData.experiments.source2{curExp};
    end
    
    %% make sure pcInitiateData was run earlier
    mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
        curFname '_s' sprintf('%02d',curTask) filesep...
        curFname '_s' sprintf('%02d',curTask) '.mat'];
    
    if ~exist(mdFname,'file')
        mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%d',curTask) filesep...
            curFname '_s' sprintf('%d',curTask) '.mat'];
    end
    
    if ~exist(mdFname,'file')
        continue;
    end
    
    MD =  MovieData.load(mdFname);
    for j = 1 : 100 : 301
        mfFname = [MD.outputDirectory_ '/MF/mf/' sprintf('%03d',j) '_mf.mat'];
        
        load(mfFname);
        
        resizedScore = imresize(scores,1.0/patchSize);
        
        imagesList = [imagesList, resizedScore];
    end
end

scale=getGaussianPSFsigmaFromData(imagesList,'Display',true);
disp(['Estimed scales: ' num2str(scale)]);

save([analysisDirname 'scalePointDetection'],'scale');

end
