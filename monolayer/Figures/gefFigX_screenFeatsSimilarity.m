function [] = gefFigX_screenFeatsSimilarity()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

controlsDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/';
load([controlsDname 'pcaParamsScreen.mat']);

allControlFeats = [outControls.speed.score(:,1:3),outControls.directional.score(:,1:3),outControls.coordination.score(:,1:3)];

[n,nfeats] = size(allControlFeats);
meanControlFeats = mean(allControlFeats,1);
stdControlFeats = std(allControlFeats,1);

allControlFeatsNorm = (allControlFeats - repmat(meanControlFeats,[n 1])) ./ repmat(stdControlFeats,[n 1]);

featureStrs = {'Speed1','Speed2','Speed3','Direct1','Direct2','Direct3','Coord1','Coord2','Coord3'};
shuffleInds = [1,4,7,2,5,8,3,6,9];
featureStrsShuffle = featureStrs(shuffleInds);
allFeatsNormShuffle = allControlFeatsNorm(:,shuffleInds);
featsSimilarityMat = pdist2(allFeatsNormShuffle',allFeatsNormShuffle','cityblock');
figure; 
imagesc(featsSimilarityMat); 
colormap 'jet';
caxis([prctile(featsSimilarityMat(~(featsSimilarityMat==0)),5),prctile(featsSimilarityMat(~(featsSimilarityMat==0)),95)]);
export_fig([controlsDname 'screenControls_SimilarityFeats.eps']);


end