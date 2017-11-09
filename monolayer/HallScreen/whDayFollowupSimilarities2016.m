function [] = whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix,caxisGenes)

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

if ~exist(figDname,'dir')
    mkdir(figDname);
end

if nargin == 0
    %    followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
    %    validateGenes = {'Y15uMp48','Y20uMp48','Blebbistatin10uMp48','RHOA','ARHGEF18','ARHGEF3','ARHGEF11','ARHGEF28','Y2763210uM','Blebbistatin25uM','Blebbistatin10uM'}; % 'Blebbistatin50uM','Y20uMp48'
    %    outputPrefix = 'RHOA_GEFs_ContractilityAccuteAnd48_';
   % 2016050516
   %    followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/dayGeneControlFollowup/';
   %    validateGenes = {'RHOA','ARHGEF18','ARHGEF11','ARHGEF28','ARHGEF3','Y2763210uM','Blebbistatin25uM','Blebbistatin10uM'}; % ,'Blebbistatin50uM'
end


properties = {'Speed','Directionality','Coordination'};

% nFeats = 10; % monolayer migration rate + 9 PCs
nFeats = 9; % 9 PCs
nGenes = length(validateGenes);




%% Control statistics
allFeatsControl = cell(1,nFeats);
allFeatsDiff = cell(1,nFeats);
allFeats = [];
N = 0;
ns = [];

for i = 1 : nFeats 
    allFeatsControl{i} = [];
end

for iGene = 1 : nGenes
    curGene = validateGenes{iGene};
    curFeats = [];
    for iProp = 1 : length(properties)
        curProp = properties{iProp};
        load([followupDname filesep curProp filesep curGene filesep curProp '_stats.mat']);                        
        
        if iProp == 1            
            N = N + length(controlRepPC1);
            ns = [ns length(controlRepPC1)];
        end
                
        diffPC1 = controlRepPC1 - geneRepPC1;
        diffPC2 = controlRepPC2 - geneRepPC2;
        diffPC3 = controlRepPC3 - geneRepPC3;
        
        j = (iProp-1)*3+1;
        allFeatsControl{j} = [allFeatsControl{j}, controlRepPC1];
        allFeatsDiff{j} = [allFeatsControl{j}, diffPC1];
        j = (iProp-1)*3+2;
        allFeatsControl{j} = [allFeatsControl{j}, controlRepPC1];
        allFeatsDiff{j} = [allFeatsControl{j}, diffPC2];
        j = (iProp-1)*3+3;
        allFeatsControl{j} = [allFeatsControl{j}, controlRepPC1];
        allFeatsDiff{j} = [allFeatsControl{j}, diffPC3];
        
        curFeats = [curFeats;diffPC1];
        curFeats = [curFeats;diffPC2];
        curFeats = [curFeats;diffPC3];
    end
    allFeats = [allFeats,curFeats];
end

cumsumNs = [0 cumsum(ns)];

meanControl = nan(1,nFeats);

stdControl = nan(1,nFeats);
for i = 1 : nFeats 
    meanControl(i) = mean(allFeatsControl{i});
    stdControl(i) = std(allFeatsControl{i});
end

%% PCA
[out] = whGetPCA(allFeats);
allFeatsNorm = out.curFeatsNorm;


%% Mean measure (normalized using diffs)
featsGene = nan(nFeats,nGenes);

for igene = 1 : nGenes
    inds = (cumsumNs(igene)+1):cumsumNs(igene+1);
    curGeneFeats = allFeatsNorm(inds,:);
    featsGene(:,igene) = mean(curGeneFeats,1)';
end
    
%% Similarity between genes

fontsize = 10;

[geneZscoreFeats,geneSimilarityMat] = getSimilarityMat(featsGene); % similarity between genes (cols)

if nargin < 5
    caxisGenes = [prctile(geneSimilarityMat(~(geneSimilarityMat==0)),10),prctile(geneSimilarityMat(~(geneSimilarityMat==0)),90)];
end

figure; 
imagesc(geneSimilarityMat); 
colormap 'jet';
% caxis([prctile(geneSimilarityMat(~(geneSimilarityMat==0)),5),prctile(geneSimilarityMat(~(geneSimilarityMat==0)),95)]);
caxis(caxisGenes);
export_fig([followupDname filesep outputPrefix 'Similarity_genes.eps']);
export_fig([figDname filesep outputPrefix 'Similarity_genes.eps']);

genClusterTreeFeats = linkage(featsGene');
D = pdist(featsGene','cityblock');
geneLeafOrderZScore = optimalleaforder(genClusterTreeFeats,D);

% FPosition = [0 0 1200 700];
% APosition = [0.1 0.2 0.7 0.75]; 

% h = figure;
% dendrogram(genClusterTreeFeats,'Labels',validateGenes,'Reorder',geneLeafOrderZScore); % ,'Reorder',geneLeafOrderZScore,'ColorThreshold','default'
% % hold on;
% % set(h,'Position',FPosition,'PaperPositionMode','auto');
% % axisHandle= findobj(h,'type','axes');
% % set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
% hold off;
% export_fig([followupDname filesep outputPrefix 'Clusters_genes.eps']);

% PCA
[outGenes] = whGetPCA(featsGene);

scoreGenes = outGenes.score;
accVarianceGenes = outGenes.accVariance;

cmap = colormap(hsv(nGenes));

figure;
hold on;
title(sprintf('PCA 1/2, var = %.2f, %.2f ',accVarianceGenes(1),accVarianceGenes(2)));
for igene = 1 : nGenes
plot(scoreGenes(igene,1),scoreGenes(igene,2),'o','MarkerEdgeColor',cmap(igene,:),'LineWidth',2,'MarkerSize',12);
end
legend(validateGenes,'Location','eastoutside');
hold off;
export_fig([followupDname filesep outputPrefix 'PCA_genes_legend.eps']);
export_fig([figDname filesep outputPrefix 'PCA_genes_legend.eps']);
legend off; 
export_fig([followupDname filesep outputPrefix 'PCA_genes.eps']);
export_fig([figDname filesep outputPrefix 'PCA_genes.eps']);


%% Similarity between measures
featureStrs = {'Speed1','Speed2','Speed3','Direct1','Direct2','Direct3','Coord1','Coord2','Coord3'};
shuffleInds = [1,4,7,2,5,8,3,6,9];
featureStrsShuffle = featureStrs(shuffleInds);
allFeatsNormShuffle = allFeatsNorm(:,shuffleInds);
featsSimilarityMat = pdist2(allFeatsNormShuffle',allFeatsNormShuffle','cityblock');
figure; 
imagesc(featsSimilarityMat); 
colormap 'jet';
caxis([prctile(featsSimilarityMat(~(featsSimilarityMat==0)),5),prctile(featsSimilarityMat(~(featsSimilarityMat==0)),95)]);
export_fig([followupDname filesep outputPrefix 'Similarity_feats.eps']);
export_fig([figDname filesep outputPrefix 'Similarity_feats.eps']);

% featsClusterTreeFeats = linkage(allFeatsNormShuffle');
% D = pdist(allFeatsNormShuffle','cityblock');
% featsLeafOrderZScore = optimalleaforder(featsClusterTreeFeats,D);
% figure();
% dendrogram(featsClusterTreeFeats,'Labels',featureStrsShuffle,'Reorder',featsLeafOrderZScore); % ,'Reorder',geneLeafOrderZScore,'ColorThreshold','default'
% export_fig([followupDname  filesep outputPrefix 'Clusters_feats.eps']);

% PCA
[outFeats] = whGetPCA(allFeatsNormShuffle);

scoreFeats = outFeats.score;
accVarianceFeats = outFeats.accVariance;

cmap = colormap(hsv(nFeats));

figure;
hold on;
title(sprintf('PCA 1/2, var = %.2f, %.2f ',accVarianceFeats(1),accVarianceFeats(2)));
for ifeat = 1 : nFeats
    plot(scoreFeats(ifeat,1),scoreFeats(ifeat,2),'o','MarkerEdgeColor',cmap(ifeat,:),'LineWidth',2,'MarkerSize',12);
end
legend(featureStrsShuffle,'Location','eastoutside');
hold off;
export_fig([followupDname filesep outputPrefix 'PCA_feats_legend.eps']);
export_fig([figDname filesep outputPrefix 'PCA_feats_legend.eps']);
legend off; 
export_fig([followupDname filesep outputPrefix 'PCA_feats.eps']);
export_fig([figDname filesep outputPrefix 'PCA_feats.eps']);


%% Sort by (1) genes, (2) features
% geneZscoreFeats_sortGenes = geneZscoreFeats(geneLeafOrderZScore,:);
% geneZscoreFeats_sortGenes_sortFeats = geneZscoreFeats_sortGenes(:,shuffleInds);
% 
% figure; 
% imagesc(geneZscoreFeats_sortGenes_sortFeats'); 
% colormap 'jet';
% caxis([-0.7,0.7]);
% export_fig([followupDname  filesep outputPrefix 'Sorted.eps']);


end

function feat = normDiffMedian(cntl,kd,meanControl,stdControl)
normCntrl = (cntl - meanControl) ./ stdControl;
normKD = (kd - meanControl) ./ stdControl;

diffs = normCntrl - normKD;

feat = median(diffs);
end


function [zscoreFeats, simMat] = getSimilarityMat(feats)
feats = feats';
meanFeats = mean(feats,1);
stdFeats = std(feats,1);
zscoreFeats = (feats - repmat(meanFeats,[size(feats,1) 1])) ./ repmat(stdFeats,[size(feats,1) 1]); % zscore(A)
simMat = pdist2(zscoreFeats,zscoreFeats,'cityblock');
end

