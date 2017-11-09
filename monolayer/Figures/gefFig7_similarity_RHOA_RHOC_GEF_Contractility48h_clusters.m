function [] = gefFig7_similarity_RHOA_RHOC_GEF_Contractility48h_clusters()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_REVISION/';

% geneClusters = {{'Y15uMp48','Y20uMp48','Blebbistatin10uMp48'},{'RHOA','ARHGEF18'},{'RHOC'},{'ARHGEF3','ARHGEF28','ARHGEF11'}};
% outputPrefix = 'RHOA_RHOC_GEF_Contractility48h_clusters';

geneClusters = {{'Y15uMp48','Y20uMp48','Blebbistatin10uMp48'},{'RHOA','ARHGEF18'},{'RHOC'},{'SMIFH2dose10uM'},{'ARHGEF3','ARHGEF28','ARHGEF11'}};
outputPrefix = 'RHOA_RHOC_GEF_Contractility48h_SMIFH10um_clusters';


clusterFeats = meanClusterFeats(followupDname,geneClusters);

clusterSimilarities(clusterFeats,figDname,outputPrefix,[6.5,16.5]);

end


%%
function clusterFeats = meanClusterFeats(followupDname,geneClusters)

properties = {'Speed','Directionality','Coordination'};
nFeats = 9;
nClusters = length(geneClusters);

clusterFeats = nan(nFeats,nClusters);

for icluster = 1 : nClusters
    curClusterGenes = geneClusters{icluster};
    nCurGenes = length(curClusterGenes);
    
    allFeatsControl = cell(1,nFeats);
    allFeatsDiff = cell(1,nFeats);
    allFeats = [];
    N = 0;
    ns = [];
    
    for i = 1 : nFeats
        allFeatsControl{i} = [];
    end
    
    for iGene = 1 : nCurGenes
        curGene = curClusterGenes{iGene};
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
    
    % Current cluster
    cumsumNs = [0 cumsum(ns)];
    
    featsGene = nan(nFeats,nCurGenes);
    
    for igene = 1 : nCurGenes
        inds = (cumsumNs(igene)+1):cumsumNs(igene+1);
        curGeneFeats = allFeats(:,inds); % curGeneFeats = allFeats(inds,:);
        featsGene(:,igene) = mean(curGeneFeats,2);
    end
    
    featsGeneCluster = mean(featsGene,2);
    
    clusterFeats(:,icluster) = featsGeneCluster;
end

end

%

function [] = clusterSimilarities(clusterFeats,figDname,outputPrefix,caxisGenes)

%% Similarity between genes

fontsize = 10;

[geneSimilarityMat] = getSimilarityMat(clusterFeats); % similarity between genes (cols)

if nargin < 4
    caxisGenes = [prctile(geneSimilarityMat(~(geneSimilarityMat==0)),10),prctile(geneSimilarityMat(~(geneSimilarityMat==0)),90)];
end

figure;
imagesc(geneSimilarityMat);
colormap 'jet';
caxis(caxisGenes);
% export_fig([followupDname filesep outputPrefix 'Similarity_genes.eps']);
export_fig([figDname filesep outputPrefix 'Similarity_genes.eps']);



end


%%
function [simMat] = getSimilarityMat(feats)
feats = feats';
meanFeats = mean(feats,1);
stdFeats = std(feats,1);
zscoreFeats = (feats - repmat(meanFeats,[size(feats,1) 1])) ./ repmat(stdFeats,[size(feats,1) 1]); % zscore(A)
simMat = pdist2(zscoreFeats,zscoreFeats,'cityblock');
end


