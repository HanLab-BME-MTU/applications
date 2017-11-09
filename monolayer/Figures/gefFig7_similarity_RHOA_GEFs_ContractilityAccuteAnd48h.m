function [] = gefFig7_similarity_RHOA_GEFs_ContractilityAccuteAnd48h()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/';
validateGenes = {'Y15uMp48','Y20uMp48','Blebbistatin10uMp48','RHOA','ARHGEF18','ARHGEF3','ARHGEF11','ARHGEF28','Y2763210uM','Blebbistatin25uM','Blebbistatin10uM'}; % 'Blebbistatin50uM','Y20uMp48'
outputPrefix = 'RHOA_GEFs_ContractilityAccuteAnd48_';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix);

end