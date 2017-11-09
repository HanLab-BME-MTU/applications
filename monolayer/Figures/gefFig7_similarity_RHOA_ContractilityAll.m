function [] = gefFig7_similarity_RHOA_ContractilityAll()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/';
validateGenes = {'RHOA','Blebbistatin10uMp48',...
    'Y20uMp48','Y2763210uM','Y15uMp48','Y2763210uMp24',...
    'Blebbistatin25uM',...
    'Blebbistatin10uMp24','Blebbistatin25uMp24',...
    'Blebbistatin10uM','Blebbistatin50uM','Y25uMp48'};

outputPrefix = 'RHOA_ContractilityAll_';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix);

end