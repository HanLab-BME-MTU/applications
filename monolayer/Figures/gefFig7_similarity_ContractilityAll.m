function [] = gefFig7_similarity_ContractilityAll()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_ContractilityAll/';
validateGenes = {'Blebbistatin10uMp48','Y15uMp48','Y20uMp48','Y2763210uM','Y2763210uMp24',...
        'Blebbistatin50uM','Blebbistatin10uM','Blebbistatin25uM',...
        'Blebbistatin10uMp24','Blebbistatin25uMp24',...
        'Y25uMp48'};

outputPrefix = 'ContractilityAll';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix);

end