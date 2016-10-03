function [] = gefFig7_similarity_RAC1_bPIX_CDC42()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_RAC1_bPIX_CDC42/';
validateGenes = {'beta-PIX','RAC1','CDC42'};

outputPrefix = 'RAC1_bPIX_CDC42_';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix);

end