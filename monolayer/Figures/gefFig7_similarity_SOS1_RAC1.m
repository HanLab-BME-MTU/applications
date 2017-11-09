function [] = gefFig7_similarity_SOS1_RAC1()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_SOS1_RAC1/';
validateGenes = {'RAC1','SOS1','GSKp96','PDp96','ERKip96'};

outputPrefix = 'SOS1_RAC1_';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix);

end