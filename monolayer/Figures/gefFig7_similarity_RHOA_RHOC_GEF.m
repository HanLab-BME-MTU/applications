function [] = gefFig7_similarity_RHOA_RHOC_GEF()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_RHOA_RHOC_GEF/';
validateGenes = {'RHOA','ARHGEF18','RHOC','ARHGEF3','ARHGEF28','ARHGEF11'};

outputPrefix = 'RHOA_RHOC_GEF';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix,[6.5,16.5]);

end