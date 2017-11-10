function [] = gefFig7_similarity_RHOA_RHOC_GEF_SMIFH10um()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/similarity_RHOA_RHOC_GEF_SMIFH10um/';
validateGenes = {'RHOA','ARHGEF18','SMIFH2dose10uM','RHOC','ARHGEF3','ARHGEF28','ARHGEF11'};

outputPrefix = 'RHOA_RHOC_GEF_SMIFH10um';

whDayFollowupSimilarities2016(followupDname,validateGenes,figDname,outputPrefix,[6.5,16.5]);

end