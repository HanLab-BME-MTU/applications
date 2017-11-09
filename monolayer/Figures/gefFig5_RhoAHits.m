function [] = gefFig5_RhoAHits()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig5/RhoAHits/';
validateGenes = {'RHOA','ARHGEF18','ARHGEF11','ARHGEF28','ARHGEF3'};

outputPrefix = 'RhoAHits_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end