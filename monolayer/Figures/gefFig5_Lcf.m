function [] = gefFig5_Lcf()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll_Lcf/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig5/Lcf/';
validateGenes = {'Lcf'};

outputPrefix = 'Lcf_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end