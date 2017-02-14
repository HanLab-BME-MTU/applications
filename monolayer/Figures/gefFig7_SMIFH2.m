function [] = gefFig7_SMIFH2()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/SMIFH2/';
validateGenes = {'SMIFH2dose5uM','SMIFH2dose10uM','SMIFH2dose15uM','SMIFH2dose25uM'};

outputPrefix = 'SMIFH2_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end