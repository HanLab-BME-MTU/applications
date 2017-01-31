function [] = gefFig7_Rhosin()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/dayGeneControlFollowup/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig7/Rhosin/';
validateGenes = {'Rhosin12dot5uM','Rhosin25uM','Rhosin50uM','Rhosin100uM','Rhosin200uM'};

outputPrefix = 'Rhosin_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end