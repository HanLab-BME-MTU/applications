function [] = gefRevision_RhoARhoC()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612/dayGeneControlFollowup/';
% figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Revision/RHOA_RHOC/';
% validateGenes = {'RHOA','RHOC'};
% outputPrefix = 'RHOARHOC_';

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Revision/RHOC_new/';
validateGenes = {'RHOC'};
outputPrefix = 'RHOCnew_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end