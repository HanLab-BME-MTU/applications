function [] = gefRevision_RhoC_hairpin3()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

followupDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612_RhoC3/dayGeneControlFollowup/';

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Revision/RHOC_hairpin3/';
validateGenes = {'RHOC'};
outputPrefix = 'RHOC_3_';

whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix);

close all;
end