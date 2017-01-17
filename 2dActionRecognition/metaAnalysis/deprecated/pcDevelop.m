function [] = pcDevelop()

addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code/utils'));
addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code/MET/woundHealing/utils'));
addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code/algs'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/texture'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/PCPhenotyping'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/Hall'));

query.baseDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/Debug/';
query.measures.lbp = true;
query.measures.localMorphDynam = true;
query.rule = 'ruleMetasVsNonMetas';%'ruleTumorInd';%

pcMetaAnalysis(query)
end