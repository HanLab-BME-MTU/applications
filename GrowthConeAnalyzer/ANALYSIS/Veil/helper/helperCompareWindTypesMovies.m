function [ output_args ] = helperMontageParamScan(projList,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%outAll= ['/project/cellbiology/gdanuser/Harvard/Maria/MethodsPaper/Control/ReRuns/20150707_WindowTypeTest/FINAL_OUTPUT']; 
  
for iProj = 1:numel(projList)
    load([projList{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat'])
   
    GCAVisualsMontagingMovie(MD,inputFinal,outC,['VeilTrackMontage_' ID]); 

end

