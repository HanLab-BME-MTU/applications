function [ output_args ] = helperCompareWindTypesMovies(projList,outAll)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%outAll= ['/project/cellbiology/gdanuser/Harvard/Maria/MethodsPaper/Control/ReRuns/20150707_WindowTypeTest/FINAL_OUTPUT']; 
  
for iProj = 1:numel(projList)
    load([projList{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat'])
    idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion Sampling')),MD.processes_));
    idxProt = idxProt(2:end); 
   
    input = arrayfun(@(i) (MD.processes_{idxProt(i)}.funParams_.OutputDirectory),1:length(idxProt),'uniformoutput',0); 
    
    for i=1:numel(input) 
       test =  searchFiles('.png',[],input{i},1,'all',1); 
      
        inputFinal{i} = upDirectory(test{1},1); 
        
    end 
   ID =  helperGCACreateID(projList{iProj});
    outC = [outAll filesep ID]; 
    if ~isdir(outC)
        mkdir(outC)
    end 
    
    GCAVisualsMontagingMovie(MD,inputFinal,outC,['VeilTrackMontage_' ID]); 

end

