function [ output_args ] = performCorticalStats(projList)
%I think I will try to keep this generic so you 

if isstruct(projList)
    projList = arrayfun(@(x) projList(x).anDir ,1:length(projList), 'uniformoutput',0);
end
for iProj = 1:numel(projList)
    if ~exist([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
        display(['No Growth Tracks in subRoi for' projList{iProj}])
        
    else
        s1 = load([projList{iProj} filesep 'meta' filesep 'projData.mat']);
        projData = s1.projData;
        
        
        s2 = load([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
        corticalData = s2.corticalData;
        
        addEndOnInfo(corticalData,3,60,0.7,0.3,projData);
        %addEndOnInfoWithClassifications(corticalData,3,60,0.7,0.3,projData); 
        
        
    end
end 

