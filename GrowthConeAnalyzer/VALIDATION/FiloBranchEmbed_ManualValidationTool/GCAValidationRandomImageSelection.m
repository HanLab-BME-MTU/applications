function [ selectedFiles ] = GCAValidationRandomImageSelection( projList,nSamples )
%GCAValidationReconstructRandomSelection
% this function randomly selects a project and frame number from a list of
% projects for validation. 
if isstruct(projList) 
    toPlot = projList; 
    projList = vertcat(toPlot.info.projList{:});
    projList = projList(:,1);
end 
% make the list of IDs 
projectsAll = repmat(projList,120,1); 
frames = arrayfun(@(i) repmat(i,length(projList),1),1:120,'uniformoutput',0); 
frames = vertcat(frames{:}); 

projectSamp = projectsAll(randsample(1:length(frames),nSamples),:); 
framesSamp = frames(randsample(1:length(frames),nSamples)); 

selectedFiles = [projectSamp  num2cell(framesSamp)];  
save('selectedFiles.mat','selectedFiles'); 

end

