function [ output_args ] = combinePartitionedTimeSeries(projList,outDir)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for iProj = 1:numel(projList) 
    
    load([projList{iProj} filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION_20150305' ...
        filesep 'Partition_Outgrowth_Trajectory' filesep 'ForMultRegression' filesep 'forMultReg.mat'])
       
    forMultFinalCell{iProj,1} = forMultReg; 
    names{iProj,1} = paramNamesAll; 
end 

% do a sanity check and make sure the name is correct for 
allNames = vertcat(names{:}); 
nVars = size(allNames,2); 

if sum(arrayfun(@(i) length(unique(allNames(:,i)))==1,1:nVars))/nVars <1
    error('You are collecting variables that are not the same!'); 
end 

forMultRegPooledMat = vertcat(forMultFinalCell{:}); 

% make sure real as I had some problems with the orientation metric on some of these (need to check) 
forMultRegPooledMat = real(forMultRegPooledMat); 

%% if filter by velocity 
vel = forMultRegPooledMat(:,end);
forMultRegPooledMat = forMultRegPooledMat(vel>0.05,:); 

%% 

save([outDir filesep 'forMultRegPooledMat'],'forMultRegPooledMat'); 
pooledMultRegCell = [names{1};num2cell(forMultRegPooledMat)];
names = allNames(1,:); 
% plot ccs 
 nParams = size(forMultRegPooledMat,2); 
 for iParam = 1:nParams
%     crosscor
%setAxis
figure
response = forMultRegPooledMat(:,end); 
descriptor = forMultRegPooledMat(:,iParam); 
      scatter(descriptor,response); 
      ylabel(strrep(names(iParam),'_',' ')); 
      xlabel(strrep(names(end),'_',' ')); 
      [R,P,RLo,Up] = corrcoef(descriptor,response,'rows','Pairwise');
      title({['CC = ' num2str(R(1,2))] ; ['P-Value ' num2str(P(1,2),2)]}); 
      output(iParam).R = R(1,2); 
      output(iParam).P = P(1,2);
      output(iParam).RLo = RLo; 
      output(iParam).RUp = RUp; 
      output(iParam).name = names{iParam}; 
      output(iParam).descriptorVal = descriptor; 
      output(iParam).responseVal = response; 
      saveas(gcf,[outDir filesep names{iParam} '.fig']); 
      saveas(gcf,[outDir filesep names{iParam} '.tif']); 
      close gcf
 end 
dataSetPooled = cell2dataset(pooledMultRegCell); 
mdl = fitlm(dataSetPooled); 
  export(dataSetPooled,'File',[outDir filesep 'VariablesCombined.csv'],'Delimiter',','); 
save([outDir filesep 'model.mat'],'mdl'); 
save([outDir filesep 'output.mat'],'output'); 
end

