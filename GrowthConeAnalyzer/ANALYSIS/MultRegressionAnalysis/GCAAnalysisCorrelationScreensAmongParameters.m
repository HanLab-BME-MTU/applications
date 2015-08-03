function [ output_args ] = GCAPlotCorrelations(toPlot,outDir)
%

projList = toPlot.info.projList{1};

for iProj = 1:size(projList,1)
    
    %% for now just searchFiles
    parameterDir = [projList{iProj} filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION_20150305' ];
    
    % collect local params for all files in group
    [localParamPaths{iProj}] = searchFiles('param_',[],[parameterDir filesep 'Descriptor'],1);
    
    %  paramNames = cellfun(@(x) strrep(x,'param_',''),localParamFiles(:,1),'uniformoutput',0);
    %  paramNames = cellfun(@(x) strrep(x,'.mat',''),paramNames,'uniformoutput',0);
    %
    %  % for now just create a structure
    %  for iParam = 1:size(localParamFiles,1)
    %      load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
    %      paramCR = reformatDataCell(paramC);
    %      localParams.(paramNames{iParam}) = paramCR;
    %  end
    
end
allParamsRun = vertcat(localParamPaths{:});
paramsRun = unique(allParamsRun(:,1));



idxCorrelate = listSelectGUI(paramsRun,[],'move');


paramsToCorrelate = paramsRun(idxCorrelate);

for iProj =1:size(projList,1)
    
    for iParam = 1:numel(paramsToCorrelate)
        parameterFile = [projList{iProj} filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION_20150305' filesep ...
            'Descriptor'];
        paramFile = searchFiles(paramsToCorrelate{iParam},[],parameterFile,1,'all',1);
        load(paramFile{1});
        forDataSetArray(iProj,iParam) = nanmedian(vertcat(paramC{:})); % collect for all frames
    end
end

x = num2cell(forDataSetArray);
varName = cellfun(@(x) strrep(x,'.mat',[]),paramsToCorrelate,'uniformoutput',0); 
varName = cellfun(@(x) strrep(x,'param_',[]),varName,'uniformoutput',0); 
varName = cellfun(@(x) strrep(x,'filo','F'),varName,'uniformoutput',0); 
varName = cellfun(@(x) strrep(x,'Length','L'),varName,'uniformoutput',0); 
varName = cellfun(@(x) strrep(x,'Intensity','I'),varName,'uniformoutput',0); 

varName = varName';

forDataSetFinal= ([varName ; x]);
paramDataSet = cell2dataset(forDataSetFinal);

[r,p] = corrplot(paramDataSet,'testR','on');
saveas(gcf,[outDir filesep 'correlationScreens']); 

end

