function [ toPlot] = GCACollectGroupData( toPlot )
%GCACollectData :
% INPUT: toPlot structure with field .info
%                                         .names : name of group
%                                         .projList : projList associated
%                                          with group with individual movie IDs
%                                         .grouping : vector of group
%                                          numbers
% OUTPUT: toPlot structure with all parameters that have been run added 
nGroups = numel(toPlot.info.names);

for iGroup = 1:nGroups
    
    % eventually it would be nice to have a matrix detailing what was run
    % already for could just check through for the entire data set.
    
    projListC = toPlot.info.projList{iGroup}(:,1);
    nMovies = size(projListC,1);
    % check for consistency among the parameters that were run.
    for iProj = 1:nMovies
        load([projListC{iProj} filesep 'ANALYSIS'  filesep 'movieData.mat']);
        parameterDir = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION' ];
        % for now just search files - redesign so that the parameters in the
        % future are more cleverly named
        
        % search all descriptor parameters.
        localParamFiles = searchFiles('param_',[],[parameterDir filesep 'Descriptor'],1);
        
        paramNamesC = cellfun(@(x) strrep(x,'param_',''),localParamFiles(:,1),'uniformoutput',0);
        paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);
        
        for iParam = 1:numel(paramNamesC)
            
            % collect the data for each parameter.
            load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
            dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{iProj} = vertcat(paramC{:});
        end
        
        
        
        
    end
%     % collect params and reformat
    paramsAll  = fieldnames(dataSetGroup);
%     

%     
     reformat = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesWholeMovie),1:numel(paramsAll),...
     'uniformoutput',0);
% 
     for iParam = 1:numel(paramsAll)
         toPlot.(paramsAll{iParam}){iGroup} = reformat{iParam};
     end
    clear reformat paramsAll dataSetGroup
end  % iGroup



end