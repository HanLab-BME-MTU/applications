function [ toPlot] = GCAGroupAnalysisCollectMeasurementsForStats( toPlot ,varargin)
%GCACollectData :
% INPUT: toPlot structure with field .info
%                                         .names : name of group
%                                         .projList : projList associated
%                                          with group with individual movie IDs
%                                         .grouping : vector of group
%                                          numbers
% OUTPUT: toPlot structure with all parameters that have been run added
% Note was GCACollectGroupData until 20151027
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');

ip.addParameter('splitMovie',false);
ip.addParameter('splitFrame', 62);  % last frame you want to include
ip.addParameter('OutputDirectory',pwd);
ip.parse(toPlot,varargin{:});
%%
nGroups = numel(toPlot.info.names);

for iGroup = 1:nGroups
    
    % eventually it would be nice to have a matrix detailing what was run
    % already for could just check through for the entire data set.
    
    projListC = toPlot.info.projList{iGroup}(:,1);
    nMovies = size(projListC,1);
    
    if ip.Results.splitMovie == true
        
        % Grouping Var1 : grouping per condition
        % create the grouping variable for pooling full group data
        % [1,(Repeated 2*nCellsProj1 times), 2(Repeated
        % 2*nCellsProj2)....[n,(Repeated 2*ncellsProjN times)]
        grpVar{iGroup} = repmat(iGroup,nMovies*2,1); %
        
        % Grouping Var2  :   grouping per cell
        % [1,1,2,2,3,3,...numCellsTotal,numCellsTotal]
        g = arrayfun(@(x) repmat(x,2,1),1:nMovies,'uniformoutput',0);
        grpVar2{iGroup} = size(vertcat(toPlot.info.projList{1:iGroup-1}),1) + vertcat(g{:});
        
        % Grouping Var3 : grouping per treatment
        g3 = repmat([1,2],1,nMovies)';
        grpVar3{iGroup} = g3 + 2*(iGroup-1);
        
    end % ip.Results.splitMovie
    
    % check for consistency among the parameters that were run.
    for iProj = 1:nMovies
        
        load([projListC{iProj} filesep 'GrowthConeAnalyzer'  filesep 'movieData.mat']);
        parameterDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ]; %
        % for now just search files - redesign so that the parameters in the
        % future are more cleverly named
        
        % might also include ylabel name and ylim for each parameter and
        % read in each Descriptor directory to keep constant.
        
        % search all descriptor parameters.
        localParamFiles = searchFiles('meas_',[],[parameterDir filesep 'Descriptor'],1);
        
        paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
        paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);
        idxOut = cellfun(@(x) strcmpi(x,'maxACFLagSpatial'),paramNamesC);
        paramNamesC(idxOut) = [];
        localParamFiles(idxOut,:) = [];
        
        for iParam = 1:numel(paramNamesC)
            
            % collect the data for each parameter.
            load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
            
            if ip.Results.splitMovie == true
                
                % currently assumes only splitting movie in two pool these
                % values
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+1} = vertcat(measC{1:ip.Results.splitFrame});
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+2} = vertcat(measC{ip.Results.splitFrame+1:end});
                
            else
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{iProj} = vertcat(measC{:});
            end
        end
    end
    %     % collect params and reformat
    paramsAll  = fieldnames(dataSetGroup);
    %
    reformat = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesWholeMovie),1:numel(paramsAll),...
        'uniformoutput',0);
    %
    for iParam = 1:numel(paramsAll)
        toPlot.(paramsAll{iParam}).dataMat{iGroup} = reformat{iParam};
    end
    clear reformat paramsAll dataSetGroup
end  % iGroup

if ip.Results.splitMovie == false
    toPlot.info.groupingPerCell = 1:size(vertcat(toPlot.info.projList{:}),1);
    toPlot.info.groupingPoolWholeMovie = toPlot.info.grouping;
else
    toPlot.info.groupingPoolWholeMovie = vertcat(grpVar{:});
    toPlot.info.groupingPerCell= vertcat(grpVar2{:});
    toPlot.info.groupingPoolBeginEndMovie = vertcat(grpVar3{:});
end

save([ip.Results.OutputDirectory filesep 'toPlotGroupMeas.mat'],'toPlot');
end