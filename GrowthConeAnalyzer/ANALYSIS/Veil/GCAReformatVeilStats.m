function [ toPlot] = GCAReformatVeilStats(toPlot,clearOldFields,saveDir)
% Here just want the groupData with the names and the relevant project
% lists for now.
if nargin<3
    saveDir =[];
end

if clearOldFields == 1
    params = fieldnames(toPlot);
    params = params(~strcmpi(params,'info')); % keep the info
    toPlot = rmfield(toPlot,params);
end


% Collect

for iGroup = 1:numel(toPlot.info.names)
    
    projListC = toPlot.info.projList{iGroup};
    % sort projList by outgrowth in 10 min
    %     delta = cellfun(@(x) x(end)-x(1),groupData.netOutgrowthTimeSeries{iGroup});
    %     [sorted,idx] = sort(delta);
    
    %     projListC = projListC(idx,:);
    %     groupData.projList{iGroup} = projListC; % save the new project list
    %
    
    
    %% if collect veil params
    % make a boxplot of the distribution for individual groups and then
    % pooled.
    % make .csv files for individual window values
    
    
    
    params{1} = 'mednVeloc';
    params{2} = 'persTime';
    
    analType{1} = 'protrusionAnalysis';
    analType{2} = 'retractionAnalysis';
    %
    %
    [nplots,~] = size(projListC);
    grpVar{iGroup} = repmat(iGroup,nplots,1); % grouping variable is a repeat matrix the number of cells long
    
    for iAnal = 1:2
        for iParam= 1:2
            [nProjs, ~]= size(projListC);
            % initialize mat
            dataMat = nan(7000,nProjs); % over initialize
            
            
            
            
            for iProj = 1:size(projListC,1)
                folder = [projListC{iProj,1} filesep 'GrowthConeAnalyzer']; 
%                 if isempty(regexp(projListC{iProj},'GrowthConeAnalyzer','ONCE')) == 1
%                     folder = [projListC{iProj} filesep 'ANALYSIS'];
%                 else
%                     folder = projListC{iProj};
%                 end
%                 
                toLoad = [folder filesep 'EdgeVelocityQuantification' filesep 'EdgeMotion.mat'  ];
                
                if exist(toLoad,'file')~=0
                    load(toLoad);
                
                    % want to collect prot persTime, retract persTime,prot mednvel in
                    % my old format for boxplot so i have the values groupable and I
                    % can do myown stats on them.
                    
                    valuesC =  analysisResults.(analType{iAnal}).total.(params{iParam}); % collect all the prot and or retract params
                    valuesC = valuesC(~isnan(valuesC));
                    if length(valuesC) >7000
                        display(['Increase the number of rows for the dataMat initialization to greater than ' num2str(length(valuesC))]);
                        
                    end
                    dataMat(1:length(valuesC),iProj) = valuesC;
                else 
                    display(['No File Found for ' projListC{iProj,1}]); 
                end
            end
            toPlot.([analType{iAnal} '_' params{iParam}]){iGroup} = dataMat;
            
            
        end % iAnal
        
    end % iParam
    
    
    % collect the larger boxplot
    
    
    
    % plot the boxplots
    %         subplot(nrows,3,iProj)
    %
    %
    %
    %         subplot(nrows,3,iProj);
    %         name = strrep(projListC{iProj,2},'_',' ');
    %
    
    
    %     saveas(gcf,[saveDir filesep 'protrusionMaps_Group' groupData.names{iGroup} '.eps'],'psc2');
    %     saveas(gcf,[saveDir filesep 'protrusionMaps_Group' groupData.names{iGroup} '.fig']);
    
    % toPlot.info.grouping = vertcat(grpVar{:});
    % toPlot.info.colors = groupData.color;
    % toPlot.info.projList = groupData.projList;
    %toPlot.info.names = groupData.names;
    if ~isempty(saveDir)
        
        save([saveDir filesep 'toPlotVeil.mat'],'toPlot');
    end
    close gcf
    
    
end





