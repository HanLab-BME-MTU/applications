function [ reRunList,toPlotHi,toPlotLo ] = GCAAnalysisCompareNeuriteOutgrowthBetweenGroups(toPlot,saveDir,emphasize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fsFigure(0.75);
nGroups = numel(toPlot.info.names);

% check that files exist
for iGroup = 1:nGroups
    projListC = toPlot.info.projList{iGroup}(:,1);
    nofile = cellfun(@(x) exist([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' filesep 'GlobalFunctional' ... 
        filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteLengthOutput.mat'])==0, projListC(:,1));
    idxNoFile = find(nofile);
    if ~isempty(nofile)
        reRunList{iGroup} = projListC(idxNoFile);
        arrayfun(@(x) display(['No Neurite Outgrowth File for ' projListC{idxNoFile(x),1}]),1:length(idxNoFile));
    else
        reRunList{iGroup} = [];
    end
end





if isempty(vertcat(reRunList{:}))
    
    setAxis('on')
    for iGroup = 1:nGroups
        
        %subplot(1,nGroups,iGroup)
        projListC = toPlot.info.projList{iGroup}(:,1);
        
        neuriteLengthStruct = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' filesep 'GlobalFunctional'  ...
            filesep 'neurite_outgrowth_measurements' ...
            filesep 'neuriteLengthOutput.mat']),projListC);
        neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0);
        deltas = cellfun(@(x) (x(end)-x(1)),neuriteLengths);
        deltaPerGroup{iGroup} = deltas;
        neuriteLengthPerGroup{iGroup} = neuriteLengths;
        colorCurrent = toPlot.info.color{iGroup};
         hold on
        h = cellfun(@(x)  plotNeuriteOutgrowthInTime(x,colorCurrent,1,1,0,[],[],0),neuriteLengths,'uniformoutput',0);
       
        forLeg{iGroup} = h(1);
        
        % 
     
        
        if ~isempty(emphasize)
               emphasizeC = emphasize{iGroup}; 
        
            %if iGroup == emphasize(1);
                plotNeuriteOutgrowthInTime(neuriteLengths{emphasizeC},colorCurrent,1,60,0,[],[],1,5);
            %end
        end
%         
%         if ~isempty(saveDir) 
%             saveas(gcf,[saveDir filesep 'neuriteLengthCompare_' toPlot.info.names{iGroup}]); 
%             
%         end 
%         clear neuriteLengths
%         close gcf
    end
    %axis([0,900,-10,35]);
    legend([forLeg{1}{1},forLeg{2}{1},forLeg{3}{1}],toPlot.info.names{1},toPlot.info.names{2},toPlot.info.names{3}, ...
        'Location','BestOutside');
    if ~isempty(saveDir)
        saveas(gcf,[saveDir filesep 'neuriteLengthCompare.fig']);
        saveas(gcf,[saveDir filesep 'neuriteLengthCompare.eps'],'psc2');
        saveas(gcf,[saveDir filesep 'neuriteLengthCompare.tif']);
    end
    
    
    % hAll = vertcat(h{iGroup});
    % allProjs = vertcat(toPlot.info.projList{:});
    % names = allProjs(:,2);
    % hLeg = legend('boxoff');
    % legend(hLeg,names);
    % fig(2) = figure;
    % hold on
    % setAxis
    % arrayfun(@(i) scatter(i*ones(length(deltaPerGroup{i}),1),deltaPerGroup{i},300 ,toPlot.info.colors{i},'filled','MarkerEdge','w'),1:nGroups);
    %  if nGroups >1
    % [h,pValue] = permTest(deltaPerGroup{1},deltaPerGroup{2});
    % title(['P-Value Per Test = ' num2str(pValue,2)]);
    %  end
    % axis([0.5,nGroups+0.5,0,25]);
    % ylabel({'Delta Neurite Outgrowth' ; '(um)'},'FontSize',30);
    % set(gca,'XTick',[1:numel(toPlot.info.names)]);
    %   set(gca,'XTickLabels',[ toPlot.info.names  ]);
    %   if ~isempty(saveDir)
    %   saveas(gcf,[saveDir filesep 'NeuriteOutgrowthCompare.fig']);
    %   saveas(gcf,[saveDir filesep 'NeuriteOutgrowthCompare.png']);
    %   end
    
    deltasCombine = vertcat(deltaPerGroup{:});
    projListCombine = vertcat(toPlot.info.projList{:});
    groupingVar = toPlot.info.grouping;
    
    nGroups = numel(toPlot.info.names);
    
    
    % if perform k-means clustering
    [idx,cCenters]  = kmeans(deltasCombine,2,'replicates',20);
    indexMin = find(cCenters == min(cCenters));
    indexMax = find(cCenters == max(cCenters));
    %  % make it such that 2 is always high and 1 is always the low cluster
    sortedIdx = zeros(length(idx),1);
    sortedIdx(idx==indexMin) = 1;
    sortedIdx(idx==indexMax) = 2;
    
    %
    
    
    
    
    
    
    
    % make project lists for the different cases
    toPlotClust.info.projList= arrayfun(@(iClust) projListCombine(sortedIdx==iClust)...
        ,1:2,'uniformoutput',0);
    
    toPlotClust.info.grouping = arrayfun(@(iClust) groupingVar(sortedIdx==iClust)...
        ,1:2, 'uniformoutput',0);
    toPlotClust.info.deltas = arrayfun(@(iClust) deltasCombine(sortedIdx==iClust)...
        , 1:2,'uniformoutput',0);
    
    colors = toPlot.info.color;
    fsFigure(0.75)
    hold on
    for iClust = 1:2
        %partition further by group
        listC = toPlotClust.info.projList{iClust};
        deltasC = toPlotClust.info.deltas{iClust};
        listByGroup{iClust} = arrayfun(@(iGroup) listC(toPlotClust.info.grouping{iClust}==iGroup),1:nGroups,'uniformoutput',0);
        deltaByGroup{iClust} = arrayfun(@(iGroup) deltasC(toPlotClust.info.grouping{iClust}==iGroup),1:nGroups,'uniformoutput',0);
        arrayfun(@(iGroup) scatter(iClust.*ones(size(listByGroup{iClust}{iGroup},1),1),deltaByGroup{iClust}{iGroup},50,toPlot.info.color{iGroup},'filled')...
            , 1:nGroups);
    end
    axis([0.5 2.5 -5 25]);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabels',{'Cluster Hi Outgrowth' ; 'Cluster Lo Outgrowth'} );
    
    ylabel('Net Outgrowth (um) in 10 min');
    saveas(gcf,'Cluster Cells by Outgrowth.fig'); 
    % get cluster 1
    % write output file that maintains grouping but give different clusters
    % for now keep toPlot Hi and toPlot Low separate. 
    
    toPlotLo.info.projList = listByGroup{1};
    toPlotLo.info.colors = toPlot.info.color; 
   % toPlotHi.info.grouping = cellfun(@(i) ones(length().*iGroup 
    %toPlotHi.info.grouping = vertcat(cellfun(@(i) size(listByGroup,2),    
    
   grpVarLo = arrayfun(@(x) repmat(x,size(listByGroup{1}{x},1),1),1:numel(listByGroup{1}),'uniformoutput',0); 
   toPlotLo.info.grouping = vertcat(grpVarLo{:}); 
    
   % for now just keep in separate groups 
   toPlotHi.info.projList = listByGroup{2}; 
   toPlotHi.info.color = toPlot.info.color; 
 grpVarHi = arrayfun(@(x) repmat(x,size(listByGroup{2}{x},1),1),1:numel(listByGroup{2}),'uniformoutput',0);
   toPlotHi.info.grouping = vertcat(grpVarHi{:}); 
   
   save('toPlotClusters.mat','toPlotLo','toPlotHi'); 
    
    %  % make new toplot file
    %  toPlotLo.info.projList{1} = projListCombine1;
    %  toPlotLo.info.grouping{1} = grouping1;
    %  toPlotLo.info.deltas{1} = deltas1;
    
    
    
    
    
    %projList
    
    
    
    
    
end

