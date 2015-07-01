function [ output_args ] = GCAAnalysisPCAMedFeat( toPlot )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
varNames  = fieldnames(toPlot); 
projList = vertcat(toPlot.info.projList{:}); 
obsNames = projList(:,2); 
varNames = varNames(~strcmpi(varNames,'info')) ;
varNamesString = strrep(varNames, '_' ,'');  
nVars = numel(varNames); 
grouping = toPlot.info.grouping; 
for iGroup  = 1:numel(toPlot.info.names)
    dataMatC = arrayfun(@(x) toPlot.(varNames{x}){iGroup}',1:nVars,'uniformoutput',0); 
    dataMatC = horzcat(dataMatC{:}); 
    dataMatAll{iGroup} = dataMatC ; 
end
dataMatZ = cellfun(@(x) zscoreMine(x),dataMatAll,'uniformoutput',0); 
dataFinal = vertcat(dataMatZ{:}); 
[coef,scores,latent,ts,exp] = pca(dataFinal); 

save('PCAData.mat','coef','scores','latent','ts','exp'); 


biplot(coef(:,1:2),'Scores',scores(:,1:2),'varlabels',varNamesString,'obslabels',obsNames); 

saveas(gcf,'PCABiPlotFirst2PCs.fig'); 
figure;
hold on
arrayfun(@(x) scatter(scores(grouping==x,1),scores(grouping==x,2),toPlot.info.colors{x},'filled'),1:3); 
hold on 

minx = min(scores(:,1));
maxx = max(scores(:,1)); 

miny = min(scores(:,2)); 
maxy = max(scores(:,2)); 

line([minx,maxx],[0,0],'color','k'); 
line([0,0],[miny,maxy],'color','k'); 

% h = legend(toPlot.info.names);
% set(h,'Box','off');
% set(h,'Location','BestOutside');
%view(45,45); 
xlabel('PCA1'); 
ylabel('PCA2');
%zlabel('PCA3'); 
saveas(gcf,'ScatterScores.fig'); 



