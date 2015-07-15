function [ output_args ] = performPCA(dataSetArray )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%dataMat = dataset2cell(dataSetArray);

% get variable names 
%obsNames = get(dataSetArray,'ObsNames'); 
%varNames = get(dataSetArray,'VarNames'); 

% take out columns
dataMat = cell2mat(dataMat(2:end,2:end)); 
setAxis('on')
h = notBoxPlot(zscore(dataMat));
%getGroupingByDay =  

% 
% set(gca,'XTick',1:numel(varNames)); 
% varNames = cellfun(@(x) strrep(x,'_',' '),varNames,'uniformoutput',0); 
% set(gca,'XTickLabels',varNames); 

%obsNames

% for i = 1:length(h(:).data)
%     
%     
% end 
close gcf
% truncate data : for now take out mean fluorescence intensity AND 
% take out net-velocity 

%dataMat = dataMat(:,1:end-1); 


names = get(dataSetArray,'VarNames'); 
names = names(1:end-1);
obsNames = get(dataSetArray,'ObsNames'); 
obsNames = cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0); 

[coef,scores,latent,ts,exp] = pca(dataMat); 

save('PCADataNoFluor.mat','coef','scores','latent','ts','exp'); 

setAxis('on')
biplot(coef(:,1:3),'Scores',scores(:,1:3),'varlabels',varNames,'obslabels',obsNames); 
saveas(gcf,'PCABiPlotFirst2PCs_NoFluoro.fig'); 

end

