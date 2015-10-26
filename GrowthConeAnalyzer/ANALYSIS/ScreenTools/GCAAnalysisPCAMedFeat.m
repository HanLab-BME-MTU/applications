function [ output_args ] = GCAAnalysisPCAMedFeat(toPlot,varargin)
% GCAAnalysisPCAMedFeat : Simple function for exploratory analysis 
%
%
%% REQUIRED 
% 
%  toPlot: (REQUIRED)  group structure for groupAnalysis with the collected data
%                      for each field
%
%% PARAMETERS 
% 
%    'interactive' (PARAM) : logical 
%       if true the user can pick which measurements to include in the
%       analysis 
% 
%    
%% Check Input 
ip = inputParser;
ip.CaseSensitive = false;
% Check input 
ip.addRequired('toPlot'); 
ip.addParameter('interactive',false);
ip.parse(toPlot,varargin{:}); 

%
varNames  = fieldnames(toPlot); 
varNames = varNames(~strcmpi(varNames,'info')) ;

if ip.Results.interactive
paramSelect  = listSelectGUI(varNames,[],'move');
varNames = varNames(paramSelect);
end 

projList = vertcat(toPlot.info.projList{:}); 
obsNames = projList(:,2); 

varNamesString = strrep(varNames, '_' ,'');  
nVars = numel(varNames); 

nGroups = numel(toPlot.info.names); 
grouping = toPlot.info.grouping; 

for iGroup  = 1:numel(toPlot.info.names)
    % transform data so that each row is a neurite and each column is a 
    % set of measurements- this will be a cell per variable. 
    dataMatC = arrayfun(@(x) nanmedian(toPlot.(varNames{x}).dataMat{iGroup},1)',1:nVars,'uniformoutput',0); 
    % put all the measurement values together (median measurement value for each neurite per column) 
    % r indicates the neurite number; 
    dataMatC = horzcat(dataMatC{:}); 
    
    % dataMatAllMeas is a rxc matrix where r is the number of neurites tested 
    % and c is the number of measurements 
    dataMatAllMeas{iGroup} = dataMatC ; % this will keep them in a cell by KD condition 
end
% try once by normalizing each measurement distribution 
% using the pooled population regardless of condition (ie control + KD) 
% this would be mainly to make sure the measurements are on the same scale.
 dataMatAllGroupsMeas = vertcat(dataMatAllMeas{:}); 
 %dataFinal = dataMatAllGroupsMeas; 
 dataFinal = zscoreForNaNs(dataMatAllGroupsMeas); % this should by default normalize along the column 
%% Plot tSNE: test  

% tsneVals = tsne(dataFinal,[],2,0.05); 
% 
% colors = vertcat(toPlot.info.color{:}); 
% grouping =toPlot.info.grouping; 
% setAxis('on')
% hold on 
% arrayfun(@(x) scatter(tsneVals(grouping ==x,1),tsneVals(grouping == x,2),50,colors(x,:),'filled'),1:size(colors,1)); 
% saveas(gcf,'tSNEPlot.fig'); 
% close gcf

%% Norm II :NOTE Do NOT use this normalization this will normalize each population per condition which is not 
%  really what we want. 
% Normalize each measurement distribution to zero mean and standard deviation of 
% 1 for each condition group- these are when each are normalized
% independently - likely not as appropriate 
% dataMatZ = cellfun(@(x) zscoreForNaNs(x),dataMatAllMeas,'uniformoutput',0); 
% dataFinal = vertcat(dataMatZ{:}); 
%% 

[coef,scores,latent,ts,exp] = pca(dataFinal); 
coef = double(coef); 
scores = double(scores); 
save('PCAData.mat','coef','scores','latent','ts','exp'); 
save('PCAZScoresFinal.mat','dataFinal'); 

%% Plot 2D data 
figure
biplot(coef(:,1:2),'varlabels',varNamesString,'obslabels',obsNames); 

saveas(gcf,'PCABiPlotFirst2PCs.fig'); 

scaleFact1 = max(max(abs(scores(:,1:2)))) ; 
coefs = coef(:,1:2); 
maxCoefLen = sqrt(max(sum(coefs.^2,2)));

hold on
arrayfun(@(x) scatter((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
    (scores(grouping==x,2)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups); 


% minx = min(scores(:,1));
% maxx = max(scores(:,1)); 
% 
% miny = min(scores(:,2)); 
% maxy = max(scores(:,2)); 
% 
% line([minx,maxx],[0,0],'color','k'); 
% line([0,0],[miny,maxy],'color','k'); 

% h = legend(toPlot.info.names);
% set(h,'Box','off');
% set(h,'Location','BestOutside');
%view(45,45); 
% xlabel('PCA1'); 
% ylabel('PCA2');
%zlabel('PCA3'); 
saveas(gcf,'ScatterScores2D.fig'); 
close gcf
%% Plot 3D data 
biplot(coef(:,1:3),'varlabels',varNamesString);
saveas(gcf,'PCABiPlotFirst3PCs.fig');
scaleFact1 = max(max(abs(scores(:,1:3)))) ;
coefs = coef(:,1:3);
maxCoefLen = sqrt(max(sum(coefs.^2,2)));

hold on
arrayfun(@(x) scatter3((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
    (scores(grouping==x,2)./scaleFact1).*maxCoefLen,(scores(grouping==x,3)./scaleFact1).*maxCoefLen,...
    50,toPlot.info.color{x},'filled'),1:nGroups);
saveas(gcf,'ScatterScores3D.fig');
%% color code by outgrowth
%figure
%biplot(coef(:,1:2),'scores', scores(:,1:2),'varlabels',varNamesString,'obslabels',obsNames);


figure
biplot(coef(:,1:2),'varlabels',varNamesString,'obslabels',obsNames);

saveas(gcf,'PCABiPlotFirst2PCs.fig');

scaleFact1 = max(max(abs(scores(:,1:2)))) ;
coefs = coef(:,1:2);
maxCoefLen = sqrt(max(sum(coefs.^2,2)));

outgrowth = gcaCollectOutgrowthDeltasPerGroup(toPlot);
outgrowth = vertcat(outgrowth{:});
cmap = brewermap(128,'RdBu');
cmap = flip(cmap,1);
% get the average velocity values for all the windows in the current frame
% and assign a color based on the the mapper.
plotValues = outgrowth;


mapper=linspace(-20,20,128)'; % lower values are red in the brewermap
D=createDistanceMatrix(plotValues,mapper);
[sD,idxCMap]=sort(abs(D),2);

%
hold on
% Plot Events Selected
for iColor = 1:length(cmap)
    if ~isempty(idxCMap(:,1)==iColor)
        scatter((scores(idxCMap(:,1) == iColor,1)./scaleFact1).*maxCoefLen,...
            (scores(idxCMap(:,1) == iColor,2)./scaleFact1).*maxCoefLen,50,cmap(iColor,:),'filled');
    end
end
saveas(gcf,'Scatter2DColoredByOutgrowth.fig');
hold on
obsNames =  cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0);

arrayfun(@(x) text((scores(x,1)./scaleFact1).*maxCoefLen,(scores(x,2)./scaleFact1).*maxCoefLen,obsNames{x}),1:length(scores(:,1)));
%% Plot names 
%figure
%biplot(coef(:,1:2),'varlabels',varNamesString,'obslabels',obsNames);


% arrayfun(@(x) scatter((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
%     (scores(grouping==x,2)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups); 

%  arrayfun(@(x) scatter((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
%      (scores(grouping==x,2)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups); 

%% Discrimination Metrics 
nCond = numel(dataMatAllMeas); 
[DB,Dunn,SC] = arrayfun(@(x) whDiscriminationMeasures(dataMatAllMeas{x}',dataMatAllMeas{1}'),1:nCond); 

names = toPlot.info.names';  
forCell  = num2cell([DB' Dunn']); 
values = [names forCell]; 
% cell2dataset(values); 
discrimValues =cell2table(values); 
save('Table_Dunn_DB','discrimValues'); 
% DB = DB; 

%% 
