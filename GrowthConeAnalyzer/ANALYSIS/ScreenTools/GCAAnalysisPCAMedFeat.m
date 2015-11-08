function [ output_args ] = GCAAnalysisPCAMedFeatGeneric(toPlot,varargin)
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

projList = vertcat(toPlot.info.projList{:}); 
obsNames = projList(:,2); 

varNamesString = strrep(varNames, '_' ,'');  
nVars = numel(varNames); 

nGroups = numel(toPlot.info.names); 
grouping = toPlot.info.grouping; 

%% User can define which measurements(variables) to include if interactive turned on 
if ip.Results.interactive
    paramSelect  = listSelectGUI(varNames,[],'move');
    varNames = varNames(paramSelect);
end
%% Collect and normalize the measurement data : Default is to take a median value for each cell
for iGroup  = 1:numel(toPlot.info.names)
    % transform data so that each row is a neurite and each column is a 
    % set of measurements- this will be a per cell variable. 
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
 
%% Collect the Outgrowth Data, Create Mapper, and Save Colorbar
fsFigure(0.75,'visible','off');
test = -2:2; 
imagesc(test);
outgrowth = gcaCollectOutgrowthDeltasPerGroup(toPlot);
outgrowth = outgrowth./10 ; % um per min 
outgrowth = vertcat(outgrowth{:});
cmap = brewermap(128,'RdBu');
cmap = flip(cmap,1);
% get the average velocity values for all the windows in the current frame
% and assign a color based on the the mapper.
plotValues = outgrowth;

mapper=linspace(-2,2,128)'; % lower values are red in the brewermap
D=createDistanceMatrix(plotValues,mapper);
[sD,idxCMap]=sort(abs(D),2);
colormap(cmap); 
colorbar
saveas(gcf,'colorbar.fig'); 
saveas(gcf,'colorbar.eps','psc2'); 

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

%% Perform the PCA
[coef,scores,latent,ts,exp] = pca(dataFinal); 
coef = double(coef); 
scores = double(scores); 
save('PCAData.mat','coef','scores','latent','ts','exp'); 
save('PCAZScoresFinal.mat','dataFinal'); 

%% Plot variance explained. 
setAxis

cumExp = cumsum(latent)./sum(latent);
PC = 1:length(cumExp); 
scatter(PC,cumExp,20,'filled','k');
ylabel('Variance Explained');
xlabel('PC Number');
forLine = PC(find(cumExp>.95,1,'first')); 
line([forLine,forLine],[0:1],'color','k'); 

saveas(gcf,'PercentVarianceExplained.fig');
saveas(gcf,'PercentVarianceExplained.eps','psc2');
saveas(gcf,'PercentVarianceExplained.png'); 

%% 2D Plots 
for iPC = 1:2; 
    %% Biplot  
       
    fsFigure(0.75)
    biplot(coef(:,iPC:iPC+1),'varlabels',varNamesString,'obslabels',obsNames);
    xlabel({['PC' num2str(iPC)] ; ['Percent Variance Explained ' num2str(exp(iPC),3) '%']});
    ylabel({['PC' num2str(iPC+1)] ; ['Percent Variance Explained ' num2str(exp(iPC+1),3) '%']});
    
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.fig']);
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.eps'],'psc2');
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.png']);
    
    %% Scores (Per Neurite) Overlay ColorCoded By Experimental Condition
    scaleFact1 = max(max(abs(scores(:,iPC:iPC+1)))) ;
    coefs = coef(:,iPC:iPC+1);
    maxCoefLen = sqrt(max(sum(coefs.^2,2)));
    
    hold on
    arrayfun(@(x) scatter((scores(grouping==x,iPC)./scaleFact1).*maxCoefLen,...
        (scores(grouping==x,iPC+1)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups);
    
    grid('off')
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.fig' '.fig']);
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.eps'],'psc2');
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.png']);
    
    %% 
%     % Plot Events Selected
%     for iColor = 1:length(cmap)
%         if ~isempty(idxCMap(:,1)==iColor)
%             scatter((scores(idxCMap(:,1) == iColor,1)./scaleFact1).*maxCoefLen,...
%                 (scores(idxCMap(:,1) == iColor,2)./scaleFact1).*maxCoefLen,50,'MarkerEdgeColor',cmap(iColor,:));
%         end
%     end
%     
    
      
    close gcf
    %% Scores (Per Neurite) Overlay ColorCoded By Outgrowth 
   
    fsFigure(0.75)
    biplot(coef(:,iPC:iPC+1),'varlabels',varNamesString,'obslabels',obsNames);
    
    hold on
    % Plot Events Selected
    for iColor = 1:length(cmap)
        if ~isempty(idxCMap(:,1)==iColor)
            scatter((scores(idxCMap(:,1) == iColor,iPC)./scaleFact1).*maxCoefLen,...
                (scores(idxCMap(:,1) == iColor,iPC+1)./scaleFact1).*maxCoefLen,50,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
        end
    end
    
    xlabel({['PC' num2str(iPC)] ; ['Percent Variance Explained ' num2str(exp(iPC),3) '%']});
    ylabel({['PC' num2str(iPC+1)] ; ['Percent Variance Explained ' num2str(exp(iPC+1),3) '%']});
    
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.fig']);
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.eps'],'psc2');
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.png']);
    
    hold on
    obsNames =  cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0);
    
    arrayfun(@(x) text((scores(x,iPC)./scaleFact1).*maxCoefLen,(scores(x,iPC+1)./scaleFact1).*maxCoefLen,obsNames{x}),1:length(scores(:,1)));
    
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.fig']);
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.eps'],'psc2');
    saveas(gcf,['PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.png']);
        
end % for iPC

%% 3D data plots 
    %% Biplot 
    fsFigure(0.75)
    biplot(coef(:,1:3),'varlabels',varNamesString);
    xlabel({'PC1' ; ['Percent Variance Explained ' num2str(exp(1),3) '%']});
    ylabel({'PC2' ; ['Percent Variance Explained ' num2str(exp(2),3) '%']});
    zlabel({'PC3' ; ['Percent Variance Explained ' num2str(exp(3),3) '%']});

    saveas(gcf,'PCABiPlotFirst3PCs.fig');
    saveas(gcf,'PCABiPlotFirst3PCs.eps','psc2');

    %% Color by Group Condition
    scaleFact1 = max(max(abs(scores(:,1:3)))) ;
    coefs = coef(:,1:3);
    maxCoefLen = sqrt(max(sum(coefs.^2,2)));

    hold on
    arrayfun(@(x) scatter3((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
        (scores(grouping==x,2)./scaleFact1).*maxCoefLen,(scores(grouping==x,3)./scaleFact1).*maxCoefLen,...
        50,toPlot.info.color{x},'filled'),1:nGroups);

    xlabel({'PC1' ; ['Percent Variance Explained ' num2str(exp(1),3) '%']});
    ylabel({'PC2' ; ['Percent Variance Explained ' num2str(exp(2),3) '%']});
    zlabel({'PC3' ; ['Percent Variance Explained ' num2str(exp(3),3) '%']});


    view(45,45);
    saveas(gcf,'PCABiPlotFirst3PCs.fig');
    saveas(gcf,'PCABiPlotFirst3PCs.eps','psc2');
    close gcf

%% Discrimination Metrics 
nCond = numel(dataMatAllMeas); 
[DB,Dunn,SC] = arrayfun(@(x) whDiscriminationMeasures(dataMatAllMeas{x}',dataMatAllMeas{1}'),1:nCond); 

names = toPlot.info.names';  
forCell  = num2cell([DB' Dunn']); 
values = [names forCell]; 
% cell2dataset(values); 
discrimValues =cell2table(values); 
save('Table_Dunn_DB','discrimValues'); 

% By PCA 

  scoresPC = arrayfun(@(x) scores(grouping==x,1:forLine),1:nGroups,'uniformoutput',0);
[DBPCA,DunnPCA] = arrayfun(@(x) whDiscriminationMeasures(scoresPC{x}',scoresPC{1}'),1:nCond); 
forCellPCA  = num2cell([DBPCA' DunnPCA']); 
valuesPCA = [names forCellPCA]; 
% cell2dataset(values); 
discrimValuesPCA =cell2table(valuesPCA); 
save('Table_Dunn_DB_PCA','discrimValuesPCA','scoresPC'); 



%% 
