function [ output_args ] =GCAGroupAnalysisScreenCorrelations(toPlotGroup,varargin)
%% GCAGroupAnalysisScreenCorrelations 
% Takes the data in a toPlot structure 
% Scans for simple linear correlations between descriptor variable and
% response variable (this case outgrowth) for each group and with the
% total 
%
% 

%% INPUT 
% toPlotGroup : toPlot structure WITH FIELDS WANT TO CORRELATE ALREADY
%               COLLECTED
%% Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = pwd;

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('Interactive',true); 
ip.addParameter('perNeuriteStatistic','nanmean'); 
ip.addParameter('matrixPlot',false); 

ip.parse(varargin{:});
figureDirectory = ip.Results.OutputDirectory; 
%% 

if ip.Results.Interactive
    groupNames = toPlotGroup.info.names;
    idxCombine = listSelectGUI(groupNames,[],'move');
else
    projListC = vertcat(toPlotGroup.info.projList{:});
    idxCombine = 1:numel(toPlotGroup.info.names);
end
nGroups = length(idxCombine);
%projListC = toPlotGroup.info.projList{1}(:,1);
%% Get all parameters 
for iGroup = 1:nGroups;
    
    projListC = toPlotGroup.info.projList{idxCombine(iGroup)};
    %toPlotGroup.info.projList = toPlotGroup.info.projList(idxCombine(iGroup)); 
    projListC = projListC(:,1); 
    neuriteLengthStruct = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' ...
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteLengthOutput.mat']),projListC);
    
    
    neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0);
    deltasPerGroup{iGroup} = cellfun(@(x) (x(end)-x(1)),neuriteLengths);
end
deltas = vertcat(deltasPerGroup{:}); 
  deltasPerGroup = deltasPerGroup'; 
  colors = toPlotGroup.info.color(idxCombine)'; 
  % 
  names = fieldnames(toPlotGroup);
  names = cellfun(@(x) ~strcmpi(x,'info'),names,'uniformoutput',0); 
  for i = 1:numel(names)
   
      
  end 
%% Screen for Correlations With Neurite Outgrowth (still a bit rough...)
[results,dataSetArray] = GCAAnalysisScreenCorrelations(deltas,toPlotGroup,'perNeuriteStatistic',...
    ip.Results.perNeuriteStatistic,'matrixPlot',ip.Results.matrixPlot);

% save the dataSetarray 
export(dataSetArray,'File',[figureDirectory filesep  'DataSet' ip.Results.perNeuriteStatistic 'WholeMovie.csv'],'Delimiter',',');

mdl = fitlm(dataSetArray); 
save([figureDirectory filesep 'ModelObjectWholeMovie.mat'],'mdl'); 
save([figureDirectory filesep 'DataSetArray.mat'],'dataSetArray'); 

% corrplotMine(dataSetArray,'testR','on'); 
% saveas(gcf,[figureDirectory filesep 'CorrelationsAll.fig']); 
% close gcf

% results of the correlations will be a simple %
%if ~isempty(resultsCorrScreen.Hit{1,1});
% Get the Indices of the Hits
idxHit = find(arrayfun(@(x) results(x).p(1,2) <0.05,1:length(results)));
idxBorderline = find(arrayfun(@(x) (results(x).p(1,2) >0.05 & results(x).p(1,2) < 0.1 ),1:length(results)));
idxNoCorr = find(arrayfun(@(x) (results(x).p(1,2) > 0.1),1:length(results)));
%% Potentially Significant Correlation With Outgrowth

sigCorr = [figureDirectory filesep 'SigCorr']; 
if ~isdir(sigCorr) 
    mkdir(sigCorr); 
end 

% Make the Correlation Plots
for iParam =1: length(idxHit)
    
    fsFigure(1);
    %  subplot(2,2,1:2);
    resultsC = results(idxHit(iParam));
    hold on
    
  
    cellfun(@(x,y,z) scatter(x',y./10,50,z,'filled'),resultsC.v,deltasPerGroup,colors); 
%     y = resultsC.values(:,1); 
%     x = resultsC.values(:,2); 
%     scatter(x,y,100,'filled'); 
%     sorty = resultsC.values(idxSortGrowth,1); % put outgrowth on the y axis
%     
%     sortx = resultsC.values(idxSortGrowth,2);
%     arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity (um/min) in 10 min','FontSize',14,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.png']); 
    close gcf
    fsFigure(1);
    %subplot(2,2,3:4);
    % sort by outgrowth
    
    
    %dataMat = toPlotGroup.(resultsC.name).dataMat{1};

   
end

%% BorderLine Correlations With Outgrowh
bordDir = [figureDirectory filesep 'BorderLineCorr']; 
if ~isdir(bordDir) 
    mkdir(bordDir); 
end 

for iParam =1: length(idxBorderline)
    
    fsFigure(1);
    hold on
    % subplot(2,2,1:2);
%     y = resultsC.values(:,1); 
%     x = resultsC.values(:,2); 
    resultsC = results(idxBorderline(iParam));
%     scatter(x,y,100,'filled'); 
cellfun(@(x,y,z) scatter(x',y./10,50,z,'filled'),resultsC.v,deltasPerGroup,colors); 
    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity (um/min) in 10 min','FontSize',14,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.png']); 
    
    

    close gcf
   
end
%%  No Correlation with Outgrowth

noCorrDir = [figureDirectory filesep 'noCorrelation']; 
if ~isdir(noCorrDir) 
    mkdir(noCorrDir)
end 

for iParam = 1:length(idxNoCorr)
    fsFigure(1);
    % subplot(2,2,1:2);
    resultsC = results(idxNoCorr(iParam));
    
%    y = resultsC.values(:,1);
%     
%    x = resultsC.values(:,2);
    hold on
    %arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
%     scatter(x,y,100,'filled'); 
cellfun(@(x,y,z) scatter(x',y./10,50,z,'filled'),resultsC.v,deltasPerGroup,colors); 

    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity (um/min) in 10 min','FontSize',14,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.png']); 
    
    close gcf
   
end

if ip.Results.matrixPlot 
    corrPlotMine(    ); 
end 

end

