function [ output_args ] =GCAGroupAnalysisScreenCorrelationsSingleGroup(toPlotGroup,varargin)
%% GCAGroupAnalysisScreenCorrelations 
% Takes the data in a toPlot structure 
% Scans for simple linear correlations between descriptor variable and
% response variable (this case outgrowth) for each group and with the
% total : 
% 
%
% 
% note was GCAGroupAnalysisScreenCorrelationsTest until 20160404
% not yet updated to handle the mutliple group plots 
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
ip.addParameter('ColorByNeuriteElongation',false); 
ip.addParameter('ColorMap',[]); 
ip.addParameter('Boxplots',true); 
ip.addParameter('perNeuriteStatistic','nanmean'); 
ip.addParameter('corrType','Pearson'); 

ip.parse(varargin{:});
figureDirectory = [ip.Results.OutputDirectory filesep ip.Results.corrType]; 
if ~isdir(figureDirectory)
    mkdir(figureDirectory)
end 
screenAll = [figureDirectory filesep 'ScreenAllCorr']; 
if ~isdir(screenAll) 
    mkdir(screenAll) 
end 

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

    toPlotGroup.info.projList = toPlotGroup.info.projList(idxCombine);
    toPlotGroup.info.names = toPlotGroup.info.names(idxCombine); 
%% ColorMap Initiation 
if isempty(ip.Results.ColorMap) 
    colorMap = brewermap(128,'spectral'); 
    colorMap = flip(colorMap,1); 
else 
    colorMap = ip.Results.ColorMap; 
end 
toPlotGroup = helperAddYLabels(toPlotGroup);  
%% Get all parameters 
for iGroup = 1:nGroups;
    
    projListC = vertcat(toPlotGroup.info.projList{:});

    projListC = projListC(:,1); 
    neuriteLengthStruct = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' ...
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteLengthOutput.mat']),projListC);
    
    
    neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0);
    deltasPerGroup{iGroup} = cellfun(@(x) (x(end)-x(1)),neuriteLengths);
end
  deltas = vertcat(deltasPerGroup{:}); 
  deltasPerGroup = deltasPerGroup'; 
  
  if ip.Results.ColorByNeuriteElongation
    mapper=linspace(0,2.5,128)';
    D=createDistanceMatrix(deltas./10,mapper);
    [sD,idxCMap]=sort(abs(D),2);
  end 
  
  % group colors 
  
  colors = toPlotGroup.info.color(idxCombine)'; 
  % 
  names = fieldnames(toPlotGroup);
  names = names(cellfun(@(x) ~strcmpi(x,'info'),names)); 
  for i = 1:numel(names)
      toPlotGroup.(names{i}).dataMat = toPlotGroup.(names{i}).dataMat(idxCombine);     
  end 
%% Screen for Correlations With Neurite Outgrowth (still a bit rough...)
[results,dataSetArray,forDataSetArray,rAll,pAll] = GCAAnalysisScreenCorrelations(deltas,toPlotGroup,'perNeuriteStatistic',ip.Results.perNeuriteStatistic,'type',ip.Results.corrType);

fsFigure(0.75); 
cMap = brewermap(128,'Rdbu'); 
cMap = flipud(cMap); 

%cMap = [[1,1,1];cMap(end,1),]; 

% set off diagonal to NaN
n = size(rAll,1); 
% for i=1:n
%     rAll(i,i:n) = NaN;
%    
%     
% end 
toPlotGroup = helperAddYLabels(toPlotGroup);
% Labels = 
imagesc(rAll,[-1,1]); % 
% for each row (ie each measurement) find pValue that is significant 
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
% set(gca, 'XTickLabel', ); % set x-axis labels
% set(gca, 'YTickLabel', L); % set y-axis labels
% get the 

L = get(dataSetArray,'VarNames'); 
LA = cellfun(@(x) strrep(x,'_',' '),L,'uniformoutput',0); 
set(gca, 'XTickLabel',LA ); % set x-axis labels
set(gca, 'YTickLabel', LA); % set y-axis labels

set(gca, 'XTickLabelRotation',45); 
title([ip.Results.corrType 'Correlation Coefficient Matrix']); 
colormap(cMap); 
colorbar

% find significant pairwise correlations 


saveas(gcf,[figureDirectory filesep ip.Results.corrType 'CorrelationMatrix' ip.Results.perNeuriteStatistic 'PerNeurite' '.fig']); 
saveas(gcf,[figureDirectory filesep ip.Results.corrType 'CorrelationMatrix' ip.Results.perNeuriteStatistic 'PerNeurite' '.eps'],'psc2'); 
saveas(gcf,[figureDirectory filesep ip.Results.corrType 'CorrelationMatrix' ip.Results.perNeuriteStatistic 'PerNeurite' '.png']); 
close gcf
% save the dataSetarray 
export(dataSetArray,'File',[figureDirectory filesep  'DataSet' ip.Results.perNeuriteStatistic '.csv'],'Delimiter',',');

 

% screen for the significant correlations and plot
for iVar = 1:n
    cDir = [screenAll filesep  L{iVar}];
    cDirName = L{iVar};
    if ~isdir(cDir)
        mkdir(cDir)
    end
    % make the directory
    %currentDir =
    
    % check for hits
    
    varToPlot = L(pAll(iVar,:)<0.05);
    idx = find(pAll(iVar,:)<0.05);
    if ~isempty(varToPlot)
        for iHit = 1:length(idx)
            % if the hits aren't empty make plot and save
            setAxis('on')
            
            hitName = varToPlot{iHit};
            
            % CurrentDir vs
            scatter(forDataSetArray(:,idx(iHit)),forDataSetArray(:,iVar),100,'k','filled');
            if isfield(toPlotGroup,hitName);
                
                xlabel(toPlotGroup.(hitName).yLabel);
            else
                xlabel({'Neurite Net Elongation Velocity in 10 min' ;'(um/min)'});
            end
            
            if isfield(toPlotGroup,cDirName);
                ylabel(toPlotGroup.(cDirName).yLabel);
            else
                ylabel({'Neurite Net Elongation Velocity in 10 min' ; '(um/min)'}); 
            end
              
            
            rC = rAll(iVar,idx(iHit)); 
            title({[ip.Results.corrType ' r  = ' num2str(rC,3)] ; ['PValue = ' num2str(pAll(iVar,idx(iHit)),3)]});  
            saveas(gcf,[cDir filesep hitName 'vs' cDirName ip.Results.corrType '.fig']);
            saveas(gcf,[cDir filesep hitName 'vs' cDirName ip.Results.corrType '.eps'],'psc2'); 
            saveas(gcf,[cDir filesep hitName 'vs' cDirName ip.Results.corrType '.png']); 
            
            close gcf
        end % iHit
        
    end % varToPlot
    
end % iVar
          
   
            
            
            
          
     


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

colorDir = [figureDirectory filesep 'CorrNeuriteElong_Color']; 
if ~isdir(colorDir) 
    mkdir(colorDir) 
end 

sigCorr = [colorDir  filesep 'SigCorr']; 
if ~isdir(sigCorr) 
    mkdir(sigCorr); 
end 

% Make the Correlation Plots
for iParam =1: length(idxHit)
    
   setAxis('on'); 
    %  subplot(2,2,1:2);
    resultsC = results(idxHit(iParam));
    hold on
    
    if ip.Results.ColorByNeuriteElongation
        outgrowthCMap = zeros(length(vertcat(deltasPerGroup{:})),3); 
        count = 1; 
        for k=1:128
            ptC  = find(idxCMap(:,1)==k);
            ptsXAll = vertcat(resultsC.v{:});
            ptsYAll = vertcat(deltasPerGroup{:});
            if ~isempty(ptC)
               
                scatter(ptsXAll(ptC),ptsYAll(ptC)./10,200,colorMap(k,:),'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',... 
                    colorMap(k,:));
                for i = 1:length(ptC)
                outgrowthCMap(count,:) = colorMap(k,:); % Quick and dirty
                count = count+1; 
                end 
            end
            hold on
        end
    else 
        cellfun(@(x,y,z) scatter(x',y./10,200,z,'filled'),resultsC.v,deltasPerGroup,colors);
        %     y = resultsC.values(:,1);
        %     x = resultsC.values(:,2);
        %     scatter(x,y,100,'filled');
        %     sorty = resultsC.values(idxSortGrowth,1); % put outgrowth on the y axis
        %
        %     sortx = resultsC.values(idxSortGrowth,2);
        %     arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
    end
    title({[ip.Results.corrType ' Correlation Coefficient ' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel({'Neurite Net Elongation Velocity in 10 min' ; '(um/min)'},'FontSize',14,'FontName','Arial');
    
    if ~isempty(regexpi(resultsC.name,'Veloc'));
%         axis([-0.5 8 -0.5 2.5]);   
     hold on 
     xValues = get(gcf,'XTick'); 
     yValues = get(gcf,'yTick'); 
     plot([min(xValues), max(xValues)], [min(yValues),max(yValues)],'Linewidth',2,'color','r'); 
    end
       xlabel(toPlotGroup.(resultsC.name).yLabel,'FontSize',14,'FontName','Arial');
    %xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[sigCorr filesep resultsC.name ip.Results.corrType '.fig']);
    saveas(gcf,[sigCorr filesep resultsC.name ip.Results.corrType '.eps'],'psc2');
    saveas(gcf,[sigCorr filesep resultsC.name ip.Results.corrType '.png']); 
    close gcf
    %'Labels',labels,'labelorientation','inline'
    if ip.Results.Boxplots
         
   
        if numel(toPlotGroup.info.names)==1 % currently only set up for 1 group
            labels = toPlotGroup.info.projList{1}(:,2);
            labels = cellfun(@(x) strrep(x,'_',' '), labels,'uniformoutput',0); 
            
            setAxis('on'); 
            [sortedDeltas,idxSort] =  sort(deltasPerGroup{iGroup});
            labels = labels(idxSort); 
            dataMatC = toPlotGroup.(resultsC.name).dataMat{1};
            dataMatSort = dataMatC(:,idxSort);
            cGroup = 1:length(sortedDeltas);
            perNeuriteStatistic = str2func(ip.Results.perNeuriteStatistic); 
            values = arrayfun(@(x) perNeuriteStatistic(dataMatSort(:,x)),1:size(dataMatSort,2)); 
            h1 =  boxplot(dataMatSort,'colorGroup',cGroup,...
                'colors',outgrowthCMap,'notch','on','outlierSize',1,...
                'symbol','+',...
                'orientation','horizontal','labels',labels);
            hold on 
            arrayfun(@(x) scatter(values(x),x,50,outgrowthCMap(x,:),'filled') ,1:length(values)); 
            set(h1(:),'Linewidth',2);
             %xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
             xlabel(toPlotGroup.(resultsC.name).yLabel,'FontSize',14,'FontName','Arial'); 
             axis([nanmin(dataMatC(:)),nanmax(dataMatC(:)),0.5,length(sortedDeltas)+0.5]); 
           set(gca,'FontName','Arial','FontSize',12);
        end
    saveas(gcf,[sigCorr filesep resultsC.name 'Boxplot.fig']);
    saveas(gcf,[sigCorr filesep resultsC.name 'Boxplot.eps'],'psc2');
    saveas(gcf,[sigCorr filesep resultsC.name 'Boxplot.png']); 
    close gcf
    end
    
 
end

%% BorderLine Correlations With Outgrowh
bordDir = [colorDir filesep 'BorderCorr']; 
if ~isdir(bordDir) 
    mkdir(bordDir); 
end 

for iParam =1: length(idxBorderline)
     
    setAxis('on')
   
    % subplot(2,2,1:2);
%     y = resultsC.values(:,1); 
%     x = resultsC.values(:,2); 
    resultsC = results(idxBorderline(iParam));
    
    if ip.Results.ColorByNeuriteElongation
        
        for k=1:128
            ptC  = find(idxCMap(:,1)==k);
            ptsXAll = vertcat(resultsC.v{:});
            ptsYAll = vertcat(deltasPerGroup{:});
            if ~isempty(ptC)
                scatter(ptsXAll(ptC),ptsYAll(ptC)./10,200,colorMap(k,:),'filled');
            end
            hold on
        end
    else
        
        %     scatter(x,y,100,'filled');
        cellfun(@(x,y,z) scatter(x',y./10,200,z,'filled'),resultsC.v,deltasPerGroup,colors);
    end
    
    
    title({[ip.Results.corrType  ' Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel({'Neurite Net Elongation Velocity in 10 min' ; '(um/min)'},'FontSize',14,'FontName','Arial');
%     xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
xlabel(toPlotGroup.(resultsC.name).yLabel,'FontSize', 14,'FontName','Arial');
    
    if ~isempty(regexpi(resultsC.name,'Veloc'));
        %axis([-0.5 8 -0.5 2.5]);
        hold on
        xValues = get(gca,'xTick');
        yValues = get(gca,'yTick');
        plot(xValues,yValues,'r');
    end
    
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[bordDir filesep resultsC.name ip.Results.corrType '.fig']);
    saveas(gcf,[bordDir filesep resultsC.name ip.Results.corrType '.eps'],'psc2');
    saveas(gcf,[bordDir filesep resultsC.name ip.Results.corrType '.png']); 
    
    

    close gcf
   
end
%%  No Correlation with Outgrowth

noCorrDir = [colorDir filesep 'noCorr']; 
if ~isdir(noCorrDir) 
    mkdir(noCorrDir)
end 

for iParam = 1:length(idxNoCorr)
    setAxis('on')
    % subplot(2,2,1:2);
    resultsC = results(idxNoCorr(iParam));
    
%    y = resultsC.values(:,1);
%     
%    x = resultsC.values(:,2);
    hold on
    %arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
%     scatter(x,y,100,'filled'); 


 if ip.Results.ColorByNeuriteElongation
        
        for k=1:128
            ptC  = find(idxCMap(:,1)==k);
            ptsXAll = vertcat(resultsC.v{:});
            ptsYAll = vertcat(deltasPerGroup{:});
            if ~isempty(ptC)
                scatter(ptsXAll(ptC),ptsYAll(ptC)./10,200,colorMap(k,:),'filled');
            end
            hold on
        end
    else


    

    cellfun(@(x,y,z) scatter(x',y./10,200,z,'filled'),resultsC.v,deltasPerGroup,colors); 
    
 end
         

    title({[ip.Results.corrType 'Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel({'Neurite Net Elongation Velocity in 10 min' ;'(um/min)'},'FontSize',14,'FontName','Arial');
%     xlabel(strrep(resultsC.name,'_',' '),'FontSize', 14,'FontName','Arial');
 xlabel(toPlotGroup.(resultsC.name).yLabel,'FontSize',14,'FontName','Arial'); 
    
    set(gca,'FontName','Arial','FontSize',12);
    
    
    saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.fig']);
    saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.eps'],'psc2');
    saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.png']); 
    
    if ~isempty(regexpi(resultsC.name,'Veloc'));
        hold on
       
        xValues = get(gca,'xTick');
        yValues = get(gca,'yTick');
        % make a line showing perfect 1 to 1 correlation with outgrowth
        % velocity 
        line([0,10],[0,10],'color','r','linewidth',2); 
        % reset the axis 
        %axis([min(xValues),max(xValues),min(yValues),max(yValues)]); 
        
        saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.fig']);
        saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.eps'],'psc2');
        saveas(gcf,[noCorrDir filesep resultsC.name ip.Results.corrType '.png']);
               
    end
    
      
    
    close gcf
   
end





end

