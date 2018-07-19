%% At the well level:
% Look at distributions of plasticity activities (delta LBP, 20 frames,
% pval = 0.01)
% Assaf Zaritsky, June. 2017
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMeta_Plasticity_LBP_t20_pval01_t120()

always = false;

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/Plasticity_LBP_t20_pval01_t120/'];

outDname = [lbpDirname filesep 'FOV/'];
accFeatsPrefix = [lbpDirname filesep 'Plasticity_LBP_t20_pval01_t120_'];

if ~exist(outDname,'dir')
    unix(sprintf('mkdir %s',outDname));
end

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

% bins = 0.025:0.05:1-0.025; % this is for the distributions
bins = 0.005:0.01:0.15; % this is for the distributions

testStr = 'fov';
for iScale = 1 : nScales % resolution (1- maximal)
    outDnameScale = [outDname filesep num2str(iScale) filesep];
    
    if ~exist(outDnameScale,'dir')
        unix(sprintf('mkdir %s',outDnameScale));
    end
    
    outDnameScaleAll = [outDnameScale filesep 'all' filesep];
    if ~exist(outDnameScaleAll,'dir')
        unix(sprintf('mkdir %s',outDnameScaleAll));
    end
    
    outDnameScaleLocations = [outDnameScaleAll 'Locations' filesep];
    if ~exist(outDnameScaleLocations,'dir')
        unix(sprintf('mkdir %s',outDnameScaleLocations));
    end
    
    accFeatsFnameAll = [accFeatsPrefix num2str(iScale) '_' testStr '_all.mat'];
    accFeatsFnameSource = [accFeatsPrefix num2str(iScale) '_' testStr '_source.mat'];
    accFeatsFnameMetastatic = [accFeatsPrefix num2str(iScale) '_' testStr '_metastatic.mat'];
    accFeatsFnameCellType = [accFeatsPrefix num2str(iScale) '_' testStr '_type.mat'];
    
    %% all
    allFeatsFname = [outDnameScale filesep 'allFeats.mat'];
    
    loggerFname = [outDnameScale filesep 'log_' num2str(iScale) '_' testStr '_all.txt'];    
    
    % take all cells to variable named "allCells"   
    load([accFeatsPrefix num2str(iScale) '_allInfo.mat']); % allInfo
    
    load(accFeatsFnameAll); % these accumulations were created by accLchFeats_LBP_Plasticity_t20_pval01_t120_All   
    allCells = allCellsFov;
    
    dataAll = allCells.accFeats;
    strAll = allInfo.strs;
    
    if ~exist(allFeatsFname,'file') || always                       
        logger = fopen(loggerFname,'w');
        
        dataLocations = allCells.locations; % has fields locationDeltaLbp and locationStr
        
        dateAll = allInfo.date;
        cellTypeAll = allInfo.cellType;
        cellTypeIndAll = allInfo.cellTypeInd; % index in metaData.cellTypes.ids
        sourceAll = allInfo.source;
        metEffAll = allInfo.metEff;
        
        clear allCells;
        n = length(dataAll);
        
        allFeats = getStats(dataAll); % an array that includes all median's dLBPs
        lowTH = prctile(allFeats,1);
        highTH = prctile(allFeats,97);
        
        [distributionsAll,distributionsLocations,meansAll,meansLocations,strLocations] = getDistributions(dataAll,dataLocations,strAll,bins,[lowTH,highTH],logger,outDnameScaleAll);
        
        fclose(logger);
        
        save(allFeatsFname,'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allFeats','lowTH','highTH','n');
        fprintf(sprintf('%s: done all: %d %s\n',testStr,iScale,'all'));
        close all;
    end
    
    load(allFeatsFname);%'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allFeats','lowTH','highTH','n');
    plotDistributionFromCellArray(dataAll,bins,lowTH,highTH,outDnameScale,'all');
    
    %% Source: cell line / primary melanoma / melanocyte
    
    cellLines = filterData(dataAll,allInfo.source,'CellLines');
    melanocytes = filterData(dataAll,allInfo.source,'Melanocytes');    
    tumors = filterData(dataAll,allInfo.source,'Tumors');
    
    fprintf(sprintf('All: %.4f\n',mean([cellLines.data,melanocytes.data,tumors.data])));
    fprintf(sprintf('Cell Lines: %.4f; Melanocytes: %.4f; Tumor: %.4f\n',...
        mean(cellLines.data),mean(melanocytes.data),mean(tumors.data)));
    
    dataSource = [cellLines.data,melanocytes.data,tumors.data];
    labelsSource = {'CellLine','Melano','Primary'};
    indsSource = [cellLines.n,cellLines.n+melanocytes.n,cellLines.n+melanocytes.n+tumors.n];
    outDirSource = outDnameScale;%[outDnameScale 'source/'];
    titleSource = 'source';
        
    doPlotFeatsDistribution(dataSource,labelsSource,indsSource,outDirSource,bins,lowTH,highTH,titleSource);
    
    sourceFeatsFname = [outDnameScale filesep 'sourceFeats.mat'];
    save(sourceFeatsFname,'dataSource','labelsSource','indsSource','outDirSource','bins','lowTH','highTH','titleSource');
    clear dataSource;
    fprintf(sprintf('%s: done source %d\n',testStr,iScale));
    close all;
    
    
    %% High vs. Low metastatic efficiency        
    highMet = filterData(dataAll,allInfo.metEff,1);
    lowMet = filterData(dataAll,allInfo.metEff,0);
    
    fprintf(sprintf('highMet: %.4f; lowMet: %.4f\n',...
        mean(highMet.data),mean(lowMet.data)));
    
    dataMetEff = [highMet.data,lowMet.data];
    labelsMetEff = {'High','Low'};
    indsMetEff = [highMet.n,highMet.n+lowMet.n];
    outDirMetEff = outDnameScale;%[outDnameScale 'metEfficiency/'];
    titleMetEff = 'metEff';
    
    doPlotFeatsDistribution(dataMetEff,labelsMetEff,indsMetEff,outDirMetEff,bins,lowTH,highTH,titleMetEff);
    
    metEffFeatsFname = [outDnameScale filesep 'metEffFeats.mat'];
    save(metEffFeatsFname,'dataMetEff','labelsMetEff','indsMetEff','outDirMetEff','bins','lowTH','highTH','titleMetEff');
    clear dataMetEff;
    fprintf(sprintf('%s: done metastatic %d\n',testStr,iScale));
    close all;
    
    %% Independent cell types        
    cellTypes = unique(allInfo.cellType);
    nCellTypes = length(cellTypes);
    
    outDirCellType = [outDnameScale 'cellType/'];
    titleCellType = 'cellType';
    
    dataCellType = [];
    labelsCellType = {};
    indsCellType = [];
    for ictype = 1 : nCellTypes
        curCellTypeStr = cellTypes{ictype};
        curCellType = filterData(dataAll,allInfo.cellType,curCellTypeStr);
        
        fprintf(sprintf('\n%s: %.4f\n',curCellTypeStr,mean(curCellType.data)));
        
        dataCellType = [dataCellType curCellType.data];
        labelsCellType{length(labelsCellType)+1} = lower(curCellTypeStr);
        if isempty(indsCellType)
            indsCellType = [curCellType.n];
        else
            indsCellType = [indsCellType, indsCellType(end)+curCellType.n];
        end
    end          
    
    cellTypeFeatsFname = [outDnameScale filesep 'cellType.mat'];
    doPlotFeatsDistribution(dataCellType,labelsCellType,indsCellType,outDirCellType,bins,lowTH,highTH,titleCellType);
    
    save(cellTypeFeatsFname,'dataCellType','labelsCellType','indsCellType','outDirCellType','bins','lowTH','highTH','titleCellType');
    clear dataCellType;
    fprintf(sprintf('%s: done cell type %d\n',testStr,iScale));
    close all;    
end
end


%%

%% goes over all experiments
function [distributionsAll,distributionsLocations,meansAll,meansLocations,strLocations] = getDistributions(data,dataLocations,strAll,bins,extremeVals,logger,outDnameScaleAll)

n = length(data);

distributionsAll = cell(1,n);
meansAll = nan(1,n);
distributionsLocations = cell(1,n);
meansLocations = cell(1,n);
strLocations = cell(1,n);

for i = 1 : n
    [distExp,distLoc,meanExp,meanLocs] = getFeatures(data{i},dataLocations{i},strAll{i},bins,extremeVals,logger,outDnameScaleAll);    
    
    distributionsAll{i} = distExp;
    distributionsLocations{i} = distLoc;
    strLocations{i} = dataLocations{i}.locationStr;
    meansAll(i) = meanExp;
    meansLocations{i} = meanLocs;
    close all;
end
end


%%
function allFeats = getStats(dataAll)
n = length(dataAll);
allFeats = [];

for i = 1 : n
    allFeats = [allFeats dataAll{i}];
end
end

%% Data from a single experiment (multiple locations)
%     locationDeltaLbp: {1x10 cell}
%          locationStr: {'01'  '02'  '03'  '04'  '05'  '06'  '07'  '08'  '09'  '10'}
function [distExp,distLoc,meanExp,meanLocs] = getFeatures(data,dataLocations,curStr,bins,extremeVals,logger,outDnameScaleAll)

outDnameScaleLocations = [outDnameScaleAll 'Locations' filesep];

lowTH = extremeVals(1);    
highTH = extremeVals(2);

dataFiltered = data((data >= lowTH) & (data <= highTH));
dataNorm = (dataFiltered - lowTH)./(highTH-lowTH);
[nelements, xcenters] = hist(dataNorm,bins);
distExp = nelements ./ sum(nelements);
meanExp = mean(dataNorm);

prcntExcludeHigh = 100*sum(data > highTH) / length(data);

fprintf(logger,sprintf('\n\n -------- %s --------- \n',curStr));
fprintf(logger,sprintf('Mean delta LBP: %.2f\n',meanExp));
tmpStr = sprintf('Excluded cells: %.1f%%%% (n = %d)\n',prcntExcludeHigh,length(data));
fprintf(logger,tmpStr);
tmpStr = sprintf('Distribution: %s\n',mat2str(round(distExp*100)));
fprintf(logger,tmpStr);
fprintf(logger,sprintf('Accumulated distribution: %s\n\n\n',mat2str(round(cumsum(distExp*100)))));

plotFeatsDistribution(distExp,bins,[outDnameScaleAll curStr '.eps']);

fprintf(logger,'Locations: \n');
nLoc = length(dataLocations.locationFeats);
meanLocs = nan(1,nLoc);
for iloc = 1 : nLoc
    curLocData = dataLocations.locationFeats{iloc};
    curLocStr = dataLocations.locationStr{iloc};
                
    curLocDataFiltered = curLocData((curLocData >= lowTH) & (curLocData <= highTH));
    locDataNorm = (curLocDataFiltered - lowTH)./(highTH-lowTH);
    [nelements, xcenters] = hist(locDataNorm,bins);
    distLoc = nelements ./ sum(nelements);
    meanLocs(iloc) = mean(locDataNorm);
    
    locPrcntExcludeHigh = 100*sum(curLocData > highTH) / length(curLocData);
    
    fprintf(logger,sprintf('%s: mean = %.2f, n = %d\n',curLocStr,meanLocs(iloc),length(curLocData)));
    fprintf(logger,sprintf('excluded = %.1f%%%%: %s\n\n',locPrcntExcludeHigh,mat2str(find(curLocData > highTH))));
    
    plotFeatsDistribution(distLoc,bins,[outDnameScaleLocations curStr '_s' curLocStr '.eps']);
end
end

%%
function [] = plotFeatsDistribution(distribution,centers,outFname)
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

h = figure;
hold on;
bar(centers,distribution,'k');
xlabel('\Delta LBP','FontSize',fontsize);
ylabel('Frequency','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,0.15]);
set(haxes,'XTick',0:0.05:0.15);
set(haxes,'XTickLabel',0:0.05:0.15);
set(haxes,'YLim',[0,0.2]);
set(haxes,'YTick',0:0.1:0.2);
set(haxes,'YTickLabel',0:0.1:0.2);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig(outFname);
end

%%
function [] = doPlotFeatsDistribution(data,labels,inds,outDir,bins,lowTH,highTH,titleStr)

close all;

if ~exist(outDir,'dir')
    unix(sprintf('mkdir %s',outDir));
end

plotDistributionFromCellArray(data,bins,lowTH,highTH,outDir,titleStr);

for i = 1 : length(labels)
    if i == 1
        istart = 1;
    else
        istart = inds(i-1) + 1;
    end
    iend = inds(i);
    
    curData = data(istart:iend);    
    plotDistributionFromCellArray(curData,bins,lowTH,highTH,outDir,[titleStr '_' labels{i}]);            
end

close all;
% plotFeatsDistribution(distribution,centers,outFname)

end

function [] = plotDistributionFromCellArray(data,bins,lowTH,highTH,outDir,titleStr)
if iscell(data)
    dataFlat = getStats(data);
else
    dataFlat = data;
end
% dataFiltered = dataFlat((dataFlat >= lowTH) & (dataFlat <= highTH));
% dataNorm = (dataFiltered - lowTH)./(highTH-lowTH);
[nelements, xcenters] = hist(dataFlat,bins);
distGroup = nelements ./ sum(nelements);
meanGroup = mean(dataFlat);

plotFeatsDistribution(distGroup,bins,[outDir titleStr '.eps']);
end

function out = filterData(data,labels,filterLabel)

n = length(data);
out.data = [];
for i = 1 : n    
    if (isscalar(filterLabel) && filterLabel == labels{i}) || (strcmpi(filterLabel,labels{i}))
        out.data = [out.data data{i}];
    end
end
out.n = length(out.data);
end


