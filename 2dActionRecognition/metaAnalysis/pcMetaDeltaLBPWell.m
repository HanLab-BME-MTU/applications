%% At the well level:
% Assaf Zaritsky, June. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaDeltaLBPWell()

always = false;

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/dLBP/'];

outDname = [lbpDirname filesep 'FOV/'];
accLbpPrefix = [lbpDirname filesep 'accumulatedDeltaLBP_'];

if ~exist(outDname,'dir')
    unix(sprintf('mkdir %s',outDname));
end

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

bins = 0.025:0.05:1-0.025; % this is for the distributions

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
    
    accLbpFnameAll = [accLbpPrefix num2str(iScale) '_' testStr '_all.mat'];
    accLbpFnameSource = [accLbpPrefix num2str(iScale) '_' testStr '_source.mat'];
    accLbpFnameMetastatic = [accLbpPrefix num2str(iScale) '_' testStr '_metastatic.mat'];
    accLbpFnameCellType = [accLbpPrefix num2str(iScale) '_' testStr '_type.mat'];
    
    %% all
    allFeatsFname = [outDnameScale filesep 'allFeats.mat'];
    
    loggerFname = [outDnameScale filesep 'log_' num2str(iScale) '_' testStr '_all.txt'];    
    
    % take all cells to variable named "allCells"   
    load([accLbpPrefix num2str(iScale) '_allInfo.mat']); % allInfo
    
    load(accLbpFnameAll); % these accumulations were created by pcAccumulateDeltaLBPWell   
    allCells = allCellsFov;
    
    dataAll = allCells.accDeltaLbp;
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
        
        allDeltaLBP = getStats(dataAll); % an array that includes all median's dLBPs
        lowTH = prctile(allDeltaLBP,1);
        highTH = prctile(allDeltaLBP,97);
        
        [distributionsAll,distributionsLocations,meansAll,meansLocations,strLocations] = getDeltaLbpDistributions(dataAll,dataLocations,strAll,bins,[lowTH,highTH],logger,outDnameScaleAll);
        
        fclose(logger);
        
        save(allFeatsFname,'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allDeltaLBP','lowTH','highTH','n');
        fprintf(sprintf('%s: done all: %d %s\n',testStr,iScale,'all'));
        close all;
    end
    
    load(allFeatsFname);%'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allDeltaLBP','lowTH','highTH','n');
    plotDeltaLbpDistributionFromCellArray(dataAll,bins,lowTH,highTH,outDnameScale,'all');
    
    %% Source: cell line / primary melanoma / melanocyte
    sourceFeatsFname = [outDnameScale filesep 'sourceFeats.mat'];
    
    if ~exist(sourceFeatsFname,'file') || always
        load(accLbpFnameSource);
        
        if exist('cellLinesFov','var')
            cellLines = cellLinesFov;
            melanocytes = melanocytesFov;
            tumors = tumorsFov;
            clear cellLinesFov melanocytesFov tumorsFov;
        end
        
        nCellLines = length(cellLines.accDeltaLbp);
        nMelanocytes = length(melanocytes.accDeltaLbp);
        nPrimaryMelanoma = length(tumors.accDeltaLbp);
        dataSource = [cellLines.accDeltaLbp,melanocytes.accDeltaLbp,tumors.accDeltaLbp];
        clear cellLines melanocytes tumors;
        labelsSource = {'CellLine','Melano','Primary'};
        indsSource = [nCellLines,nCellLines+nMelanocytes,nCellLines+nMelanocytes+nPrimaryMelanoma];
        outDirSource = outDnameScale;%[outDnameScale 'source/'];
        titleSource = 'source';
        
        doMetaDeltaLBP(dataSource,labelsSource,indsSource,outDirSource,bins,lowTH,highTH,titleSource);
        
        save(sourceFeatsFname,'dataSource','labelsSource','indsSource','outDirSource','bins','lowTH','highTH','titleSource');
        clear dataSource;
        fprintf(sprintf('%s: done source %d\n',testStr,iScale));
        close all;
    end
    
    %% High vs. Low metastatic efficiency
    metEffFeatsFname = [outDnameScale filesep 'metEffFeats.mat'];
    
    if ~exist(metEffFeatsFname,'file') || always
        load(accLbpFnameMetastatic);
        
        if exist('tumorHighFov','var')
            tumorHigh = tumorHighFov;
            tumorLow = tumorLowFov;
            clear tumorHighFov tumorLowFov;
        end
        
        nHighMet = length(tumorHigh.accDeltaLbp);
        nLowMet = length(tumorLow.accDeltaLbp);
        dataMetEff = [tumorHigh.accDeltaLbp,tumorLow.accDeltaLbp];
        clear tumorHigh tumorLow;
        labelsMetEff = {'High','Low'};
        indsMetEff = [nHighMet,nHighMet+nLowMet];
        outDirMetEff = outDnameScale;%[outDnameScale 'metEfficiency/'];
        titleMetEff = 'metEff';
        
        doMetaDeltaLBP(dataMetEff,labelsMetEff,indsMetEff,outDirMetEff,bins,lowTH,highTH,titleMetEff);
        
        save(metEffFeatsFname,'dataMetEff','labelsMetEff','indsMetEff','outDirMetEff','bins','lowTH','highTH','titleMetEff');
        clear dataMetEff;
        fprintf(sprintf('%s: done metastatic %d\n',testStr,iScale));
        close all;
    end
    
    %% Independent cell types
    cellTypeFeatsFname = [outDnameScale filesep 'cellType.mat'];
    
    if ~exist(cellTypeFeatsFname,'file') || always
        
        load(accLbpFnameCellType);
        
        if exist('cellTypesFov','var')
            cellTypes = cellTypesFov;
            clear cellTypesFov;
        end
        
        nCellTypes = length(cellTypes);
        
        dataCellType = {};
        labelsCellType = {};
        indsCellType = [];
        for i = 1 : nCellTypes
            if isempty(cellTypes{i}.accDeltaLbp)
                continue;
            end
            curI = length(labelsCellType) + 1;
            dataCellType = [dataCellType, cellTypes{i}.accDeltaLbp];
            labelsCellType{curI} = lower(cellTypesStr{i});
            if curI == 1
                curInd = 0;
            else
                curInd = indsCellType(end);
            end
            indsCellType = [indsCellType, (curInd+length(cellTypes{i}.accDeltaLbp))];
        end
        
        outDirCellType = [outDnameScale 'cellType/'];
        titleCellType = 'cellType';
        
        doMetaDeltaLBP(dataCellType,labelsCellType,indsCellType,outDirCellType,bins,lowTH,highTH,titleCellType);
        
        save(cellTypeFeatsFname,'dataCellType','labelsCellType','indsCellType','outDirCellType','bins','lowTH','highTH','titleCellType');
        clear dataCellType;
        fprintf(sprintf('%s: done cell type %d\n',testStr,iScale));
        close all;
    end
end
end


%%

%% goes over all experiments
function [distributionsAll,distributionsLocations,meansAll,meansLocations,strLocations] = getDeltaLbpDistributions(data,dataLocations,strAll,bins,extremeVals,logger,outDnameScaleAll)

n = length(data);

distributionsAll = cell(1,n);
meansAll = nan(1,n);
distributionsLocations = cell(1,n);
meansLocations = cell(1,n);
strLocations = cell(1,n);

for i = 1 : n
    [distExp,distLoc,meanExp,meanLocs] = getDeltaLbpFeatures(data{i},dataLocations{i},strAll{i},bins,extremeVals,logger,outDnameScaleAll);    
    
    distributionsAll{i} = distExp;
    distributionsLocations{i} = distLoc;
    strLocations{i} = dataLocations{i}.locationStr;
    meansAll(i) = meanExp;
    meansLocations{i} = meanLocs;
    close all;
end
end


%%
function allDeltaLBP = getStats(dataAll)
n = length(dataAll);
allDeltaLBP = [];

for i = 1 : n
    allDeltaLBP = [allDeltaLBP dataAll{i}];
end
end

%% Data from a single experiment (multiple locations)
%     locationDeltaLbp: {1x10 cell}
%          locationStr: {'01'  '02'  '03'  '04'  '05'  '06'  '07'  '08'  '09'  '10'}
function [distExp,distLoc,meanExp,meanLocs] = getDeltaLbpFeatures(data,dataLocations,curStr,bins,extremeVals,logger,outDnameScaleAll)

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

plotDeltaLbpDistribution(distExp,bins,[outDnameScaleAll curStr '.eps']);

fprintf(logger,'Locations: \n');
nLoc = length(dataLocations.locationDeltaLbp);
meanLocs = nan(1,nLoc);
for iloc = 1 : nLoc
    curLocData = dataLocations.locationDeltaLbp{iloc};
    curLocStr = dataLocations.locationStr{iloc};
                
    curLocDataFiltered = curLocData((curLocData >= lowTH) & (curLocData <= highTH));
    locDataNorm = (curLocDataFiltered - lowTH)./(highTH-lowTH);
    [nelements, xcenters] = hist(locDataNorm,bins);
    distLoc = nelements ./ sum(nelements);
    meanLocs(iloc) = mean(locDataNorm);
    
    locPrcntExcludeHigh = 100*sum(curLocData > highTH) / length(curLocData);
    
    fprintf(logger,sprintf('%s: mean = %.2f, n = %d\n',curLocStr,meanLocs(iloc),length(curLocData)));
    fprintf(logger,sprintf('excluded = %.1f%%%%: %s\n\n',locPrcntExcludeHigh,mat2str(find(curLocData > highTH))));
    
    plotDeltaLbpDistribution(distLoc,bins,[outDnameScaleLocations curStr '_s' curLocStr '.eps']);
end
end

%%
function [] = plotDeltaLbpDistribution(distribution,centers,outFname)
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

h = figure;
hold on;
bar(centers,distribution,'k');
xlabel('\Delta LBP','FontSize',fontsize);
ylabel('Frequency','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,1]);
set(haxes,'XTick',0:0.5:1);
set(haxes,'XTickLabel',0:0.5:1);
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
function [] = doMetaDeltaLBP(data,labels,inds,outDir,bins,lowTH,highTH,titleStr)

close all;

if ~exist(outDir,'dir')
    unix(sprintf('mkdir %s',outDir));
end

plotDeltaLbpDistributionFromCellArray(data,bins,lowTH,highTH,outDir,titleStr);

for i = 1 : length(labels)
    if i == 1
        istart = 1;
    else
        istart = inds(i-1) + 1;
    end
    iend = inds(i);
    
    curData = data(istart:iend);    
    plotDeltaLbpDistributionFromCellArray(curData,bins,lowTH,highTH,outDir,[titleStr '_' labels{i}]);            
end

close all;
% plotDeltaLbpDistribution(distribution,centers,outFname)

end

function [] = plotDeltaLbpDistributionFromCellArray(data,bins,lowTH,highTH,outDir,titleStr)
dataFlat = getStats(data);
dataFiltered = dataFlat((dataFlat >= lowTH) & (dataFlat <= highTH));
dataNorm = (dataFiltered - lowTH)./(highTH-lowTH);
[nelements, xcenters] = hist(dataNorm,bins);
distGroup = nelements ./ sum(nelements);
meanGroup = mean(dataNorm);

plotDeltaLbpDistribution(distGroup,bins,[outDir titleStr '.eps']);
end


