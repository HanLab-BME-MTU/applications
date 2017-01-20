%% Mixing labeling of high/low met efficiency and examining kymographs
% Assaf Zaritsky, July. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaDeltaLBPWell_mixMetEffLabels()

always = true;

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/dLBP/'];

outDnameFov = [lbpDirname filesep 'FOV/'];
accLbpPrefix = [lbpDirname filesep 'accumulatedLBP_'];

if ~exist(outDnameFov,'dir')
    unix(sprintf('mkdir %s',outDnameFov));
end

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

outDirs = {outDnameFov}; 
testStrs = {'fov'}; 
nTest = length(testStrs);

for iTest = 1 : nTest 
    outDname = outDirs{iTest};
    testStr = testStrs{iTest}; % field of view 
    for iScale = 1 : nScales % resolution (1- maximal)
        outDnameScale = [outDname filesep num2str(iScale) filesep];
        
        if ~exist(outDnameScale,'dir')
            unix(sprintf('mkdir %s',outDnameScale));
        end
        
        accLbpFnameAll = [accLbpPrefix num2str(iScale) '_' testStr '_all.mat'];        
                        
        tic;
        load([accLbpPrefix num2str(iScale) '_allInfo.mat']); % allInfo
        
        load(accLbpFnameAll); % these accumulations were created by pcAccumulatedLBPWell
        
        allCells = allCellsFov;
        clear allCellsFov;
                
        tt = toc;
        fprintf(sprintf('%s: done loading %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        
        tic;
        dataAll = allCells.accDeltaLbp;
        
        sourceAll = allInfo.source;
        metEffAll = metEff2Str(allInfo.metEff);
        
        clear allCells;
        n = length(dataAll);
        
        outDnameScaleAll = [outDnameScale filesep 'mixingLabelsMaps' filesep];
        if ~exist(outDnameScaleAll,'dir')
            unix(sprintf('mkdir %s',outDnameScaleAll));
        end
        
        if ~exist([outDnameScaleAll 'allDists.mat'],'file')
            [allDistMetEff] = getLbpDistributionsMixed(dataAll,metEffAll,generalStats);XXX
            [allDistSource] = getLbpDistributionsMixed(dataAll,sourceAll,generalStats);
            
            printLbpMaps(allMapsMetEff,outDnameScaleAll,'MetEff',unique(metEffAll));
            printLbpMaps(allMapsSource,outDnameScaleAll,'Source',unique(sourceAll));
            
            clear dataAll;
            save([outDnameScaleAll 'allMaps.mat'],'allMapsMetEff','allMapsSource');
            tt = toc;
            fprintf(sprintf('%s: done all: %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
            close all;
        end
        
        load([outDnameScaleAll 'allMaps.mat']); % allMapsMetEff,allMapsSource
        calcDistances(allMapsMetEff,[2,3],outDnameScaleAll,'D_MetEff_HL');
        calcDistances(allMapsSource,[1,3],outDnameScaleAll,'D_Source_CT');
        
    end
end
end

%%
function [] = calcDistances(allMaps,indsMaps,ourDir,prefix)
nRep = length(allMaps);

nConds = length(indsMaps);
D = nan(nRep*nConds);

for irep = 1 : nRep
    curMaps = allMaps{irep};
    for icond = 1 : nConds
        curMap = curMaps{indsMaps(icond)};
        i = (irep-1) * nConds + icond;
        
        for irep1 = 1 : nRep
            curMaps1 = allMaps{irep1};
            for icond1 = 1 : nConds
                curMap1 = curMaps1{indsMaps(icond1)};
                i1 = (irep1-1) * nConds + icond1;
                
                diff = abs(curMap - curMap1);
                D(i,i1) = sum(diff(:));
                D(i1,i) = D(i,i1);
            end
        end                
    end
end

fontsize = 24;
h = figure;
imagesc(D);
hold on;
title(prefix);
colormap('jet');
colorbar;caxis([0,2]);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',[]);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc([ourDir prefix '.eps']);

end

%%
function [dataNormGeneral,scoreGeneral] = ...
    getPCAs(data,generalStats)

dataNormGeneral = (data - repmat(generalStats.meanLbp,[1,size(data,2)]))./...
    repmat(generalStats.stdLbp,[1,size(data,2)]);

scoreGeneral = dataNormGeneral' * generalStats.coeffLbp;
end

%%
function [pcsMap] = getPcaMaps(score,generalStats)

nBins = 21;

pctl1Low = generalStats.PCLimitLow(1);
pctl1High = generalStats.PCLimitHigh(1);
pctl2Low = generalStats.PCLimitLow(2);
pctl2High = generalStats.PCLimitHigh(2);

xBins = pctl1Low : (pctl1High - pctl1Low)/(nBins-1) : pctl1High;
yBins = pctl2Low : (pctl2High - pctl2Low)/(nBins-1) : pctl2High;

xData = score(:,1);
yData = score(:,2);

[pcsMap] = getScatterQuantization(xData,yData,xBins,yBins);
end


%%
function [map] = getScatterQuantization(xData,yData,xBins,yBins)
N = length(xData);
sizeX = length(xBins) - 1;
sizeY = length(yBins) - 1;
map = zeros(sizeY,sizeX);
for x = 1 : sizeX
    for y = 1 : sizeY
        count = xData >= xBins(x) & xData < xBins(x+1) & yData >= yBins(y) & yData < yBins(y+1);
        count = sum(count) ./ double(N);
        map(sizeY - y + 1,x) = count;
    end
end
map = map ./ sum(map(:)); % normalize
end

%% allMaps2D{i} - maps2D; maps2D{i} is a map per label
function [allMaps2D] = getLbpDistributionsMixed(data,labels,generalStats)

n = length(data);

uniqueLabels = unique(labels);
nUniqueLabels = length(uniqueLabels);

nRep = 5; % first is actual data and rest with scrambeled labels
allMaps2D = cell(1,nRep);

nsLabels = nan(1,nUniqueLabels);
for ilabel = 1 : nUniqueLabels
    indsCellArray = strfind(labels,uniqueLabels{ilabel});
    inds = find(not(cellfun('isempty', indsCellArray)));
    nsLabels(ilabel) = length(inds); 
end

pLabels = nsLabels ./ sum(nsLabels);
pAccLabels = cumsum(pLabels);

for irep = 1 : nRep
    
    maps2D = cell(1,nUniqueLabels);
    
    dataUniqueLabels = cell(1,nUniqueLabels);
    for i = 1 : nUniqueLabels
        dataUniqueLabels{i} = [];
    end
    
    for i = 1 : n
        curData = data{i};
        
        if isempty(labels{i})
            continue;
        end
        
        if irep == 1
            curLabel = labels{i};
        else
            curLabel = getPropLabel(pAccLabels,uniqueLabels);
        end        
        
        for ilabel = 1 : nUniqueLabels
            if strcmp(curLabel,uniqueLabels{ilabel})
                dataUniqueLabels{ilabel} = [dataUniqueLabels{ilabel} curData];
            end
        end
    end
    for ilabel = 1 : nUniqueLabels
        if isempty(dataUniqueLabels{ilabel})
            continue;
        end
        [dataNormGeneral,scoreGeneral] = getPCAs(dataUniqueLabels{ilabel},generalStats);
        [curMap] = getPcaMaps(scoreGeneral,generalStats);
        maps2D{ilabel} = curMap;
    end
    allMaps2D{irep} = maps2D;
end
end

function [] = printLbpMaps(allMapsMetEff,outDnameScaleAll,repStr,labels)

nRep = length(allMapsMetEff);

for irep = 1 : nRep
    curMaps = allMapsMetEff{irep};
    n = length(curMaps);        
    
    for i = 1 : n       
        titleStr = [repStr '_' labels{i}];
        if ~isempty(curMaps{i})
            printCurMap(curMaps{i},titleStr,[outDnameScaleAll titleStr '_' num2str(irep) '.eps']);
        end
        close all;
    end
end
end


function [] = printCurMap(map,titleStr,outFname)
fontsize = 24;
h = figure;
imagesc(map);
hold on;
title(titleStr);
colormap('jet');
colorbar;caxis([0,0.02]);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',[]);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
xlabel('PC1','FontSize',fontsize); ylabel('PC2','FontSize',fontsize);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc(outFname);
end

%%
function metEffAll = metEff2Str(metEff)
n = length(metEff);
metEffAll = cell(1,n);
for i = 1 : n
    curMetEff = metEff{i};
    if isnan(curMetEff)
        metEffAll{i} = '';
    else
        if curMetEff
            metEffAll{i} = 'High';
        else
            metEffAll{i} = 'Low';
        end
    end
end
end

%%
function curLabel = getPropLabel(pAccLabels,uniqueLabels)
curLabel = nan;
r = rand();
for i = 1 : length(uniqueLabels) 
    if pAccLabels(i) >= r
        curLabel = uniqueLabels{i};
        break;
    end
end
assert(~isempty(curLabel));
end