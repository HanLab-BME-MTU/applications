%% At the well level: visualizing LBP maps, collecting features
% Assaf Zaritsky, April. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaLBPWell()

always = true;

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBPWell/'];

outDnameFov = [lbpDirname filesep 'FOV/'];
outDnameBck = [lbpDirname filesep 'BCK/'];
outDnameFwd = [lbpDirname filesep 'FWD/'];
accLbpPrefix = [lbpDirname filesep 'accumulatedLBP_'];

if ~exist(outDnameFov,'dir')
    unix(sprintf('mkdir %s',outDnameFov));
end

if ~exist(outDnameBck,'dir')
    unix(sprintf('mkdir %s',outDnameBck));
end

if ~exist(outDnameFwd,'dir')
    unix(sprintf('mkdir %s',outDnameFwd));
end

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);
% nScales = 1; % Debug 201604
% scales = 1; % Debug 201604

outDirs = {outDnameFov,outDnameBck,outDnameFwd}; 
% outDirs = {outDnameFov}; % Debug 201604

testStrs = {'fov','bck','fwd'}; 
% testStrs = {'fov'}; % Debug 201604

% outDirs = {outDnameBck,outDnameFov};
% testStrs = {'bck','fov'};
assert(length(outDirs) == length(testStrs));
nTest = length(testStrs);

for iTest = 1 : nTest 
    outDname = outDirs{iTest};
    testStr = testStrs{iTest}; % field of view (all fov, background or forground)
    for iScale = 1 : nScales % resolution (1- maximal)
        outDnameScale = [outDname filesep num2str(iScale) filesep];
        
        if ~exist(outDnameScale,'dir')
            unix(sprintf('mkdir %s',outDnameScale));
        end
        
        accLbpFnameAll = [accLbpPrefix num2str(iScale) '_' testStr '_all.mat'];
        accLbpFnameSource = [accLbpPrefix num2str(iScale) '_' testStr '_source.mat'];
        accLbpFnameMetastatic = [accLbpPrefix num2str(iScale) '_' testStr '_metastatic.mat'];
        accLbpFnameCellType = [accLbpPrefix num2str(iScale) '_' testStr '_type.mat'];
        
        %% all
        allFeatsFname = [outDnameScale filesep 'allFeats.mat'];
        
        % take all cells to variable named "allCells" (for the specific field of view: fov, bck, fwd)
        tic;
        load([accLbpPrefix num2str(iScale) '_allInfo.mat']); % allInfo
        
        load(accLbpFnameAll); % these accumulations were created by pcAccumulateLBPWell
        if exist('allCellsBck','var')
            allCells = allCellsBck;
            clear allCellsBck;
        else if exist('allCellsFov','var')
                allCells = allCellsFov;
                clear allCellsFov;
            else if exist('allCellsFwd','var')
                    allCells = allCellsFwd;
                    clear allCellsFwd;
                end
            end
        end
        
        generalStats = allCells.pca;
        
        tt = toc;
        fprintf(sprintf('%s: done loading %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        
        if ~exist(allFeatsFname,'file') || always                                                
            tic;
            dataAll = allCells.accLbp;
            strAll = allInfo.strs;
            
            dataLocations = allCells.locations;%{curAll}.locationLbp{iLocation}
            
            dateAll = allInfo.date;
            cellTypeAll = allInfo.cellType;
            cellTypeIndAll = allInfo.cellTypeInd; % index in metaData.cellTypes.ids
            sourceAll = allInfo.source;
            metEffAll = allInfo.metEff;
            
            clear allCells;
            n = length(dataAll);
            
            [mapsAll,featsAll, mapsLocations,strLocations] = getLbpDistributionsAndFeatures(dataAll,dataLocations,generalStats);
            % all 2D distributions by experiment
            
            outDnameScaleAll = [outDnameScale filesep 'all' filesep];
            if ~exist(outDnameScaleAll,'dir')
                unix(sprintf('mkdir %s',outDnameScaleAll));
            end
            
%             unix(sprintf('rm %s*.eps',[outDnameScale filesep 'all']));
            
            printLbpMaps(mapsAll,strAll,mapsLocations,strLocations,outDnameScaleAll);
            clear dataAll;
            save(allFeatsFname,'mapsAll','featsAll','mapsLocations','strAll','dateAll','cellTypeAll','cellTypeIndAll','sourceAll','metEffAll','n');
            tt = toc;
            fprintf(sprintf('%s: done all: %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
            close all;
        end
        
%         %% cell line / primary melanoma / melanocyte
%         sourceFeatsFname = [outDnameScale filesep 'sourceFeats.mat'];                
%         
%         tic;
%         load(accLbpFnameSource);
%         
%         if exist('cellLinesBck','var')
%             cellLines = cellLinesBck;
%             melanocytes = melanocytesBck;
%             tumors = tumorsBck;
%             clear cellLinesBck melanocytesBck tumorsBck;
%         else if exist('cellLinesFov','var')
%                 cellLines = cellLinesFov;
%                 melanocytes = melanocytesFov;
%                 tumors = tumorsFov;
%                 clear cellLinesFov melanocytesFov tumorsFov;
%             else if exist('cellLinesFwd','var')
%                     cellLines = cellLinesFwd;
%                     melanocytes = melanocytesFwd;
%                     tumors = tumorsFwd;
%                     clear cellLinesFwd melanocytesFwd tumorsFwd;
%                 end
%             end
%         end
%         
%         nCellLines = length(cellLines.accLbp);
%         nMelanocytes = length(melanocytes.accLbp);
%         nPrimaryMelanoma = length(tumors.accLbp);
%         dataSource = [cellLines.accLbp,melanocytes.accLbp,tumors.accLbp];
%         clear cellLines melanocytes tumors;
%         labelsSource = {'CellLine','Melano','Primary'};
%         indsSource = [nCellLines,nCellLines+nMelanocytes,nCellLines+nMelanocytes+nPrimaryMelanoma];
%         outDirSource = [outDnameScale 'source/'];
%         titleSource = 'source';
%         
%         doMetaLBP(dataSource,labelsSource,indsSource,outDirSource,generalStats,titleSource);
%         
%         clear dataSource;
%         tt = toc;
%         fprintf(sprintf('%s: done source %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
%         close all;
%         
%         %% High vs. Low metastatic efficiency
%         tic;
%         load(accLbpFnameMetastatic);
%         
%         if exist('tumorHighBck','var')
%             tumorHigh = tumorHighBck;
%             tumorLow = tumorLowBck;            
%             clear tumorHighBck tumorLowBck;
%         else if exist('tumorHighFov','var')
%                 tumorHigh = tumorHighFov;
%                 tumorLow = tumorLowFov;                
%                 clear tumorHighFov tumorLowFov;
%             else if exist('tumorHighFwd','var')
%                     tumorHigh = tumorHighFwd;
%                     tumorLow = tumorLowFwd;
%                     clear tumorHighFwd tumorLowFwd;
%                 end
%             end
%         end
%                
%         nHighMet = size(tumorHigh.accLbp,2);
%         nLowMet = size(tumorLow.accLbp,2);
%         dataMetEff = [tumorHigh.accLbp,tumorLow.accLbp];
%         clear tumorHigh tumorLow;
%         labelsMetEff = {'High','Low'};
%         indsMetEff = [nHighMet,nHighMet+nLowMet];
%         outDirMetEff = [outDnameScale 'metEfficiency/'];
%         titleMetEff = 'metEff';
%         
%         doMetaLBP(dataMetEff,labelsMetEff,indsMetEff,outDirMetEff,generalStats,titleMetEff);
%         
%         clear dataMetEff;
%         tt = toc;
%         fprintf(sprintf('%s: done metastatic %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
%         close all;
%         
%         %% Independent cell types
%         tic;
%         load(accLbpFnameCellType);
%         
%         if exist('cellTypesBck','var')
%             cellTypes = cellTypesBck;
%             clear cellTypesBck;
%         else if exist('cellTypesFov','var')
%                 cellTypes = cellTypesFov;
%                 clear cellTypesFov;
%             else if exist('cellTypesFwd','var')
%                     cellTypes = cellTypesFwd;
%                     clear cellTypesFwd;
%                 end
%             end
%         end
%         
%         nCellTypes = length(cellTypes.accLbp);
%         
%         dataCellType = [];
%         labelsCellType = {};
%         indsCellType = [];
%         for i = 1 : nCellTypes
%             if isempty(cellTypes.accLbp{i})
%                 continue;
%             end
%             curI = length(labelsCellType) + 1;
%             dataCellType = [dataCellType, cellTypes.accLbp{i}];
%             labelsCellType{curI} = lower(cellTypes.strs{i});
%             if curI == 1
%                 curInd = 0;
%             else
%                 curInd = indsCellType(end);
%             end
%             indsCellType = [indsCellType (curInd+size(cellTypes.accLbp{i},2))];
%         end
%         
%         outDirCellType = [outDnameScale 'cellType/'];
%         titleCellType = 'cellType';
%         
%         doMetaLBP(dataCellType,labelsCellType,indsCellType,outDirCellType,generalStats,titleCellType);
%         clear dataCellType;
%         
%         tt = toc;
%         fprintf(sprintf('%s: done cell type %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
%         close all;
    end
end
end

%%


%%
function [dataNormGeneral,scoreGeneral,scoreLocationGeneral] = ...
    getPCAs(data,generalStats,dataLocations)

dataNormGeneral = (data - repmat(generalStats.meanLbp,[1,size(data,2)]))./...
    repmat(generalStats.stdLbp,[1,size(data,2)]);

scoreGeneral = dataNormGeneral' * generalStats.coeffLbp;

scoreLocationGeneral = cell(1,length(dataLocations.locationLbp));
for iLocation = 1 : length(scoreLocationGeneral)
    curData = dataLocations.locationLbp{iLocation};
    curDataNormGeneral = (curData - repmat(generalStats.meanLbp,[1,size(curData,2)]))./...
        repmat(generalStats.stdLbp,[1,size(curData,2)]);
    scoreLocationGeneral{iLocation} = curDataNormGeneral' * generalStats.coeffLbp;
end

end

%%
function [pcsMap, locationPcsMaps] = getPcaMaps(score,scoreLocationGeneral,generalStats)

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

% --------
nLocations = length(scoreLocationGeneral);
locationPcsMaps = cell(1,nLocations);

for iLocation = 1 : nLocations
    curScoreLocations = scoreLocationGeneral{iLocation};
    xLocationData = curScoreLocations(:,1);
    yLocationData = curScoreLocations(:,2);
    
    [curLocationPcsMap] = getScatterQuantization(xLocationData,yLocationData,xBins,yBins);
    locationPcsMaps{iLocation} = curLocationPcsMap;
end
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

function [maps2D,feats,maps2DLocations,strLocations] = getLbpDistributionsAndFeatures(data,dataLocations,generalStats)
% featW = 21/3;
% featH = 21/3;
n = length(data);

maps2D = cell(1,n);
maps2DLocations = cell(1,n);
strLocations = cell(1,n);

% feats = nan(featW*featH,n);
feats = [];
for i = 1 : n
    [dataNormGeneral,scoreGeneral,scoreLocationGeneral] = getPCAs(data{i},generalStats,dataLocations{i});
    [curMap, curLocationMaps] = getPcaMaps(scoreGeneral,scoreLocationGeneral,generalStats); % TODO: input - score by location, output - maps for locations 
    
    maps2D{i} = curMap;
    maps2DLocations{i} = curLocationMaps;  
    strLocations{i} = dataLocations{i}.locationStr;
    
    %     curMapLowRes = imresize(maps2D{i},[featH,featW]);
    %     feats(:,i) =  curMapLowRes(:);
    curMapLowRes = imresize(curMap,0.5);
    feats = [feats,curMapLowRes(:)];
    %     feats = [feats,curMap(:)];
end
end

function [] = printLbpMaps(maps,strs,mapsLocations,strLocations,outDir)
n = length(maps);

locationsOutDir = [outDir filesep 'Locations' filesep];
if ~exist(locationsOutDir,'dir')
    unix(sprintf('mkdir %s',locationsOutDir));
end

for i = 1 : n
    titleStr = strs{i};
    printCurMap(maps{i},titleStr,[outDir titleStr '.eps']);
    
    curMapsLocations = mapsLocations{i};
    curStrLocations = strLocations{i};        
    
    nLocations = length(curMapsLocations);
    for iLocation = 1 : nLocations
        curMapLocation = curMapsLocations{iLocation};
        printCurMap(curMapLocation,[titleStr '_' curStrLocations{iLocation}],[locationsOutDir titleStr '_s' curStrLocations{iLocation} '.eps']);
    end    
    close all;
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