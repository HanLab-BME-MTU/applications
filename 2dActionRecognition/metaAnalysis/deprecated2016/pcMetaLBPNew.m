%% Differences (visualization) comparing between different cell types
% Assaf Zaritsky, Jan. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaLBPNew()

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBP/'];

outDnameFov = [lbpDirname filesep 'FOV/'];
outDnameBck = [lbpDirname filesep 'BCK/'];
outDnameFwd = [lbpDirname filesep 'FWD/'];
accLbpPrefix = [lbpDirname filesep 'accumulatedLBP_'];

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);
% nScales = 1; % Debug 201604
% scales = 1; % Debug 201604

outDirs = {outDnameBck,outDnameFwd,outDnameFov}; 
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
        
        % take all cells to variable named "allCells" (for the specific field of view: fov, bck, fwd)
        tic;
        load(accLbpFnameAll);
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
        
        limitPercentile = 2;
        
        generalStats.mean = allCells.meansLbp;
        generalStats.std = allCells.stdsLbp;
        generalStats.coeff = allCells.accLbpCoeff;
        generalStats.PCLimitLow = prctile(allCells.accLbpScore,limitPercentile,1);
        generalStats.PCLimitHigh = prctile(allCells.accLbpScore,100-limitPercentile,1);
        generalStats.limitPrecentile = limitPercentile;
        generalStats.low = prctile(allCells.accLbpNorm(:),limitPercentile);
        generalStats.high = prctile(allCells.accLbpNorm(:),100-limitPercentile);
        generalStats.lowCoeff = prctile(allCells.accLbpCoeff(:),limitPercentile);
        generalStats.highCoeff = prctile(allCells.accLbpCoeff(:),100-limitPercentile);
        generalStats.latent = allCells.accLbpLatent;
        
        tt = toc;
        fprintf(sprintf('%s: done loading %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        
        %% all
        tic;
        dataAll = allCells.accLbp;
        clear allCells;
        labelsAll = {'all'};
        indsAll = size(dataAll,2);
        outDirAll = [outDnameScale 'all/'];
        titleAll = 'all';
        
        doMetaLBP(dataAll,labelsAll,indsAll,outDirAll,generalStats,titleAll);
        clear dataAll;
        tt = toc;
        fprintf(sprintf('%s: done all: %d %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        close all;
        
        %% cell line / primary melanoma / melanocite
        tic;
        load(accLbpFnameSource);
        
        if exist('cellLinesBck','var')
            cellLines = cellLinesBck;
            melanocytes = melanocytesBck;
            tumors = tumorsBck;
            clear cellLinesBck melanocytesBck tumorsBck;
        else if exist('cellLinesFov','var')
                cellLines = cellLinesFov;
                melanocytes = melanocytesFov;
                tumors = tumorsFov;
                clear cellLinesFov melanocytesFov tumorsFov;
            else if exist('cellLinesFwd','var')
                    cellLines = cellLinesFwd;
                    melanocytes = melanocytesFwd;
                    tumors = tumorsFwd;
                    clear cellLinesFwd melanocytesFwd tumorsFwd;
                end
            end
        end
        
        nCellLines = size(cellLines.accLbp,2);
        nMelanocytes = size(melanocytes.accLbp,2);
        nPrimaryMelanoma = size(tumors.accLbp,2);
        dataSource = [cellLines.accLbp,melanocytes.accLbp,tumors.accLbp];
        clear cellLines melanocytes tumors;
        labelsSource = {'CellLine','Melano','Primary'};
        indsSource = [nCellLines,nCellLines+nMelanocytes,nCellLines+nMelanocytes+nPrimaryMelanoma];
        outDirSource = [outDnameScale 'source/'];
        titleSource = 'source';
        
        doMetaLBP(dataSource,labelsSource,indsSource,outDirSource,generalStats,titleSource);
        
        clear dataSource;
        tt = toc;
        fprintf(sprintf('%s: done source %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        close all;
        
        %% High vs. Low metastatic efficiency
        tic;
        load(accLbpFnameMetastatic);
        
        if exist('tumorHighBck','var')
            tumorHigh = tumorHighBck;
            tumorLow = tumorLowBck;            
            clear tumorHighBck tumorLowBck;
        else if exist('tumorHighFov','var')
                tumorHigh = tumorHighFov;
                tumorLow = tumorLowFov;                
                clear tumorHighFov tumorLowFov;
            else if exist('tumorHighFwd','var')
                    tumorHigh = tumorHighFwd;
                    tumorLow = tumorLowFwd;
                    clear tumorHighFwd tumorLowFwd;
                end
            end
        end
               
        nHighMet = size(tumorHigh.accLbp,2);
        nLowMet = size(tumorLow.accLbp,2);
        dataMetEff = [tumorHigh.accLbp,tumorLow.accLbp];
        clear tumorHigh tumorLow;
        labelsMetEff = {'High','Low'};
        indsMetEff = [nHighMet,nHighMet+nLowMet];
        outDirMetEff = [outDnameScale 'metEfficiency/'];
        titleMetEff = 'metEff';
        
        doMetaLBP(dataMetEff,labelsMetEff,indsMetEff,outDirMetEff,generalStats,titleMetEff);
        
        clear dataMetEff;
        tt = toc;
        fprintf(sprintf('%s: done metastatic %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        close all;
        
        %% Independent cell types
        tic;
        load(accLbpFnameCellType);
        
        if exist('cellTypesBck','var')
            cellTypes = cellTypesBck;
            clear cellTypesBck;
        else if exist('cellTypesFov','var')
                cellTypes = cellTypesFov;
                clear cellTypesFov;
            else if exist('cellTypesFwd','var')
                    cellTypes = cellTypesFwd;
                    clear cellTypesFwd;
                end
            end
        end
        
        nCellTypes = length(cellTypes.accLbp);
        
        dataCellType = [];
        labelsCellType = {};
        indsCellType = [];
        for i = 1 : nCellTypes
            if isempty(cellTypes.accLbp{i})
                continue;
            end
            curI = length(labelsCellType) + 1;
            dataCellType = [dataCellType, cellTypes.accLbp{i}];
            labelsCellType{curI} = lower(cellTypes.strs{i});
            if curI == 1
                curInd = 0;
            else
                curInd = indsCellType(end);
            end
            indsCellType = [indsCellType (curInd+size(cellTypes.accLbp{i},2))];
        end
        
        outDirCellType = [outDnameScale 'cellType/'];
        titleCellType = 'cellType';
        
        doMetaLBP(dataCellType,labelsCellType,indsCellType,outDirCellType,generalStats,titleCellType);
        clear dataCellType;
        
        tt = toc;
        fprintf(sprintf('%s: done cell type %d: %s (%d min.)\n',testStr,iScale,'all',round(tt/60)));
        close all;
    end
end
end

%%

function [] = doMetaLBP(data,labels,inds,outDir,generalStats,titleStr)

close all;

if ~exist(outDir,'dir')
    unix(sprintf('mkdir %s',outDir));
end

fontsize = 6;

[dataNormGeneral,scoreGeneral,...
    dataNorm,coeff,score,latent] = getPCAs(data,generalStats);

specificStats.PCLimitLow = prctile(score,generalStats.limitPrecentile,1);
specificStats.PCLimitHigh = prctile(score,100-generalStats.limitPrecentile,1);
specificStats.lowCoeff = prctile(coeff(:),generalStats.limitPrecentile);
specificStats.highCoeff = prctile(coeff(:),100-generalStats.limitPrecentile);
specificStats.latent = latent;
specificStats.coeff = coeff;

plotLBP(dataNormGeneral,dataNorm,labels,inds,generalStats,outDir,titleStr,fontsize);

[pca2DGeneral,pca2D] = plotPCA(scoreGeneral,score,generalStats,specificStats,outDir,titleStr,fontsize);
plotCoeff(generalStats.coeff,generalStats,outDir,[titleStr '_general_coeff'],fontsize);
plotCoeff(specificStats.coeff,specificStats,outDir,[titleStr '_specific_coeff'],fontsize);


for i = 1 : length(labels)
    if i == 1
        istart = 1;
    else
        istart = inds(i-1) + 1;
    end
    iend = inds(i);
    [pca2DGeneral,pca2D] = plotPCA(scoreGeneral(istart:iend,:),score(istart:iend,:),generalStats,specificStats,outDir,[titleStr '_' labels{i}],fontsize);    
end

% distsPCA(pca2DGeneral,pca2D,outDir,titleStr,fontsize);

end


%%
function [dataNormGeneral,scoreGeneral,dataNorm,coeff,score,latent] = ...
    getPCAs(data,generalStats)

dataNormGeneral = (data - repmat(generalStats.mean,[1,size(data,2)]))./...
    repmat(generalStats.std,[1,size(data,2)]);

scoreGeneral = dataNormGeneral' * generalStats.coeff;


[dataNorm,coeff,score,latent] = normAndPCA(data);
end

%%
function [] = plotLBP(dataNormGeneral,dataNorm,labels,inds,generalStats,outDir,titleStr,fontsize)

h = figure;
imagesc(dataNormGeneral);
hold on;
colormap('jet');
colorbar;
caxis([generalStats.low,generalStats.high]);
title(strrep(['LBP (general): ' titleStr],'_',' '),'FontSize',fontsize);
% caxis(xaxisVals);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',inds);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',labels);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
xlabel('Cells','FontSize',fontsize); ylabel('LBP','FontSize',fontsize);
% colorBarYTick = [0,0.02];
% hColorbar = colorbar;
% set(hColorbar,'YTick',colorBarYTick);
% set(hColorbar,'YTickLabel',colorBarYTick);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc([outDir 'lbpGeneral_' titleStr '.eps']);

h = figure;
imagesc(dataNorm);
hold on;
colormap('jet');
colorbar;
caxis([generalStats.low,generalStats.high]);
title(strrep(['LBP (general): ' titleStr],'_',' '),'FontSize',fontsize);
% caxis(xaxisVals);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',inds);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',labels);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
xlabel('Cells','FontSize',fontsize); ylabel('LBP','FontSize',fontsize);
% colorBarYTick = [0,0.02];
% hColorbar = colorbar;
% set(hColorbar,'YTick',colorBarYTick);
% set(hColorbar,'YTickLabel',colorBarYTick);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc([outDir titleStr '_lbp.eps']);

end


function [pca2DGeneral,pca2D] = plotPCA(scoreGeneral,score,generalStats,specificStats,outDir,titleStr,fontsize)
pca2D = plotPCAHelper(scoreGeneral,generalStats,outDir,[titleStr '_pca2DGeneral'],fontsize);
pca2DGeneral = plotPCAHelper(score,specificStats,outDir,[titleStr '_pca2D_'],fontsize);
end

function PCsMap = plotPCAHelper(score,generalStats,outDir,titleStr,fontsize)

% percentile = 2;

nBins = 20;

pctl1Low = generalStats.PCLimitLow(1); %prctile(score(:,1),percentile);
pctl1High = generalStats.PCLimitHigh(1); % pctl1High = prctile(score(:,1),100-percentile);
pctl2Low = generalStats.PCLimitLow(2); % pctl2Low = prctile(score(:,2),percentile);
pctl2High = generalStats.PCLimitHigh(2); % pctl2High = prctile(score(:,2),100-percentile);

xBins = pctl1Low : (pctl1High - pctl1Low)/(nBins-1) : pctl1High;
yBins = pctl2Low : (pctl2High - pctl2Low)/(nBins-1) : pctl2High;

xData = score(:,1);
yData = score(:,2);

[PCsMap] = getScatterQuantization(xData,yData,xBins,yBins);

h = figure;
imagesc(PCsMap);
hold on;
title(strrep([titleStr ' variance:' num2str(round(generalStats.latent(1)*10)) '%, ' num2str(round(sum(generalStats.latent(1:2))*10)) '%'],'_',' '),'FontSize',fontsize);
colormap('jet');
colorbar;caxis([0,0.02]);
% caxis(xaxisVals);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',[]);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
xlabel('PC1','FontSize',fontsize); ylabel('PC2','FontSize',fontsize);
% colorBarYTick = [0,0.02];
% hColorbar = colorbar;
% set(hColorbar,'YTick',colorBarYTick);
% set(hColorbar,'YTickLabel',colorBarYTick);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc([outDir titleStr '.eps']);

end

function [] = plotCoeff(coeff,generalStats,outDir,titleStr,fontsize)

h = figure;
imagesc(coeff);
hold on;
caxis([generalStats.lowCoeff,generalStats.highCoeff]);
title(strrep([titleStr ' variance:' num2str(round(generalStats.latent(1)*10)) '%, ' num2str(round(sum(generalStats.latent(1:2))*10)) '%'],'_',' '),'FontSize',fontsize);
colormap('jet');
colorbar;
% caxis(xaxisVals);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',[]);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(haxes,'FontSize',fontsize);
xlabel('PCs','FontSize',fontsize); ylabel('Coefficients','FontSize',fontsize);
% colorBarYTick = [0,0.02];
% hColorbar = colorbar;
% set(hColorbar,'YTick',colorBarYTick);
% set(hColorbar,'YTickLabel',colorBarYTick);
set(h,'Color','w','PaperPositionMode','auto');
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(haxes,'XLabel'),'FontSize',fontsize); set(get(haxes,'YLabel'),'FontSize',fontsize);
export_fig_biohpc([outDir titleStr '.eps']);
end

%%

%% Normalize
function [dataNorm,pcaCoeff,pcaScore,pcaLatent] = normAndPCA(data)
meanVal = mean(data,2);
stdVal = std(data')';
dataNorm = (data - repmat(meanVal,[1,size(data,2)]))./...
    repmat(stdVal,[1,size(data,2)]);
[pcaCoeff,pcaScore,pcaLatent] = princomp(dataNorm');
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


% nCellTypes = length(metaData.cellTypes.ids);
%
% if ~exist(accLbpFname,'file')
%
%     %% init
%     accLbpFov.accLbpFov = [];
%
%     accLbpFov.cellTypes.strs = metaData.cellTypes.ids;
%     accLbpFov.cellTypes.accLbpFov = cell(1,nCellTypes);
%
%     accLbpFov.melanocytes.accLbpFov = [];
%     accLbpFov.cellLines.accLbpFov = [];
%     accLbpFov.tumorAll.accLbpFov = [];
%
%     accLbpFov.tumorHigh.accLbpFov = [];
%     accLbpFov.tumorLow.accLbpFov = [];
%

%
%     %% Accumulate
%     for itask = 1 : metaData.tasks.N
%         curExp = metaData.tasks.exps(itask);
%         curTask = metaData.tasks.tasks(itask);
%         curFname = metaData.experiments.fnames{curExp};
%         if curTask <= metaData.experiments.n1{curExp}
%             curSource = metaData.experiments.source1{curExp};
%             curCellType = metaData.experiments.cellType1{curExp};
%         else
%             curSource = metaData.experiments.source2{curExp};
%             curCellType = metaData.experiments.cellType1{curExp};
%         end
%         cellTypeInd = find(strcmpi(curCellType,metaData.cellTypes.ids));
%         curMetEff = metaData.cellTypes.metastaticEfficiency(cellTypeInd);
%         fprintf(sprintf('Accumulate LBP for %s_s%02d\n',curFname,curTask));
%
%         lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
%             curFname '_s' sprintf('%02d',curTask) filesep...
%             'tracking' filesep 'lbpData.mat']; % lbpData.fov
%
%         if ~exist(lbpFname,'file')
%             lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
%             curFname '_s' sprintf('%d',curTask) filesep...
%             'tracking' filesep 'lbpData.mat']; % lbpData.fov
%         end
%
%         if ~exist(lbpFname,'file')
%             error(['LBP file' lbpFname '(either 1 or 01) does not exist']);
%         end
%
%         load(lbpFname);
%
%         %% Actual accumulation
%         accLbpFov.accLbpFov = [accLbpFov.accLbpFov,accumulatedFovLbp]; % lbpData.fov{icell}.lbp(curT,:)
%
%         if strcmp(curSource,'Tumors')
%             accLbpFov.tumorAll.accLbpFov = [accLbpFov.tumorAll.accLbpFov,accumulatedFovLbp];
%         else if strcmp(curSource,'CellLines')
%                 accLbpFov.cellLines.accLbpFov = [accLbpFov.cellLines.accLbpFov,accumulatedFovLbp];
%             else if strcmp(curSource,'Melanocytes')
%                     accLbpFov.melanocytes.accLbpFov = [accLbpFov.melanocytes.accLbpFov,accumulatedFovLbp];
%                 end
%             end
%         end
%
%         if ~isnan(curMetEff)
%             if curMetEff
%                 accLbpFov.tumorHigh.accLbpFov = [accLbpFov.tumorHigh.accLbpFov,accumulatedFovLbp];
%             else
%                 accLbpFov.tumorLow.accLbpFov = [accLbpFov.tumorLow.accLbpFov,accumulatedFovLbp];
%             end
%         end
%
%         % Cell Type
%         accLbpFov.cellTypes.accLbpFov{cellTypeInd} = [accLbpFov.cellTypes.accLbpFov{cellTypeInd},accumulatedFovLbp];
%
%         %     ncells = length(lbpData.fov);
%         %     for icell = 1 : ncells
%         %         accLbpFov = [accLbpFov, lbpData.fov{8}.lbp'];
%         %     end
%     end
%
%     save(accLbpFname,'accLbpFov');
% else
%     load(accLbpFname);
% end
%
% accLbpFov.meansLbpFov = mean(accLbpFov.accLbpFov,2);
% accLbpFov.stdsLbpFov = std(accLbpFov.accLbpFov')';
% accLbpFov.accLbpFovNorm = (accLbpFov.accLbpFov - ...
%     repmat(accLbpFov.meansLbpFov,[1,size(accLbpFov.accLbpFov,2)]))./...
%     repmat(accLbpFov.stdsLbpFov,[1,size(accLbpFov.accLbpFov,2)]); % normalizing
% [accLbpCoeff, accLbpScore, accLbpLatent] = princomp(accLbpFov.accLbpFovNorm');
%
% accLbpFov.accLbpCoeff = accLbpCoeff;
% accLbpFov.accLbpScore = accLbpScore;
% accLbpFov.accLbpLatent = accLbpLatent;
%
% %%
%
% accLbpFov.cellTypes.accLbpFovPCA = cell(1,nCellTypes);
%
% for i = 1 : nCellTypes
%     if ~isempty(accLbpFov.cellTypes.accLbpFov{i})
%         [accLbpFov.cellTypes.accLbpFovPCA{i}.pcaCoeff,...
%             accLbpFov.cellTypes.accLbpFovPCA{i}.pcaScore,...
%             accLbpFov.cellTypes.accLbpFovPCA{i}.pcaLatent] = normAndPCA(accLbpFov.cellTypes.accLbpFov{i});
%     end
% end
%
% [accLbpFov.melanocytes.pcaCoeff, accLbpFov.melanocytes.pcaScore, accLbpFov.melanocytes.pcaLatent] = normAndPCA(accLbpFov.melanocytes.accLbpFov);
% [accLbpFov.cellLines.pcaCoeff, accLbpFov.cellLines.pcaScore, accLbpFov.cellLines.pcaLatent] = normAndPCA(accLbpFov.cellLines.accLbpFov);
% [accLbpFov.tumorAll.pcaCoeff, accLbpFov.tumorAll.pcaScore, accLbpFov.tumorAll.pcaLatent] = normAndPCA(accLbpFov.tumorAll.accLbpFov);
%
% [accLbpFov.tumorHigh.pcaCoeff, accLbpFov.tumorHigh.pcaScore, accLbpFov.tumorHigh.pcaLatent] = normAndPCA(accLbpFov.tumorHigh.accLbpFov);
% [accLbpFov.tumorLow.pcaCoeff,accLbpFov.tumorLow.pcaScore,accLbpFov.tumorLow.pcaLatent] = normAndPCA(accLbpFov.tumorLow.accLbpFov);
%
% %%
%
% save(accLbpFname,'accLbpFov');
%
% fprintf('Accumulated variance:\n');
% accLbpFov.accLbpLatent
%
% %% Save 2D distribution of LBP space
% percentile = 2;
% nBins = 20;
%
% pctl1Low = prctile(accLbpFov.accLbpScore(:,1),percentile);
% pctl1High = prctile(accLbpFov.accLbpScore(:,1),100-percentile);
% pctl2Low = prctile(accLbpFov.accLbpScore(:,2),percentile);
% pctl2High = prctile(accLbpFov.accLbpScore(:,2),100-percentile);
%
% xBins = pctl1Low : (pctl1High - pctl1Low)/(nBins-1) : pctl1High;
% yBins = pctl2Low : (pctl2High - pctl2Low)/(nBins-1) : pctl2High;
%
% xData = accLbpFov.accLbpScore(:,1);
% yData = accLbpFov.accLbpScore(:,2);
%
% [PCsMap] = getScatterQuantization(xData,yData,xBins,yBins);
%
% accLbpFov.PCsMap = PCsMap;
% accLbpFov.xBins = xBins;
% accLbpFov.yBins = yBins;
% accLbpFov.percentile = percentile;
% accLbpFov.nBins = nBins;
% save(accLbpFname,'accLbpFov');
%
%
% end
%
% %% ff
%
