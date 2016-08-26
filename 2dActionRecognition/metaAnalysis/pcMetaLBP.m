%% Differences (visualization) comparing between different cell types
% Assaf Zaritsky, Jan. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaLBP()

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
accLbpFname = [analysisDirname 'metaAnalysis/accumulatedLBP.mat'];
outDname = [analysisDirname 'metaAnalysis/LBP/FOV/'];

tic
load(accLbpFname); % accLbpFov
toc

if exist('tmp','var')
    accLbpFov = tmp;
end

generalMean = accLbpFov.meansLbpFov;
generalStd = accLbpFov.stdsLbpFov;
generalCoeff = accLbpFov.accLbpCoeff;

%% all
dataAll = accLbpFov.accLbpFov;
labelsAll = {'all'};
indsAll = [size(dataAll,2)];
outDirAll = [outDname 'all/'];
titleAll = 'all';

doMetaLBP(dataAll,labelsAll,indsAll,outDirAll,generalMean,generalStd,generalCoeff,titleAll);

fprintf('Finished all\n'); close all;

%% cell line / primary melanoma / melanocite
nCellLines = size(accLbpFov.cellLines.accLbpFov,2);
nMelanocytes = size(accLbpFov.melanocytes.accLbpFov,2);
nPrimaryMelanoma = size(accLbpFov.tumorAll.accLbpFov,2);
dataSource = [accLbpFov.cellLines.accLbpFov,accLbpFov.melanocytes.accLbpFov,accLbpFov.tumorAll.accLbpFov];
labelsSource = {'CellLine','Melano','Primary'};
indsSource = [nCellLines,nCellLines+nMelanocytes,nCellLines+nMelanocytes+nPrimaryMelanoma];
outDirSource = [outDname 'source/'];
titleSource = 'source';

doMetaLBP(dataSource,labelsSource,indsSource,outDirSource,generalMean,generalStd,generalCoeff,titleSource);

fprintf('Finished cell line / primary melanoma / melanocyte\n'); close all;

%% High vs. Low metastatic efficiency
nHighMet = size(accLbpFov.tumorHigh.accLbpFov,2);
nLowMet = size(accLbpFov.tumorLow.accLbpFov,2);
dataMetEff = [accLbpFov.tumorHigh.accLbpFov,accLbpFov.tumorLow.accLbpFov];
labelsMetEff = {'High','Low'};
indsMetEff = [nHighMet,nHighMet+nLowMet];
outDirMetEff = [outDname 'metEfficiency/'];
titleMetEff = 'metEff';

doMetaLBP(dataMetEff,labelsMetEff,indsMetEff,outDirMetEff,generalMean,generalStd,generalCoeff,titleMetEff);

fprintf('Finished metastatic efficiency \n'); close all;

%% Independent cell types
nCellTypes = length(accLbpFov.cellTypes.accLbpFov);

dataCellType = [];
labelsCellType = {};
indsCellType = [];
for i = 1 : nCellTypes    
    if isempty(accLbpFov.cellTypes.accLbpFov{i})
        continue;
    end
    curI = length(labelsCellType) + 1;
    dataCellType = [dataCellType, accLbpFov.cellTypes.accLbpFov{i}];
    labelsCellType{curI} = lower(accLbpFov.cellTypes.strs{i});
    if curI == 1
        curInd = 0;
    else
        curInd = indsCellType(end);
    end
    indsCellType = [indsCellType (curInd+size(accLbpFov.cellTypes.accLbpFov{i},2))];
end


outDirCellType = [outDname 'cellType/'];
titleCellType = 'cellType';

doMetaLBP(dataCellType,labelsCellType,indsCellType,outDirCellType,generalMean,generalStd,generalCoeff,titleCellType);

fprintf('Finished cell type \n'); close all;
end

%%

function [] = doMetaLBP(data,labels,inds,outDir,generalMean,generalStd,generalCoeff,titleStr)

close all;

if ~exist(outDir,'dir')
    unix(sprintf('mkdir %s',outDir));
end

fontsize = 6;

[dataNormGeneral,scoreGeneral,...
    dataNorm,coeff,score,latent] = getPCAs(data,generalMean,generalStd,generalCoeff);

plotLBP(dataNormGeneral,dataNorm,labels,inds,outDir,titleStr,fontsize);

[pca2DGeneral,pca2D] = plotPCA(scoreGeneral,coeff,score,latent,outDir,titleStr,fontsize);

for i = 1 : length(labels)
    if i == 1
        istart = 1;
    else
        istart = inds(i-1) + 1;
    end
    iend = inds(i);
    [pca2DGeneral,pca2D] = plotPCA(scoreGeneral(istart:iend,:),coeff,score(istart:iend,:),latent,outDir,[titleStr '_' labels{i}],fontsize);
end

% distsPCA(pca2DGeneral,pca2D,outDir,titleStr,fontsize);

end


%% 
function [dataNormGeneral,scoreGeneral,dataNorm,coeff,score,latent] = ...
    getPCAs(data,generalMean,generalStd,generalCoeff)

dataNormGeneral = (data - repmat(generalMean,[1,size(data,2)]))./...
    repmat(generalStd,[1,size(data,2)]);

scoreGeneral = dataNormGeneral' * generalCoeff;


[dataNorm,coeff,score,latent] = normAndPCA(data);
end

%%
function [] = plotLBP(dataNormGeneral,dataNorm,labels,inds,outDir,titleStr,fontsize)

h = figure;
imagesc(dataNormGeneral);
hold on;
colormap('jet');
colorbar;
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
export_fig([outDir 'lbpGeneral_' titleStr '.eps']);

h = figure;
imagesc(dataNorm);
hold on;
colormap('jet');
colorbar;
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
export_fig([outDir titleStr '_lbp.eps']);

end


function [pca2DGeneral,pca2D] = plotPCA(scoreGeneral,coeff,score,latent,outDir,titleStr,fontsize)
pca2D = plotPCAHelper(scoreGeneral,latent,outDir,[titleStr '_pca2DGeneral'],fontsize);
pca2DGeneral = plotPCAHelper(score,latent,outDir,[titleStr '_pca2D_'],fontsize);

plotCoeff(coeff,latent,outDir,[titleStr '_coeff'],fontsize);
end

function PCsMap = plotPCAHelper(score,latent,outDir,titleStr,fontsize)

percentile = 2;
nBins = 20;

pctl1Low = prctile(score(:,1),percentile);
pctl1High = prctile(score(:,1),100-percentile);
pctl2Low = prctile(score(:,2),percentile);
pctl2High = prctile(score(:,2),100-percentile);

xBins = pctl1Low : (pctl1High - pctl1Low)/(nBins-1) : pctl1High;
yBins = pctl2Low : (pctl2High - pctl2Low)/(nBins-1) : pctl2High;

xData = score(:,1);
yData = score(:,2);

[PCsMap] = getScatterQuantization(xData,yData,xBins,yBins);

h = figure;
imagesc(PCsMap);
hold on;
title(strrep([titleStr ' variance:' num2str(round(latent(1)*10)) '%, ' num2str(round(sum(latent(1:2))*10)) '%'],'_',' '),'FontSize',fontsize);
colormap('jet');
colorbar;
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
export_fig([outDir titleStr '.eps']);

end

function [] = plotCoeff(coeff,latent,outDir,titleStr,fontsize)

h = figure;
imagesc(coeff);
hold on;
title(strrep([titleStr ' variance:' num2str(round(latent(1)*10)) '%, ' num2str(round(sum(latent(1:2))*10)) '%'],'_',' '),'FontSize',fontsize);
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
export_fig([outDir titleStr '.eps']);
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
        