% Collect single cell LBP
function [] = singleCellLBP(dbFname,analysisDname)

close all;

%% File names
if nargin == 0
    dbFname = 'singleCellLabelDB.mat';
end

if nargin < 2
    analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
end

dbFname = [analysisDname 'SingleCellAnalysis/' dbFname];
outDname = [analysisDname 'SingleCellAnalysis/LBP/'];

load(dbFname); % cellInfoDB,TYXDB

accLbpFname = [analysisDname 'metaAnalysis/LBP/accumulatedLBP_1_all.mat'];
load(accLbpFname); % accLbpFov

ncells = length(cellInfoDB);
for icell = 1 : ncells    
    recordLbp(cellInfoDB{icell},analysisDname,allCells.fov.LbpFov,outDname);
end

end

function [] = recordLbp(cellInfo,analysisDname,accLbpFov,outDname)
close all;
lbpFname = [analysisDname filesep 'Data' filesep cellInfo.source filesep cellInfo.fname(1:end-4) filesep ...
    cellInfo.fname filesep 'tracking/lbpData.mat'];
load(lbpFname); % lbpData.fov{1}.lbp time x 10
cellLbpFov = lbpData.fov{cellInfo.serialNum}.lbp;

% cellLbpFovNorm = cellLbpFov - repmat(accLbpFov.meansLbpFov,[size(cellLbpFov,1),1])./repmat(accLbpFov.stdsLbpFov,[size(cellLbpFov,1),1]); % normalizing

cellLbpFovNorm = (cellLbpFov' - repmat(accLbpFov.meansLbpFov,[1,size(cellLbpFov,1)]))./...
    repmat(accLbpFov.stdsLbpFov,[1,size(cellLbpFov,1)]);

cellScoreGeneral = cellLbpFovNorm' * accLbpFov.accLbpCoeff;
pcs1 = cellScoreGeneral(:,1);
pcs2 = cellScoreGeneral(:,2);

figure;
subplot(3,1,1);
imagesc(accLbpFov.PCsMap);
hold on;
[binPCs1,binPCs2] = get2DPcBins(pcs1,pcs2,accLbpFov.xBins,accLbpFov.yBins);
plot(binPCs1,binPCs2,'k-','LineWidth',2);
hold off;

subplot(3,1,2);
imagesc(cellLbpFovNorm);
hold on;
hold off;

subplot(3,1,3);
imagesc(cellLbpFovNorm - repmat(cellLbpFovNorm(:,1),[1,size(cellLbpFovNorm,2)]));
hold on;
hold off;

outFname = [outDname filesep num2str(cellInfo.serialNum) '_' cellInfo.fname '_t' num2str(cellInfo.ts(1)) '_y' num2str(round(cellInfo.ys(1))) '_x' num2str(round(cellInfo.xs(1))) '.eps'];
export_fig(outFname);
end

function [binPCs1,binPCs2] = get2DPcBins(pcs1,pcs2,xBins,yBins)
n = length(pcs1);
binPCs1 = nan(1,n); binPCs2 = nan(1,n);
for i = 1 : length(pcs1)
    [binPCs1(i), binPCs2(i)] = get2DPcBin(pcs1(i),pcs2(i),xBins,yBins);
end
end

function [binPC1,binPC2] = get2DPcBin(pc1,pc2,xBins,yBins)
binPC1 = find(xBins>=pc1,1);
binPC2 = find(yBins>=pc2,1);

if isempty(binPC1)
    binPC1 = length(xBins);
end


if isempty(binPC2)
    binPC2 = length(yBins);
end

end

