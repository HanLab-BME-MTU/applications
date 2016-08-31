%% Accumulates LBP from all experiments for all time points (lbp vector for every cell at each time point):
% all data, cluster by cell type, cluster by Tumor/Cell Line/Melanocyte,
% cluster by high/low metastatic efficiency
% Assaf Zaritsky, Jan. 2016

function [] = pcAccumulateLBP_highMemory()

close all;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
accLbpFname = [analysisDirname 'metaAnalysis/accumulatedLBP.mat'];
metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);

nCellTypes = length(metaData.cellTypes.ids);

if ~exist(accLbpFname,'file')
    
    %% init
    accLbpFov.accLbpFov = cell(1,nScales);
    
    accLbpFov.cellTypes.strs = metaData.cellTypes.ids;
    accLbpFov.cellTypes.accLbpFov = cell(1,nCellTypes);
    
    accLbpFov.melanocytes.accLbpFov = cell(1,nScales);
    accLbpFov.cellLines.accLbpFov = cell(1,nScales);
    accLbpFov.tumorAll.accLbpFov = cell(1,nScales);
    
    accLbpFov.tumorHigh.accLbpFov = cell(1,nScales);
    accLbpFov.tumorLow.accLbpFov = cell(1,nScales);
    
    for iCellType = 1 : nCellTypes
        accLbpFov.cellTypes.accLbpFov{iCellType}.accLbpFov = cell(1,nScales);
    end
    
    for iScale = 1 : nScales
        accLbpFov.accLbpFov{iScale} = [];
        %         accLbpFov.cellTypes.accLbpFov{iScale} = [];
        
        accLbpFov.melanocytes.accLbpFov{iScale} = [];
        accLbpFov.cellLines.accLbpFov{iScale} = [];
        accLbpFov.tumorAll.accLbpFov{iScale} = [];
        
        accLbpFov.tumorHigh.accLbpFov{iScale} = [];
        accLbpFov.tumorLow.accLbpFov{iScale} = [];
        
        for iCellType = 1 : nCellTypes
            accLbpFov.cellTypes.accLbpFov{iCellType}.accLbpFov{iScale} = [];
        end
    end
    
    %% Accumulate 
    for itask = 1 : metaData.tasks.N
        curExp = metaData.tasks.exps(itask);
        curTask = metaData.tasks.tasks(itask);
        curFname = metaData.experiments.fnames{curExp};
        if curTask <= metaData.experiments.n1{curExp}
            curSource = metaData.experiments.source1{curExp};
            curCellType = metaData.experiments.cellType1{curExp};                        
        else
            curSource = metaData.experiments.source2{curExp};
            curCellType = metaData.experiments.cellType1{curExp};            
        end
        cellTypeInd = find(strcmpi(curCellType,metaData.cellTypes.ids));
        curMetEff = metaData.cellTypes.metastaticEfficiency(cellTypeInd);
        fprintf(sprintf('Accumulate LBP for %s_s%02d\n',curFname,curTask));
        
        lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%02d',curTask) filesep...
            'tracking' filesep 'lbpData.mat']; 
        
        if ~exist(lbpFname,'file')
            lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%d',curTask) filesep...
            'tracking' filesep 'lbpData.mat']; 
        end
        
        if ~exist(lbpFname,'file')
            error(['LBP file' lbpFname '(either 1 or 01) does not exist']);
        end
        
        load(lbpFname); % pyramidLBP
        
        %% Actual accumulation
        for iScale = 1 : nScales
            accLbpFov.accLbpFov{iScale} = [accLbpFov.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}]; %#ok<USENS> % lbpData.fov{icell}.lbp(curT,:)
            
            if strcmp(curSource,'Tumors')
                accLbpFov.tumorAll.accLbpFov{iScale} = [accLbpFov.tumorAll.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
            else if strcmp(curSource,'CellLines')
                    accLbpFov.cellLines.accLbpFov{iScale} = [accLbpFov.cellLines.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
                else if strcmp(curSource,'Melanocytes')
                        accLbpFov.melanocytes.accLbpFov{iScale} = [accLbpFov.melanocytes.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
                    end
                end
            end
            
            if ~isnan(curMetEff)
                if curMetEff
                    accLbpFov.tumorHigh.accLbpFov{iScale} = [accLbpFov.tumorHigh.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
                else
                    accLbpFov.tumorLow.accLbpFov{iScale} = [accLbpFov.tumorLow.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
                end
            end
            
            % Cell Type
            accLbpFov.cellTypes.accLbpFov{cellTypeInd}.accLbpFov{iScale} = [accLbpFov.cellTypes.accLbpFov{cellTypeInd}.accLbpFov{iScale},accumulatedFovPyramidLBP{iScale}];
            
            %     ncells = length(lbpData.fov);
            %     for icell = 1 : ncells
            %         accLbpFov = [accLbpFov, lbpData.fov{8}.lbp'];
            %     end
        end
    end
    
    save(accLbpFname,'accLbpFov','-v7.3');
else
    tic
    load(accLbpFname);
    toc
end

accLbpFov.meansLbpFov = mean(accLbpFov.accLbpFov,2); 
accLbpFov.stdsLbpFov = std(accLbpFov.accLbpFov')';
accLbpFov.accLbpFovNorm = (accLbpFov.accLbpFov - ...
    repmat(accLbpFov.meansLbpFov,[1,size(accLbpFov.accLbpFov,2)]))./...
    repmat(accLbpFov.stdsLbpFov,[1,size(accLbpFov.accLbpFov,2)]); % normalizing
[accLbpCoeff, accLbpScore, accLbpLatent] = princomp(accLbpFov.accLbpFovNorm');

accLbpFov.accLbpCoeff = accLbpCoeff;
accLbpFov.accLbpScore = accLbpScore;
accLbpFov.accLbpLatent = accLbpLatent;

%%

accLbpFov.cellTypes.accLbpFovPCA = cell(1,nCellTypes);

for i = 1 : nCellTypes
    if ~isempty(accLbpFov.cellTypes.accLbpFov{i})
        [accLbpFov.cellTypes.accLbpFovPCA{i}.pcaCoeff,...
            accLbpFov.cellTypes.accLbpFovPCA{i}.pcaScore,...
            accLbpFov.cellTypes.accLbpFovPCA{i}.pcaLatent] = normAndPCA(accLbpFov.cellTypes.accLbpFov{i});
    end
end

[accLbpFov.melanocytes.pcaCoeff, accLbpFov.melanocytes.pcaScore, accLbpFov.melanocytes.pcaLatent] = normAndPCA(accLbpFov.melanocytes.accLbpFov);
[accLbpFov.cellLines.pcaCoeff, accLbpFov.cellLines.pcaScore, accLbpFov.cellLines.pcaLatent] = normAndPCA(accLbpFov.cellLines.accLbpFov);
[accLbpFov.tumorAll.pcaCoeff, accLbpFov.tumorAll.pcaScore, accLbpFov.tumorAll.pcaLatent] = normAndPCA(accLbpFov.tumorAll.accLbpFov);

[accLbpFov.tumorHigh.pcaCoeff, accLbpFov.tumorHigh.pcaScore, accLbpFov.tumorHigh.pcaLatent] = normAndPCA(accLbpFov.tumorHigh.accLbpFov);
[accLbpFov.tumorLow.pcaCoeff,accLbpFov.tumorLow.pcaScore,accLbpFov.tumorLow.pcaLatent] = normAndPCA(accLbpFov.tumorLow.accLbpFov);

%%

save(accLbpFname,'accLbpFov','-v7.3');

fprintf('Accumulated variance:\n');
accLbpFov.accLbpLatent

%% Save 2D distribution of LBP space
percentile = 2;
nBins = 20;

pctl1Low = prctile(accLbpFov.accLbpScore(:,1),percentile);
pctl1High = prctile(accLbpFov.accLbpScore(:,1),100-percentile);
pctl2Low = prctile(accLbpFov.accLbpScore(:,2),percentile);
pctl2High = prctile(accLbpFov.accLbpScore(:,2),100-percentile);

xBins = pctl1Low : (pctl1High - pctl1Low)/(nBins-1) : pctl1High;
yBins = pctl2Low : (pctl2High - pctl2Low)/(nBins-1) : pctl2High;

xData = accLbpFov.accLbpScore(:,1);
yData = accLbpFov.accLbpScore(:,2);

[PCsMap] = getScatterQuantization(xData,yData,xBins,yBins);

accLbpFov.PCsMap = PCsMap;
accLbpFov.xBins = xBins;
accLbpFov.yBins = yBins;
accLbpFov.percentile = percentile;
accLbpFov.nBins = nBins;
save(accLbpFname,'accLbpFov','-v7.3');


end

%% ff

%% Normalize
function [pcaCoeff,pcaScore,pcaLatent] = normAndPCA(data)
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
        