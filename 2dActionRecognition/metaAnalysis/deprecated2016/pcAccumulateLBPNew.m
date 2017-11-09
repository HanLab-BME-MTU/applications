%% Accumulates LBP from all tasks for all time points (lbp vector for every cell at each time point - no use of time information):
% all data, cluster by cell type, cluster by Tumor/Cell Line/Melanocyte,
% cluster by high/low metastatic efficiency
% Assaf Zaritsky, Jan. 2016

function [] = pcAccumulateLBPNew()

close all;

always = true;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBP'];
accLbpPrefix = [lbpDirname filesep 'accumulatedLBP_'];
metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);%metaData

for iScale = 1 : nScales

    accLbpFnameAllFov = [accLbpPrefix num2str(iScale) '_fov_all.mat'];
    accLbpFnameSourceFov = [accLbpPrefix num2str(iScale) '_fov_source.mat'];
    accLbpFnameMetastaticFov = [accLbpPrefix num2str(iScale) '_fov_metastatic.mat'];
    accLbpFnameCellTypeFov = [accLbpPrefix num2str(iScale) '_fov_type.mat'];
    
    accLbpFnameAllBck = [accLbpPrefix num2str(iScale) '_bck_all.mat'];
    accLbpFnameSourceBck = [accLbpPrefix num2str(iScale) '_bck_source.mat'];
    accLbpFnameMetastaticBck = [accLbpPrefix num2str(iScale) '_bck_metastatic.mat'];
    accLbpFnameCellTypeBck = [accLbpPrefix num2str(iScale) '_bck_type.mat'];
    
    accLbpFnameAllFwd = [accLbpPrefix num2str(iScale) '_fwd_all.mat'];
    accLbpFnameSourceFwd = [accLbpPrefix num2str(iScale) '_fwd_source.mat'];
    accLbpFnameMetastaticFwd = [accLbpPrefix num2str(iScale) '_fwd_metastatic.mat'];
    accLbpFnameCellTypeFwd = [accLbpPrefix num2str(iScale) '_fwd_type.mat'];
    
    %% Accumulate 
    if ~exist(accLbpFnameAllFwd,'file') || always
        tic;
        out = accumulateLBP(metaData,analysisDirname,iScale,'all');
        allCellsFov = out.allCells.fov;
        allCellsBck = out.allCells.bck;
        allCellsFwd = out.allCells.fwd;
        % do FWD
        clear out;
        save(accLbpFnameAllFov,'allCellsFov','-v7.3');
        save(accLbpFnameAllBck,'allCellsBck','-v7.3');
        save(accLbpFnameAllFwd,'allCellsFwd','-v7.3');
        clear allCellsFov allCellsBck allCellsFwd;
        tt = toc;    
        fprintf(sprintf('done %d %s (%d min.)\n',iScale,'all',round(tt/60)));
    end        
    
    if ~exist(accLbpFnameCellTypeFwd,'file') || always
        tic;
        out = accumulateLBP(metaData,analysisDirname,iScale,'type');
        cellTypesFov = out.cellTypes.fov;
        cellTypesBck = out.cellTypes.bck;
        cellTypesFwd = out.cellTypes.fwd;
        % do FWD
        clear out;
        save(accLbpFnameCellTypeFov,'cellTypesFov','-v7.3');
        clear cellTypesFov
        save(accLbpFnameCellTypeBck,'cellTypesBck','-v7.3');
        clear cellTypesBck;
        save(accLbpFnameCellTypeFwd,'cellTypesFwd','-v7.3');
        clear cellTypesFwd;
        tt = toc;
        fprintf(sprintf('done %d %s\n',iScale,'type',round(tt/60)));
    end
    
    if ~exist(accLbpFnameSourceFwd,'file') || always
        tic;
        out = accumulateLBP(metaData,analysisDirname,iScale,'source');
        melanocytesFov = out.melanocytes.fov;
        cellLinesFov = out.cellLines.fov;
        tumorsFov = out.tumors.fov;
        melanocytesBck = out.melanocytes.bck;
        cellLinesBck = out.cellLines.bck;
        tumorsBck = out.tumors.bck;
        melanocytesFwd = out.melanocytes.fwd;
        cellLinesFwd = out.cellLines.fwd;
        tumorsFwd = out.tumors.fwd;
        % do FWD
        clear out;
        save(accLbpFnameSourceFov,'melanocytesFov','cellLinesFov','tumorsFov','-v7.3');
        clear melanocytesFov cellLinesFov tumorsFov;
        save(accLbpFnameSourceBck,'melanocytesBck','cellLinesBck','tumorsBck','-v7.3');
        clear melanocytesBck cellLinesBck tumorsBck;
        save(accLbpFnameSourceFwd,'melanocytesFwd','cellLinesFwd','tumorsFwd','-v7.3');
        clear melanocytesFwd cellLinesFwd tumorsFwd;
        tt = toc;   
        fprintf(sprintf('done %d %s\n',iScale,'source',round(tt/60)));
    end
    
    if ~exist(accLbpFnameMetastaticFwd,'file') || always
        tic;
        out = accumulateLBP(metaData,analysisDirname,iScale,'metastatic');
        tumorHighFov = out.tumorHigh.fov;
        tumorLowFov = out.tumorLow.fov;
        tumorHighBck = out.tumorHigh.bck;
        tumorLowBck = out.tumorLow.bck;
        tumorHighFwd = out.tumorHigh.fwd;
        tumorLowFwd = out.tumorLow.fwd;
        % do FWD
        clear out;
        save(accLbpFnameMetastaticFov,'tumorHighFov','tumorLowFov','-v7.3');
        clear tumorHighFov tumorLowFov;
        save(accLbpFnameMetastaticBck,'tumorHighBck','tumorLowBck','-v7.3');
        clear tumorHighBck tumorLowBck;
        save(accLbpFnameMetastaticFwd,'tumorHighFwd','tumorLowFwd','-v7.3');
        clear tumorHighFwd tumorLowFwd;
        tt = toc;   
        fprintf(sprintf('done %d %s\n',iScale,'metastatic',round(tt/60)));
    end
% end
end
end

%%
function out = accumulateLBP(metaData,analysisDirname,iScale,strLabel)
%% init
nCellTypes = length(metaData.cellTypes.ids);
out.allCells.fov.accLbp = [];
out.allCells.bck.accLbp = [];
out.allCells.fwd.accLbp = [];

out.cellTypes.fov.strs = [];
out.cellTypes.bck.strs = [];
out.cellTypes.fwd.strs = [];

out.cellTypes.fov.accLbp = cell(1,nCellTypes);
out.cellTypes.bck.accLbp = cell(1,nCellTypes);
out.cellTypes.fwd.accLbp = cell(1,nCellTypes);

out.melanocytes.fov.accLbp = [];
out.cellLines.fov.accLbp = [];
out.tumors.fov.accLbp = [];
out.melanocytes.bck.accLbp = [];
out.cellLines.bck.accLbp = [];
out.tumors.bck.accLbp = [];
out.melanocytes.fwd.accLbp = [];
out.cellLines.fwd.accLbp = [];
out.tumors.fwd.accLbp = [];

out.tumorHigh.fov.accLbp = [];
out.tumorLow.fov.accLbp = [];
out.tumorHigh.bck.accLbp = [];
out.tumorLow.bck.accLbp = [];
out.tumorHigh.fwd.accLbp = [];
out.tumorLow.fwd.accLbp = [];

for iCellType = 1 : nCellTypes
    out.cellTypes.fov.accLbp{iCellType} = [];
    out.cellTypes.bck.accLbp{iCellType} = [];
    out.cellTypes.fwd.accLbp{iCellType} = [];
end

condInd = find(strcmp(strLabel,{'all','type','source','metastatic'}));

%% Accumulation
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
    
    load(lbpFname); % pyramidLBP (?), accumulatedFovPyramidLBP
        
    %% Actual accumulation        
    curLBPsFov = accumulatedFovPyramidLBP{iScale}; %#ok<USENS>
    curLBPsBck = accumulatedBckPyramidLBP{iScale}; %#ok<USENS>
    curLBPsFwd = accumulatedFwdPyramidLBP{iScale}; %#ok<USENS>
    
    if condInd == 1 % all
        out.allCells.fov.accLbp = [out.allCells.fov.accLbp,curLBPsFov];
        out.allCells.bck.accLbp = [out.allCells.bck.accLbp,curLBPsBck];
        out.allCells.fwd.accLbp = [out.allCells.fwd.accLbp,curLBPsFwd];
    end
    
    if condInd == 2 % type
        out.cellTypes.fov.accLbp{cellTypeInd} = [out.cellTypes.fov.accLbp{cellTypeInd},curLBPsFov];
        out.cellTypes.bck.accLbp{cellTypeInd} = [out.cellTypes.bck.accLbp{cellTypeInd},curLBPsBck];
        out.cellTypes.fwd.accLbp{cellTypeInd} = [out.cellTypes.fwd.accLbp{cellTypeInd},curLBPsBck];
    end
    
    if condInd == 3 % source
        if strcmp(curSource,'Tumors')
            out.tumors.fov.accLbp = [out.tumors.fov.accLbp,curLBPsFov];
            out.tumors.bck.accLbp = [out.tumors.bck.accLbp,curLBPsBck];
            out.tumors.fwd.accLbp = [out.tumors.fwd.accLbp,curLBPsFwd];
        else if strcmp(curSource,'CellLines')
                out.cellLines.fov.accLbp = [out.cellLines.fov.accLbp,curLBPsFov];
                out.cellLines.bck.accLbp = [out.cellLines.bck.accLbp,curLBPsBck];
                out.cellLines.fwd.accLbp = [out.cellLines.fwd.accLbp,curLBPsFwd];
            else if strcmp(curSource,'Melanocytes')
                    out.melanocytes.fov.accLbp = [out.melanocytes.fov.accLbp,curLBPsFov];
                    out.melanocytes.bck.accLbp = [out.melanocytes.bck.accLbp,curLBPsBck];
                    out.melanocytes.fwd.accLbp = [out.melanocytes.fwd.accLbp,curLBPsFwd];
                end
            end
        end
    end
    
    if condInd == 4 % metastatic
        if ~isnan(curMetEff)
            if curMetEff
                out.tumorHigh.fov.accLbp = [out.tumorHigh.fov.accLbp,curLBPsFov];
                out.tumorHigh.bck.accLbp = [out.tumorHigh.bck.accLbp,curLBPsBck];
                out.tumorHigh.fwd.accLbp = [out.tumorHigh.fwd.accLbp,curLBPsFwd];
            else
                out.tumorLow.fov.accLbp = [out.tumorLow.fov.accLbp,curLBPsFov];
                out.tumorLow.bck.accLbp = [out.tumorLow.bck.accLbp,curLBPsBck];
                out.tumorLow.fwd.accLbp = [out.tumorLow.fwd.accLbp,curLBPsBck];
            end
        end
    end
end

%% general stats - FOV
if condInd == 1 % all
    out.allCells.fov.meansLbp = mean(out.allCells.fov.accLbp,2);
    out.allCells.fov.stdsLbp = std(out.allCells.fov.accLbp')';
    out.allCells.fov.accLbpNorm = (out.allCells.fov.accLbp - ...
        repmat(out.allCells.fov.meansLbp,[1,size(out.allCells.fov.accLbp,2)]))./...
        repmat(out.allCells.fov.stdsLbp,[1,size(out.allCells.fov.accLbp,2)]); % normalizing
    [accLbpCoeff, accLbpScore, accLbpLatent] = princomp(out.allCells.fov.accLbpNorm');
    
    out.allCells.fov.accLbpCoeff = accLbpCoeff;
    out.allCells.fov.accLbpScore = accLbpScore;
    out.allCells.fov.accLbpLatent = accLbpLatent;
end


if condInd == 2 % type
    out.cellTypes.fov.accLbpPCA = cell(1,nCellTypes);    
    for i = 1 : nCellTypes
        if ~isempty(out.cellTypes.fov.accLbp{i})
            [out.cellTypes.fov.accLbpPCA{i}.pcaCoeff,...
                out.cellTypes.fov.accLbpPCA{i}.pcaScore,...
                out.cellTypes.fov.accLbpPCA{i}.pcaLatent] = normAndPCA(out.cellTypes.fov.accLbp{i});
        end
    end
    out.cellTypes.fov.strs = metaData.cellTypes.ids;
    out.cellTypes.bck.strs = metaData.cellTypes.ids;
    out.cellTypes.fwd.strs = metaData.cellTypes.ids;
end


if condInd == 3 % source
    [out.melanocytes.fov.pcaCoeff, out.melanocytes.fov.pcaScore, out.melanocytes.fov.pcaLatent] = normAndPCA(out.melanocytes.fov.accLbp);
    [out.cellLines.fov.pcaCoeff, out.cellLines.fov.pcaScore, out.cellLines.fov.pcaLatent] = normAndPCA(out.cellLines.fov.accLbp);
    [out.tumors.fov.pcaCoeff, out.tumors.fov.pcaScore, out.tumors.fov.pcaLatent] = normAndPCA(out.tumors.fov.accLbp);
end

if condInd == 4 % metastatic
    [out.tumorHigh.fov.pcaCoeff, out.tumorHigh.fov.pcaScore, out.tumorHigh.fov.pcaLatent] = normAndPCA(out.tumorHigh.fov.accLbp);
    [out.tumorLow.fov.pcaCoeff,out.tumorLow.fov.pcaScore,out.tumorLow.fov.pcaLatent] = normAndPCA(out.tumorLow.fov.accLbp);    
end

%% general stats - BCK
if condInd == 1 % all
    out.allCells.bck.meansLbp = mean(out.allCells.bck.accLbp,2);
    out.allCells.bck.stdsLbp = std(out.allCells.bck.accLbp')';
    out.allCells.bck.accLbpNorm = (out.allCells.bck.accLbp - ...
        repmat(out.allCells.bck.meansLbp,[1,size(out.allCells.bck.accLbp,2)]))./...
        repmat(out.allCells.bck.stdsLbp,[1,size(out.allCells.bck.accLbp,2)]); % normalizing
    [accLbpCoeff, accLbpScore, accLbpLatent] = princomp(out.allCells.bck.accLbpNorm');
    
    out.allCells.bck.accLbpCoeff = accLbpCoeff;
    out.allCells.bck.accLbpScore = accLbpScore;
    out.allCells.bck.accLbpLatent = accLbpLatent;
end


if condInd == 2 % type
    out.cellTypes.bck.accLbpPCA = cell(1,nCellTypes);    
    for i = 1 : nCellTypes
        if ~isempty(out.cellTypes.bck.accLbp{i})
            [out.cellTypes.bck.accLbpPCA{i}.pcaCoeff,...
                out.cellTypes.bck.accLbpPCA{i}.pcaScore,...
                out.cellTypes.bck.accLbpPCA{i}.pcaLatent] = normAndPCA(out.cellTypes.bck.accLbp{i});
        end
    end
    out.cellTypes.strs = metaData.cellTypes.ids;
end


if condInd == 3 % source
    [out.melanocytes.bck.pcaCoeff, out.melanocytes.bck.pcaScore, out.melanocytes.bck.pcaLatent] = normAndPCA(out.melanocytes.bck.accLbp);
    [out.cellLines.bck.pcaCoeff, out.cellLines.bck.pcaScore, out.cellLines.bck.pcaLatent] = normAndPCA(out.cellLines.bck.accLbp);
    [out.tumors.bck.pcaCoeff, out.tumors.bck.pcaScore, out.tumors.bck.pcaLatent] = normAndPCA(out.tumors.bck.accLbp);
end

if condInd == 4 % metastatic
    [out.tumorHigh.bck.pcaCoeff, out.tumorHigh.bck.pcaScore, out.tumorHigh.bck.pcaLatent] = normAndPCA(out.tumorHigh.bck.accLbp);
    [out.tumorLow.bck.pcaCoeff,out.tumorLow.bck.pcaScore,out.tumorLow.bck.pcaLatent] = normAndPCA(out.tumorLow.bck.accLbp);    
end

% general stats - Fwd
if condInd == 1 % all
    out.allCells.fwd.meansLbp = mean(out.allCells.fwd.accLbp,2);
    out.allCells.fwd.stdsLbp = std(out.allCells.fwd.accLbp')';
    out.allCells.fwd.accLbpNorm = (out.allCells.fwd.accLbp - ...
        repmat(out.allCells.fwd.meansLbp,[1,size(out.allCells.fwd.accLbp,2)]))./...
        repmat(out.allCells.fwd.stdsLbp,[1,size(out.allCells.fwd.accLbp,2)]); % normalizing
    [accLbpCoeff, accLbpScore, accLbpLatent] = princomp(out.allCells.fwd.accLbpNorm');
    
    out.allCells.fwd.accLbpCoeff = accLbpCoeff;
    out.allCells.fwd.accLbpScore = accLbpScore;
    out.allCells.fwd.accLbpLatent = accLbpLatent;
end


if condInd == 2 % type
    out.cellTypes.fwd.accLbpPCA = cell(1,nCellTypes);    
    for i = 1 : nCellTypes
        if ~isempty(out.cellTypes.fwd.accLbp{i})
            [out.cellTypes.fwd.accLbpPCA{i}.pcaCoeff,...
                out.cellTypes.fwd.accLbpPCA{i}.pcaScore,...
                out.cellTypes.fwd.accLbpPCA{i}.pcaLatent] = normAndPCA(out.cellTypes.fwd.accLbp{i});
        end
    end
    out.cellTypes.strs = metaData.cellTypes.ids;
end


if condInd == 3 % source
    [out.melanocytes.fwd.pcaCoeff, out.melanocytes.fwd.pcaScore, out.melanocytes.fwd.pcaLatent] = normAndPCA(out.melanocytes.fwd.accLbp);
    [out.cellLines.fwd.pcaCoeff, out.cellLines.fwd.pcaScore, out.cellLines.fwd.pcaLatent] = normAndPCA(out.cellLines.fwd.accLbp);
    [out.tumors.fwd.pcaCoeff, out.tumors.fwd.pcaScore, out.tumors.fwd.pcaLatent] = normAndPCA(out.tumors.fwd.accLbp);
end

if condInd == 4 % metastatic
    [out.tumorHigh.fwd.pcaCoeff, out.tumorHigh.fwd.pcaScore, out.tumorHigh.fwd.pcaLatent] = normAndPCA(out.tumorHigh.fwd.accLbp);
    [out.tumorLow.fwd.pcaCoeff,out.tumorLow.fwd.pcaScore,out.tumorLow.fwd.pcaLatent] = normAndPCA(out.tumorLow.fwd.accLbp);    
end


end

%% Normalize
function [pcaCoeff,pcaScore,pcaLatent] = normAndPCA(data)
meanVal = mean(data,2);
stdVal = std(data')';
dataNorm = (data - repmat(meanVal,[1,size(data,2)]))./...
    repmat(stdVal,[1,size(data,2)]); 
[pcaCoeff,pcaScore,pcaLatent] = princomp(dataNorm');


end
        