%% Accumulates LBP for taks by experiment (well)
% (lbp vector for every cell at each time point - no use of time information):
% all data, cluster by cell type, cluster by Tumor/Cell Line/Melanocyte,
% cluster by high/low metastatic efficiency

% allCells include all the information that is needed for control-analysis
% / classifications / distnace maps

% Assaf Zaritsky, April. 2016

function [] = pcAccumulateLBPWell()

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

always = true;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBPWell'];
accLbpPrefix = [lbpDirname filesep 'accumulatedLBP_'];
metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);%metaData

for iScale = 1 : nScales
    
    % these are different from the pcAccumulateLBPNew, holding an array of
    % features, one per well
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
    
    pcaStatsGeneralPrefix = [lbpDirname '/../LBP/'];
    
    %% Accumulate
    if ~exist(accLbpFnameAllFwd,'file') || always
        tic;
        out = accumulateLBPWell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,'all');
        allCellsFov = out.allCells.fov;
        allCellsBck = out.allCells.bck;
        allCellsFwd = out.allCells.fwd;
        rmfield(out.allCells,'fov');
        rmfield(out.allCells,'bck');
        rmfield(out.allCells,'fwd');
        allInfo = out.allCells;
        % do FWD
        clear out;
        save(accLbpFnameAllFov,'allCellsFov','-v7.3');
        save(accLbpFnameAllBck,'allCellsBck','-v7.3');
        save(accLbpFnameAllFwd,'allCellsFwd','-v7.3');
        save([accLbpPrefix num2str(iScale) '_allInfo.mat'],'allInfo','-v7.3');
        clear allCellsFov allCellsBck allCellsFwd;
        tt = toc;
        fprintf(sprintf('done %d %s (%d min.)\n',iScale,'all',round(tt/60)));
    end
    
    if ~exist(accLbpFnameCellTypeFwd,'file') || always
        tic;
        out = accumulateLBPWell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,'type');
        cellTypesFov = cell(1,length(out.cellTypes));
        cellTypesBck = cell(1,length(out.cellTypes));
        cellTypesFwd = cell(1,length(out.cellTypes));
        for ict = 1 : length(out.cellTypes)
            cellTypesFov{ict} = out.cellTypes{ict}.fov; % remember that now accLbp has multiple entries, one for each experiment
            cellTypesBck{ict} = out.cellTypes{ict}.bck;
            cellTypesFwd{ict} = out.cellTypes{ict}.fwd;
        end
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
        out = accumulateLBPWell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,'source');
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
        out = accumulateLBPWell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,'metastatic');
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

%% Here we do the accumulation for every well!
function out = accumulateLBPWell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,strLabel)

%% PCA data is taken from the single cell, single time point data! (only general)
[out.allCells.fov.pca, out.allCells.bck.pca, out.allCells.fwd.pca] = ...
    getPcaGeneral(pcaStatsGeneralPrefix,iScale); 


%% init
nCellTypes = length(metaData.cellTypes.ids);
out.allCells.fov.accLbp = {};
out.allCells.bck.accLbp = {};
out.allCells.fwd.accLbp = {};
out.allCells.strs = {};
out.allCells.date = {};
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

out.cellTypes = cell(1,nCellTypes);
out.cellTypes = cell(1,nCellTypes);
out.cellTypes = cell(1,nCellTypes);

out.melanocytes.fov.accLbp = {};
out.cellLines.fov.accLbp = {};
out.tumors.fov.accLbp = {};
out.melanocytes.bck.accLbp = {};
out.cellLines.bck.accLbp = {};
out.tumors.bck.accLbp = {};
out.melanocytes.fwd.accLbp = {};
out.cellLines.fwd.accLbp = {};
out.tumors.fwd.accLbp = {};
out.melanocytes.strs = {};
out.cellLines.strs = {};
out.tumors.strs = {};

out.tumorHigh.fov.accLbp = {};
out.tumorLow.fov.accLbp = {};
out.tumorHigh.bck.accLbp = {};
out.tumorLow.bck.accLbp = {};
out.tumorHigh.fwd.accLbp = {};
out.tumorLow.fwd.accLbp = {};
out.tumorHigh.strs = {};
out.tumorLow.strs = {};

for iCellType = 1 : nCellTypes
    out.cellTypes{iCellType}.fov.accLbp = {};
    out.cellTypes{iCellType}.bck.accLbp = {};
    out.cellTypes{iCellType}.fwd.accLbp = {};
    out.cellTypes{iCellType}.strs = {};
end

condInd = find(strcmp(strLabel,{'all','type','source','metastatic'}));

%% Accumulation
for iexp = 1 : metaData.experiments.N
    curFname = metaData.experiments.fnames{iexp};
    curDate = curFname(1:6);
    for in = 1 : 2
        %% count the next open spot for each condition
        curAll = length(out.allCells.fov.accLbp) + 1;
        
        curTypes = nan(1,nCellTypes);
        for iCellType = 1 : nCellTypes
            curTypes(iCellType) = length(out.cellTypes{iCellType}.fov.accLbp) + 1;
        end
        
        curMelanocytes = length(out.melanocytes.fov.accLbp) + 1;
        curCellLines = length(out.cellLines.fov.accLbp) + 1;
        curTumors = length(out.tumors.fov.accLbp) + 1;
        
        curTumorHigh = length(out.tumorHigh.fov.accLbp) + 1;
        curTumorLow = length(out.tumorLow.fov.accLbp) + 1;
        %%
        
        if in == 1
            curSource = metaData.experiments.source1{iexp};
            curCellType = metaData.experiments.cellType1{iexp};
            tasksItr = 1 : metaData.experiments.n1{iexp};
        else
            curSource = metaData.experiments.source2{iexp};
            curCellType = metaData.experiments.cellType2{iexp};
            tasksItr = (metaData.experiments.n1{iexp} + 1) : (metaData.experiments.n1{iexp} + metaData.experiments.n2{iexp});
        end
        
        cellTypeInd = find(strcmpi(curCellType,metaData.cellTypes.ids));
        curMetEff = metaData.cellTypes.metastaticEfficiency(cellTypeInd);
        
        for itask = tasksItr
            
            % in exclude list
            if ismember(itask,metaData.experiments.exclude{iexp});
                continue;
            end
            
            lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                curFname '_s' sprintf('%02d',itask) filesep...
                'tracking' filesep 'lbpData.mat'];
            
            if ~exist(lbpFname,'file')
                lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                    curFname '_s' sprintf('%d',itask) filesep...
                    'tracking' filesep 'lbpData.mat'];
            end
            
            if ~exist(lbpFname,'file')
                warning(['LBP file' lbpFname ' (either 1 or 01) does not exist']);
                continue;
                error(['LBP file' lbpFname ' (either 1 or 01) does not exist']);
            end
            
            load(lbpFname); % pyramidLBP (?), accumulatedFovPyramidLBP
            
            %% Actual accumulation
            curLBPsFov = accumulatedFovPyramidLBP{iScale}; %#ok<USENS>
            curLBPsBck = accumulatedBckPyramidLBP{iScale}; %#ok<USENS>
            curLBPsFwd = accumulatedFwdPyramidLBP{iScale}; %#ok<USENS>
            
            if condInd == 1 % all                
                if length(out.allCells.fov.accLbp) < curAll
                    out.allCells.fov.accLbp{curAll} = [];
                    out.allCells.bck.accLbp{curAll} = [];
                    out.allCells.fwd.accLbp{curAll} = [];
                    % loaction 
                    out.allCells.fov.locations{curAll}.locationLbp = {};
                    out.allCells.bck.locations{curAll}.locationLbp = {};
                    out.allCells.fwd.locations{curAll}.locationLbp = {};
                    out.allCells.fov.locations{curAll}.locationStr = {};
                    out.allCells.bck.locations{curAll}.locationStr = {};
                    out.allCells.fwd.locations{curAll}.locationStr = {};
                end                                
                
                % loaction
                iLocation = length(out.allCells.fov.locations{curAll}.locationLbp) + 1;
                out.allCells.fov.locations{curAll}.locationLbp{iLocation} = curLBPsFov;
                out.allCells.bck.locations{curAll}.locationLbp{iLocation} = curLBPsBck;
                out.allCells.fwd.locations{curAll}.locationLbp{iLocation} = curLBPsFwd;
                out.allCells.fov.locations{curAll}.locationStr{iLocation} = sprintf('%02d',itask);
                out.allCells.bck.locations{curAll}.locationStr{iLocation} = sprintf('%02d',itask);
                out.allCells.fwd.locations{curAll}.locationStr{iLocation} = sprintf('%02d',itask);
                
                % accumulation
                out.allCells.fov.accLbp{curAll} = [out.allCells.fov.accLbp{curAll},curLBPsFov];
                out.allCells.bck.accLbp{curAll} = [out.allCells.bck.accLbp{curAll},curLBPsBck];
                out.allCells.fwd.accLbp{curAll} = [out.allCells.fwd.accLbp{curAll},curLBPsFwd];
                out.allCells.strs{curAll} = [curCellType '_' curDate]; % just repeats for every time
                out.allCells.date{curAll} = curDate;
                out.allCells.cellType{curAll} = curCellType;
                out.allCells.cellTypeInd{curAll} = cellTypeInd; % index in metaData.cellTypes.ids
                assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
                out.allCells.source{curAll} = curSource;
                out.allCells.metEff{curAll} = curMetEff;
            end
            
            %             if condInd == 2 % type
            %                 if length(out.cellTypes{cellTypeInd}.fov.accLbp) < curTypes(cellTypeInd)
            %                     out.cellTypes{cellTypeInd}.fov.accLbp{curTypes(cellTypeInd)} = [];
            %                     out.cellTypes{cellTypeInd}.bck.accLbp{curTypes(cellTypeInd)} = [];
            %                     out.cellTypes{cellTypeInd}.fwd.accLbp{curTypes(cellTypeInd)} = [];
            %                 end
            %                 out.cellTypes{cellTypeInd}.fov.accLbp{curTypes(cellTypeInd)} = [out.cellTypes{cellTypeInd}.fov.accLbp{curTypes(cellTypeInd)},curLBPsFov];
            %                 out.cellTypes{cellTypeInd}.bck.accLbp{curTypes(cellTypeInd)} = [out.cellTypes{cellTypeInd}.bck.accLbp{curTypes(cellTypeInd)},curLBPsBck];
            %                 out.cellTypes{cellTypeInd}.fwd.accLbp{curTypes(cellTypeInd)} = [out.cellTypes{cellTypeInd}.fwd.accLbp{curTypes(cellTypeInd)},curLBPsBck];
            %                 out.cellTypes{cellTypeInd}.strs{curAll} = [curCellType '_' curDate];
            %             end
            %
            %             if condInd == 3 % source
            %                 if strcmp(curSource,'Tumors')
            %                     if length(out.tumors.fov.accLbp) < curTumors
            %                         out.tumors.fov.accLbp{curTumors} = [];
            %                         out.tumors.bck.accLbp{curTumors} = [];
            %                         out.tumors.fwd.accLbp{curTumors} = [];
            %                     end
            %                     out.tumors.fov.accLbp{curTumors} = [out.tumors.fov.accLbp{curTumors},curLBPsFov];
            %                     out.tumors.bck.accLbp{curTumors} = [out.tumors.bck.accLbp{curTumors},curLBPsBck];
            %                     out.tumors.fwd.accLbp{curTumors} = [out.tumors.fwd.accLbp{curTumors},curLBPsFwd];
            %                     out.tumors.strs{curAll} = [curCellType '_' curDate];
            %                 else if strcmp(curSource,'CellLines')
            %                         if length(out.cellLines.fov.accLbp) < curCellLines
            %                             out.cellLines.fov.accLbp{curCellLines} = [];
            %                             out.cellLines.bck.accLbp{curCellLines} = [];
            %                             out.cellLines.fwd.accLbp{curCellLines} = [];
            %                         end
            %                         out.cellLines.fov.accLbp{curCellLines} = [out.cellLines.fov.accLbp{curCellLines},curLBPsFov];
            %                         out.cellLines.bck.accLbp{curCellLines} = [out.cellLines.bck.accLbp{curCellLines},curLBPsBck];
            %                         out.cellLines.fwd.accLbp{curCellLines} = [out.cellLines.fwd.accLbp{curCellLines},curLBPsFwd];
            %                         out.cellLines.strs{curAll} = [curCellType '_' curDate];
            %                     else if strcmp(curSource,'Melanocytes')
            %                             if length(out.melanocytes.fov.accLbp) < curMelanocytes
            %                                 out.melanocytes.fov.accLbp{curMelanocytes} = [];
            %                                 out.melanocytes.bck.accLbp{curMelanocytes} = [];
            %                                 out.melanocytes.fwd.accLbp{curMelanocytes} = [];
            %                             end
            %                             out.melanocytes.fov.accLbp{curMelanocytes} = [out.melanocytes.fov.accLbp{curMelanocytes},curLBPsFov];
            %                             out.melanocytes.bck.accLbp{curMelanocytes} = [out.melanocytes.bck.accLbp{curMelanocytes},curLBPsBck];
            %                             out.melanocytes.fwd.accLbp{curMelanocytes} = [out.melanocytes.fwd.accLbp{curMelanocytes},curLBPsFwd];
            %                             out.melanocytes.strs{curAll} = [curCellType '_' curDate];
            %                         end
            %                     end
            %                 end
            %             end
            %
            %             if condInd == 4 % metastatic
            %                 if ~isnan(curMetEff)
            %                     if curMetEff
            %                         if length(out.tumorHigh.fov.accLbp) < curTumorHigh
            %                             out.tumorHigh.fov.accLbp{curTumorHigh} = [];
            %                             out.tumorHigh.bck.accLbp{curTumorHigh} = [];
            %                             out.tumorHigh.fwd.accLbp{curTumorHigh} = [];
            %                         end
            %                         out.tumorHigh.fov.accLbp{curTumorHigh} = [out.tumorHigh.fov.accLbp{curTumorHigh},curLBPsFov];
            %                         out.tumorHigh.bck.accLbp{curTumorHigh} = [out.tumorHigh.bck.accLbp{curTumorHigh},curLBPsBck];
            %                         out.tumorHigh.fwd.accLbp{curTumorHigh} = [out.tumorHigh.fwd.accLbp{curTumorHigh},curLBPsFwd];
            %                         out.tumorHigh.strs{curAll} = [curCellType '_' curDate];
            %                     else
            %                         if length(out.tumorLow.fov.accLbp) < curTumorLow
            %                             out.tumorLow.fov.accLbp{curTumorLow} = [];
            %                             out.tumorLow.bck.accLbp{curTumorLow} = [];
            %                             out.tumorLow.fwd.accLbp{curTumorLow} = [];
            %                         end
            %                         out.tumorLow.fov.accLbp{curTumorLow} = [out.tumorLow.fov.accLbp{curTumorLow},curLBPsFov];
            %                         out.tumorLow.bck.accLbp{curTumorLow} = [out.tumorLow.bck.accLbp{curTumorLow},curLBPsBck];
            %                         out.tumorLow.fwd.accLbp{curTumorLow} = [out.tumorLow.fwd.accLbp{curTumorLow},curLBPsBck];
            %                         out.tumorLow.strs{curAll} = [curCellType '_' curDate];
            %                     end
            %                 end
            %             end
        end % task
    end % first or second cell type in well
end % number of experiments

end

%% Use precalculated PCA
function [fovpca, bckpca, fwdpca] = ...
    getPcaGeneral(pcaStatsGeneralPrefix,iScale)

load([pcaStatsGeneralPrefix 'accumulatedLBP_' num2str(iScale) '_fov_all.mat']);
fovpca.limitPercentile = 2;
fovpca.meanLbp = allCellsFov.meansLbp;
fovpca.stdLbp = allCellsFov.stdsLbp;
fovpca.normLbp = allCellsFov.accLbpNorm; % normalizing
fovpca.coeffLbp = allCellsFov.accLbpCoeff;
fovpca.scoreLbp = allCellsFov.accLbpScore;
fovpca.latentLbp = allCellsFov.accLbpLatent;
fovpca.PCLimitLow = prctile(fovpca.scoreLbp,fovpca.limitPercentile,1);
fovpca.PCLimitHigh = prctile(fovpca.scoreLbp,100-fovpca.limitPercentile,1);
fovpca.low = prctile(fovpca.normLbp(:),fovpca.limitPercentile);
fovpca.high = prctile(fovpca.normLbp(:),100-fovpca.limitPercentile);
fovpca.lowCoeff = prctile(fovpca.coeffLbp(:),fovpca.limitPercentile);
fovpca.highCoeff = prctile(fovpca.coeffLbp(:),100-fovpca.limitPercentile);
clear out;

load([pcaStatsGeneralPrefix 'accumulatedLBP_' num2str(iScale) '_bck_all.mat']);
bckpca.limitPercentile = 2;
bckpca.meanLbp = allCellsBck.meansLbp;
bckpca.stdLbp = allCellsBck.stdsLbp;
bckpca.normLbp = allCellsBck.accLbpNorm; % normalizing
bckpca.coeffLbp = allCellsBck.accLbpCoeff;
bckpca.scoreLbp = allCellsBck.accLbpScore;
bckpca.latentLbp = allCellsBck.accLbpLatent;
bckpca.PCLimitLow = prctile(bckpca.scoreLbp,bckpca.limitPercentile,1);
bckpca.PCLimitHigh = prctile(bckpca.scoreLbp,100-bckpca.limitPercentile,1);
bckpca.low = prctile(bckpca.normLbp(:),bckpca.limitPercentile);
bckpca.high = prctile(bckpca.normLbp(:),100-bckpca.limitPercentile);
bckpca.lowCoeff = prctile(bckpca.coeffLbp(:),bckpca.limitPercentile);
bckpca.highCoeff = prctile(bckpca.coeffLbp(:),100-bckpca.limitPercentile);
clear out;

load([pcaStatsGeneralPrefix 'accumulatedLBP_' num2str(iScale) '_fwd_all.mat']);
fwdpca.limitPercentile = 2;
fwdpca.meanLbp = allCellsFwd.meansLbp;
fwdpca.stdLbp = allCellsFwd.stdsLbp;
fwdpca.normLbp = allCellsFwd.accLbpNorm; % normalizing
fwdpca.coeffLbp = allCellsFwd.accLbpCoeff;
fwdpca.scoreLbp = allCellsFwd.accLbpScore;
fwdpca.latentLbp = allCellsFwd.accLbpLatent;
fwdpca.PCLimitLow = prctile(fwdpca.scoreLbp,fwdpca.limitPercentile,1);
fwdpca.PCLimitHigh = prctile(fwdpca.scoreLbp,100-fwdpca.limitPercentile,1);
fwdpca.low = prctile(fwdpca.normLbp(:),fwdpca.limitPercentile);
fwdpca.high = prctile(fwdpca.normLbp(:),100-fwdpca.limitPercentile);
fwdpca.lowCoeff = prctile(fwdpca.coeffLbp(:),fwdpca.limitPercentile);
fwdpca.highCoeff = prctile(fwdpca.coeffLbp(:),100-fwdpca.limitPercentile);
clear out;

end