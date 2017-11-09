%% Accumulates LBP for taks by experiment (well) for each cell
% (lbp vector for every cell at each time point, grouped by cell):
% all data (no clusters for cell type, source or metastatic efficiency)
% The main change is that the accumulation occurs at the single cell level.
% 

% allCells include all the information that is needed for control-analysis
% / classifications / distnace maps

% Assaf Zaritsky, Nov. 2016

function [] = pcAccumulateLBPWellCell()

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

always = true;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBPWellCell'];
cellLbpPrefix = [lbpDirname filesep 'cellLBP_'];
metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);%metaData

for iScale = 1 : nScales
    
    % these are different from the pcAccumulateLBPNew, holding an array of
    % features, one per well
    cellLbpFnameAllFov = [cellLbpPrefix num2str(iScale) '_fov_all.mat'];       
    cellLbpFnameAllBck = [cellLbpPrefix num2str(iScale) '_bck_all.mat'];        
    cellLbpFnameAllFwd = [cellLbpPrefix num2str(iScale) '_fwd_all.mat'];    
    
    pcaStatsGeneralPrefix = [lbpDirname '/../LBP/'];
    
    %% Accumulate
    if ~exist(cellLbpFnameAllFwd,'file') || always
        tic;
        out = accumulateLBPWellPerCell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,'all');
        allCellsFov = out.allCells.fov;
        allCellsBck = out.allCells.bck;
        allCellsFwd = out.allCells.fwd;
        rmfield(out.allCells,'fov');
        rmfield(out.allCells,'bck');
        rmfield(out.allCells,'fwd');
        allInfo = out.allCells;
        % do FWD
        clear out;
        save(cellLbpFnameAllFov,'allCellsFov','-v7.3');
        save(cellLbpFnameAllBck,'allCellsBck','-v7.3');
        save(cellLbpFnameAllFwd,'allCellsFwd','-v7.3');
        save([cellLbpPrefix num2str(iScale) '_allInfo.mat'],'allInfo','-v7.3');
        clear allCellsFov allCellsBck allCellsFwd;
        tt = toc;
        fprintf(sprintf('done %d %s (%d min.)\n',iScale,'all',round(tt/60)));
    end    
end
end

%% Here we do the accumulation for every well!
function out = accumulateLBPWellPerCell(metaData,analysisDirname,pcaStatsGeneralPrefix,iScale,strLabel)

%% PCA data is taken from the single cell, single time point data! (only general)
[out.allCells.fov.pca, out.allCells.bck.pca, out.allCells.fwd.pca] = ...
    getPcaGeneral(pcaStatsGeneralPrefix,iScale); 


%% init
out.allCells.fov.cellLbp = {};
out.allCells.bck.cellLbp = {};
out.allCells.fwd.cellLbp = {};
out.allCells.strs = {};
out.allCells.date = {};
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

condInd = find(strcmp(strLabel,{'all','type','source','metastatic'}));

%% Accumulation
for iexp = 1 : metaData.experiments.N
    curFname = metaData.experiments.fnames{iexp};
    curDate = curFname(1:6);
    for in = 1 : 2
        %% count the next open spot for each condition
        curAll = length(out.allCells.fov.cellLbp) + 1; % # of the current experiment / well
        
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
            
            load(lbpFname); % lbpData
            
            trajectoriesT0 = lbpData.t0; % NEW
            % lbpData.fov{136 - cell #}.pyramidLBP{scale #}
            
            %% Actual accumulation
            curLBPsFov = accumulatedFovPyramidLBP{iScale}; %#ok<USENS>
            curLBPsBck = accumulatedBckPyramidLBP{iScale}; %#ok<USENS>
            curLBPsFwd = accumulatedFwdPyramidLBP{iScale}; %#ok<USENS>
            
            if condInd == 1 % all                
                if length(out.allCells.fov.cellLbp) < curAll
                    out.allCells.fov.cellLbp{curAll} = [];
                    out.allCells.bck.cellLbp{curAll} = [];
                    out.allCells.fwd.cellLbp{curAll} = [];
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
                out.allCells.fov.cellLbp{curAll} = [out.allCells.fov.cellLbp{curAll},curLBPsFov];
                out.allCells.bck.cellLbp{curAll} = [out.allCells.bck.cellLbp{curAll},curLBPsBck];
                out.allCells.fwd.cellLbp{curAll} = [out.allCells.fwd.cellLbp{curAll},curLBPsFwd];
                out.allCells.strs{curAll} = [curCellType '_' curDate]; % just repeats for every time
                out.allCells.date{curAll} = curDate;
                out.allCells.cellType{curAll} = curCellType;
                out.allCells.cellTypeInd{curAll} = cellTypeInd; % index in metaData.cellTypes.ids
                assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
                out.allCells.source{curAll} = curSource;
                out.allCells.metEff{curAll} = curMetEff;
            end           
        end % task
    end % first or second cell type in well
end % number of experiments

end

%% Use precalculated PCA
function [fovpca, bckpca, fwdpca] = ...
    getPcaGeneral(pcaStatsGeneralPrefix,iScale)

load([pcaStatsGeneralPrefix 'cellLBP_' num2str(iScale) '_fov_all.mat']);
fovpca.limitPercentile = 2;
fovpca.meanLbp = allCellsFov.meansLbp;
fovpca.stdLbp = allCellsFov.stdsLbp;
fovpca.normLbp = allCellsFov.cellLbpNorm; % normalizing
fovpca.coeffLbp = allCellsFov.cellLbpCoeff;
fovpca.scoreLbp = allCellsFov.cellLbpScore;
fovpca.latentLbp = allCellsFov.cellLbpLatent;
fovpca.PCLimitLow = prctile(fovpca.scoreLbp,fovpca.limitPercentile,1);
fovpca.PCLimitHigh = prctile(fovpca.scoreLbp,100-fovpca.limitPercentile,1);
fovpca.low = prctile(fovpca.normLbp(:),fovpca.limitPercentile);
fovpca.high = prctile(fovpca.normLbp(:),100-fovpca.limitPercentile);
fovpca.lowCoeff = prctile(fovpca.coeffLbp(:),fovpca.limitPercentile);
fovpca.highCoeff = prctile(fovpca.coeffLbp(:),100-fovpca.limitPercentile);
clear out;

load([pcaStatsGeneralPrefix 'cellLBP_' num2str(iScale) '_bck_all.mat']);
bckpca.limitPercentile = 2;
bckpca.meanLbp = allCellsBck.meansLbp;
bckpca.stdLbp = allCellsBck.stdsLbp;
bckpca.normLbp = allCellsBck.cellLbpNorm; % normalizing
bckpca.coeffLbp = allCellsBck.cellLbpCoeff;
bckpca.scoreLbp = allCellsBck.cellLbpScore;
bckpca.latentLbp = allCellsBck.cellLbpLatent;
bckpca.PCLimitLow = prctile(bckpca.scoreLbp,bckpca.limitPercentile,1);
bckpca.PCLimitHigh = prctile(bckpca.scoreLbp,100-bckpca.limitPercentile,1);
bckpca.low = prctile(bckpca.normLbp(:),bckpca.limitPercentile);
bckpca.high = prctile(bckpca.normLbp(:),100-bckpca.limitPercentile);
bckpca.lowCoeff = prctile(bckpca.coeffLbp(:),bckpca.limitPercentile);
bckpca.highCoeff = prctile(bckpca.coeffLbp(:),100-bckpca.limitPercentile);
clear out;

load([pcaStatsGeneralPrefix 'cellLBP_' num2str(iScale) '_fwd_all.mat']);
fwdpca.limitPercentile = 2;
fwdpca.meanLbp = allCellsFwd.meansLbp;
fwdpca.stdLbp = allCellsFwd.stdsLbp;
fwdpca.normLbp = allCellsFwd.cellLbpNorm; % normalizing
fwdpca.coeffLbp = allCellsFwd.cellLbpCoeff;
fwdpca.scoreLbp = allCellsFwd.cellLbpScore;
fwdpca.latentLbp = allCellsFwd.cellLbpLatent;
fwdpca.PCLimitLow = prctile(fwdpca.scoreLbp,fwdpca.limitPercentile,1);
fwdpca.PCLimitHigh = prctile(fwdpca.scoreLbp,100-fwdpca.limitPercentile,1);
fwdpca.low = prctile(fwdpca.normLbp(:),fwdpca.limitPercentile);
fwdpca.high = prctile(fwdpca.normLbp(:),100-fwdpca.limitPercentile);
fwdpca.lowCoeff = prctile(fwdpca.coeffLbp(:),fwdpca.limitPercentile);
fwdpca.highCoeff = prctile(fwdpca.coeffLbp(:),100-fwdpca.limitPercentile);
clear out;

end