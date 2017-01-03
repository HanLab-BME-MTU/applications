%% Accumulates dLBP for tasks by experiment (well)
% median dLBP for every cell in every well
% all data, cluster by cell type, cluster by Tumor/Cell Line/Melanocyte,
% cluster by high/low metastatic efficiency

% WHAT DOES THIS MEAN? (allCells include all the information that is needed for control-analysis
% / classifications / distnace maps)

% Assaf Zaritsky, June. 2016

function [] = pcAccumulateDeltaLBPWell()

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

always = false;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/dLBP'];
accLbpPrefix = [lbpDirname filesep 'accumulatedDeltaLBP_'];
metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);%metaData

for iScale = 1 : nScales
    
    % these are different from the pcAccumulateLBPNew, holding an array of
    % features, one per well
    accLbpFnameAllFov = [accLbpPrefix num2str(iScale) '_fov_all.mat'];
    accLbpFnameSourceFov = [accLbpPrefix num2str(iScale) '_fov_source.mat'];
    accLbpFnameMetastaticFov = [accLbpPrefix num2str(iScale) '_fov_metastatic.mat'];
    accLbpFnameCellTypeFov = [accLbpPrefix num2str(iScale) '_fov_type.mat'];        
    
    %% Accumulate
    if ~exist(accLbpFnameAllFov,'file') || always
        tic;
        out = accumulateDeltaLBPWell(metaData,analysisDirname,iScale,'all');
        allCellsFov = out.allCells.fov;        
        allInfo = out.allCells;        
        clear out;
        save(accLbpFnameAllFov,'allCellsFov','-v7.3');        
        save([accLbpPrefix num2str(iScale) '_allInfo.mat'],'allInfo','-v7.3');
        clear allCellsFov;
        tt = toc;
        fprintf(sprintf('done %d %s (%d min.)\n',iScale,'all',round(tt/60)));
    end
    
    if ~exist(accLbpFnameCellTypeFov,'file') || always
        tic;
        out = accumulateDeltaLBPWell(metaData,analysisDirname,iScale,'type');
        cellTypesFov = cell(1,length(out.cellTypes));        
        cellTypesStr = cell(1,length(out.cellTypes));        
        for ict = 1 : length(out.cellTypes)
            cellTypesFov{ict} = out.cellTypes{ict}.fov; % remember that now accLbp has multiple entries, one for each experiment            
            if ~isempty(cellTypesFov{ict}.accDeltaLbp)
                cellTypesStr{ict} = out.cellTypes{ict}.strs{:};
            end
        end
        
        clear out;
        save(accLbpFnameCellTypeFov,'cellTypesFov','cellTypesStr','-v7.3');
        clear cellTypesFov        
        tt = toc;
        fprintf(sprintf('done %d %s\n',iScale,'type',round(tt/60)));
    end
    
    if ~exist(accLbpFnameSourceFov,'file') || always
        tic;
        out = accumulateDeltaLBPWell(metaData,analysisDirname,iScale,'source');
        melanocytesFov = out.melanocytes.fov;
        cellLinesFov = out.cellLines.fov;
        tumorsFov = out.tumors.fov;        
        clear out;
        save(accLbpFnameSourceFov,'melanocytesFov','cellLinesFov','tumorsFov','-v7.3');
        clear melanocytesFov cellLinesFov tumorsFov;        
        tt = toc;
        fprintf(sprintf('done %d %s\n',iScale,'source',round(tt/60)));
    end
    
    if ~exist(accLbpFnameMetastaticFov,'file') || always
        tic;
        out = accumulateDeltaLBPWell(metaData,analysisDirname,iScale,'metastatic');
        tumorHighFov = out.tumorHigh.fov;
        tumorLowFov = out.tumorLow.fov;        
        clear out;
        save(accLbpFnameMetastaticFov,'tumorHighFov','tumorLowFov','-v7.3');
        clear tumorHighFov tumorLowFov;        
        tt = toc;
        fprintf(sprintf('done %d %s\n',iScale,'metastatic',round(tt/60)));
    end
    % end
end
end

%% Here we do the accumulation for every well!
function out = accumulateDeltaLBPWell(metaData,analysisDirname,iScale,strLabel)

dLbpDname = [analysisDirname 'Cells/dLBP/'];

%% init
nCellTypes = length(metaData.cellTypes.ids);
out.allCells.fov.accDeltaLbp = {};
out.allCells.strs = {};
out.allCells.date = {};
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

out.cellTypes = cell(1,nCellTypes);


out.melanocytes.fov.accDeltaLbp = {};
out.cellLines.fov.accDeltaLbp = {};
out.tumors.fov.accDeltaLbp = {};
out.melanocytes.strs = {};
out.cellLines.strs = {};
out.tumors.strs = {};

out.tumorHigh.fov.accDeltaLbp = {};
out.tumorLow.fov.accDeltaLbp = {};
out.tumorHigh.strs = {};
out.tumorLow.strs = {};

for iCellType = 1 : nCellTypes
    out.cellTypes{iCellType}.fov.accDeltaLbp = {};    
    out.cellTypes{iCellType}.strs = {};
end

condInd = find(strcmp(strLabel,{'all','type','source','metastatic'}));

%% Accumulation
for iexp = 1 : metaData.experiments.N
    curFname = metaData.experiments.fnames{iexp};
    curDate = curFname(1:6);
    for in = 1 : 2
        %% count the next open spot for each condition
        curAll = length(out.allCells.fov.accDeltaLbp) + 1;
        
        curTypes = nan(1,nCellTypes);
        for iCellType = 1 : nCellTypes
            curTypes(iCellType) = length(out.cellTypes{iCellType}.fov.accDeltaLbp) + 1;
        end
        
        curMelanocytes = length(out.melanocytes.fov.accDeltaLbp) + 1;
        curCellLines = length(out.cellLines.fov.accDeltaLbp) + 1;
        curTumors = length(out.tumors.fov.accDeltaLbp) + 1;
        
        curTumorHigh = length(out.tumorHigh.fov.accDeltaLbp) + 1;
        curTumorLow = length(out.tumorLow.fov.accDeltaLbp) + 1;
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
            
            deltaLbpFname = [dLbpDname num2str(iScale) filesep curFname sprintf('_s%02d_dLBP.mat',itask)];
            if ~exist(deltaLbpFname,'file')
                error('dLBP file %s does not exist',[curFname sprintf('_s%02d_dLBP.mat',itask)]);
            end
                        
            load(deltaLbpFname); % dLbpWell: 1 x nCells, each had deltaLbpMedian field - median dLBP for a cell 
            curDeltaLBPFov = getDeltaLbps(dLbpWell);
            
            
            %% Actual accumulation
            if condInd == 1 % all                
                if length(out.allCells.fov.accDeltaLbp) < curAll
                    out.allCells.fov.accDeltaLbp{curAll} = [];                    
                    % loaction 
                    out.allCells.fov.locations{curAll}.locationDeltaLbp = {};                    
                    out.allCells.fov.locations{curAll}.locationStr = {};                    
                end                                
                
                % loaction
                iLocation = length(out.allCells.fov.locations{curAll}.locationDeltaLbp) + 1;
                out.allCells.fov.locations{curAll}.locationDeltaLbp{iLocation} = curDeltaLBPFov;                
                out.allCells.fov.locations{curAll}.locationStr{iLocation} = sprintf('%02d',itask);                
                
                % accumulation
                out.allCells.fov.accDeltaLbp{curAll} = [out.allCells.fov.accDeltaLbp{curAll},curDeltaLBPFov];                
                out.allCells.strs{curAll} = [curCellType '_' curDate]; % just repeats for every time
                out.allCells.date{curAll} = curDate;
                out.allCells.cellType{curAll} = curCellType;
                out.allCells.cellTypeInd{curAll} = cellTypeInd; % index in metaData.cellTypes.ids
                assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
                out.allCells.source{curAll} = curSource;
                out.allCells.metEff{curAll} = curMetEff;
            end
            
            if condInd == 2 % type
                if length(out.cellTypes{cellTypeInd}.fov.accDeltaLbp) < curTypes(cellTypeInd)
                    out.cellTypes{cellTypeInd}.fov.accDeltaLbp{curTypes(cellTypeInd)} = [];
                end
                out.cellTypes{cellTypeInd}.fov.accDeltaLbp{curTypes(cellTypeInd)} = [out.cellTypes{cellTypeInd}.fov.accDeltaLbp{curTypes(cellTypeInd)},curDeltaLBPFov];
                %                 out.cellTypes{cellTypeInd}.strs{curAll} = [curCellType '_' curDate];
                out.cellTypes{cellTypeInd}.strs{curAll} = curCellType;
            end
            
            if condInd == 3 % source
                if strcmp(curSource,'Tumors')
                    if length(out.tumors.fov.accDeltaLbp) < curTumors
                        out.tumors.fov.accDeltaLbp{curTumors} = [];
                    end
                    out.tumors.fov.accDeltaLbp{curTumors} = [out.tumors.fov.accDeltaLbp{curTumors},curDeltaLBPFov];
                    out.tumors.strs{curAll} = [curCellType '_' curDate];
                else if strcmp(curSource,'CellLines')
                        if length(out.cellLines.fov.accDeltaLbp) < curCellLines
                            out.cellLines.fov.accDeltaLbp{curCellLines} = [];
                        end
                        out.cellLines.fov.accDeltaLbp{curCellLines} = [out.cellLines.fov.accDeltaLbp{curCellLines},curDeltaLBPFov];
                        out.cellLines.strs{curAll} = [curCellType '_' curDate];
                    else if strcmp(curSource,'Melanocytes')
                            if length(out.melanocytes.fov.accDeltaLbp) < curMelanocytes
                                out.melanocytes.fov.accDeltaLbp{curMelanocytes} = [];
                            end
                            out.melanocytes.fov.accDeltaLbp{curMelanocytes} = [out.melanocytes.fov.accDeltaLbp{curMelanocytes},curDeltaLBPFov];
                            out.melanocytes.strs{curAll} = [curCellType '_' curDate];
                        end
                    end
                end
            end
            
            if condInd == 4 % metastatic
                if ~isnan(curMetEff)
                    if curMetEff
                        if length(out.tumorHigh.fov.accDeltaLbp) < curTumorHigh
                            out.tumorHigh.fov.accDeltaLbp{curTumorHigh} = [];
                        end
                        out.tumorHigh.fov.accDeltaLbp{curTumorHigh} = [out.tumorHigh.fov.accDeltaLbp{curTumorHigh},curDeltaLBPFov];
                        out.tumorHigh.strs{curAll} = [curCellType '_' curDate];
                    else
                        if length(out.tumorLow.fov.accDeltaLbp) < curTumorLow
                            out.tumorLow.fov.accDeltaLbp{curTumorLow} = [];
                        end
                        out.tumorLow.fov.accDeltaLbp{curTumorLow} = [out.tumorLow.fov.accDeltaLbp{curTumorLow},curDeltaLBPFov];
                        out.tumorLow.strs{curAll} = [curCellType '_' curDate];
                    end
                end
            end
        end % task
    end % first or second cell type in well
end % number of experiments

end

%%
function deltaLBPs = getDeltaLbps(dLbpWell)
n = length(dLbpWell);
deltaLBPs = nan(1,n);
for i = 1 : n
    deltaLBPs(i) = dLbpWell{i}.deltaLbpMedian;
end
end