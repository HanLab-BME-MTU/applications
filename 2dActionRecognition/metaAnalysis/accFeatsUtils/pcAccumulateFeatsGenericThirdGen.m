%% Generic accumulation of features from cell type + condition (over all locations)
%   Params:
%        featsInDname   % input directory
%        featsOutDname  % output directory
%        metaDataFname  % meta data (from excel DB)
%        fGetFeats      % rule for extracting the features (e.g.,different time points, combination of features) 
%                                   also retruns the ID of cells (from the raw feature file) for cell explorer
%        featsStrIn     % feature string for the input files
%        featsStrOut    % feature string for the output files

function [] = pcAccumulateFeatsGenericThirdGen(params)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;
always = false;

featsInDname = params.featsInDname;
featsOutDname = params.featsOutDname;
featsStrIn = params.featsStrIn; % prefix for the input of the features
featsStrOut = params.featsStrOut; % prefix for output of the features
metaDataFname = params.metaDataFname; % meta data (from excel DB)
fGetFeats = params.fGetFeats;

nScales = 4; % 4 scales (1 to 1/8)
scales = 1.0./2.^((1:nScales)-1);

if ~exist(featsOutDname,'dir')
    unix(sprintf('mkdir %s',featsOutDname));
end

load(metaDataFname);%metaData

for iScale = 1 : nScales
    
    % features
    accFeatsFnameAllFov = [featsOutDname filesep featsStrOut '_' num2str(iScale) '_fov_all.mat'];
    
    %% Accumulate
    if ~exist(accFeatsFnameAllFov,'file') || always
        tic;
        out = accumulateFeatsGenericThirdGen(metaData,featsInDname,featsStrIn,iScale,fGetFeats);
        allCellsFov = out.allCells.fov;
        allInfo = out.allCells;
        clear out;
        save(accFeatsFnameAllFov,'allCellsFov','-v7.3');
        save([featsOutDname filesep featsStrOut '_' num2str(iScale) '_allInfo.mat'],'allInfo','-v7.3');
        clear allCellsFov;
        tt = toc;
        fprintf(sprintf('done %d %s (%d min.)\n',iScale,'all',round(tt/60)));
    end
end
end

%% Here we do the accumulation for every well!
function out = accumulateFeatsGenericThirdGen(metaData,featsInDname,featsStrIn,iScale,fGetFeats)

%% init
out.allCells.fov.accFeats = {};
out.allCells.fov.accLocation = {};
out.allCells.fov.accCellID = {};
                
out.allCells.strs = {};
out.allCells.date = {};
out.allCells.task = {}; % needed to go back to the cell movie for cell explorer!
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

%% Accumulation
for iexp = 1 : metaData.experiments.N
    expFname = metaData.experiments.fnames{iexp};
    expDate = expFname(1:6);
    %     curSources = metaData.experiments.source{iexp};
    
    expCellTypes = metaData.experiments.cellTypes{iexp};
    uniqueExpCellTypes = unique(expCellTypes);
    nExpCellTypes = length(uniqueExpCellTypes);
    
    for iCellType = 1 : nExpCellTypes
        
        curCellType = uniqueExpCellTypes(iCellType);
        tasksItr = find(strcmp(curCellType,expCellTypes));
        
        %% count the next open spot for each condition
        nextPosition = length(out.allCells.fov.accFeats) + 1;
        
        for itask = tasksItr
            
            % in exclude list
            if ismember(itask,metaData.experiments.exclude{iexp});
                continue;
            end
                        
            cellTypeInd = find(strcmpi(curCellType,metaData.cellTypes.ids));
            curSource = metaData.cellTypes.source{cellTypeInd};
            curMetEff = metaData.cellTypes.metastaticEfficiency(cellTypeInd);
            
            featsInFname = [featsInDname num2str(iScale) filesep expFname sprintf('_s%02d_%s.mat',itask,featsStrIn)];
            if ~exist(featsInFname,'file')
                error('file %s does not exist',featsInFname);
            end
            
            cellData = getCellData(featsInFname);
            [curFeats,curIDs] = getLchFeats(fGetFeats,cellData);
            nCurCells = size(curFeats,2);
            
            
            %% Actual accumulation
            
            if length(out.allCells.fov.accFeats) < nextPosition
                out.allCells.fov.accFeats{nextPosition} = [];
                out.allCells.fov.accLocation{nextPosition} = [];
                out.allCells.fov.accCellID{nextPosition} = [];
                
                % loaction
                out.allCells.fov.locations{nextPosition}.locationDeltaLbp = {};
                out.allCells.fov.locations{nextPosition}.locationStr = {};
            end
            
            % loaction
            iLocation = length(out.allCells.fov.locations{nextPosition}.locationDeltaLbp) + 1;
            out.allCells.fov.locations{nextPosition}.locationDeltaLbp{iLocation} = curFeats;
            out.allCells.fov.locations{nextPosition}.locationStr{iLocation} = sprintf('%02d',itask);
            
            % accumulation
            out.allCells.fov.accFeats{nextPosition} = [out.allCells.fov.accFeats{nextPosition},curFeats];
            out.allCells.fov.accLocation{nextPosition} = [out.allCells.fov.accLocation{nextPosition},ones(1,nCurCells)*itask];
            out.allCells.fov.accCellID{nextPosition} = [out.allCells.fov.accCellID{nextPosition},curIDs];
            
            out.allCells.strs{nextPosition} = [curCellType '_' expDate]; % just repeats for every time
            out.allCells.date{nextPosition} = expDate;            
            out.allCells.cellType{nextPosition} = curCellType;
            out.allCells.cellTypeInd{nextPosition} = cellTypeInd; % index in metaData.cellTypes.ids
            assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
            out.allCells.source{nextPosition} = curSource;
            out.allCells.metEff{nextPosition} = curMetEff;
        end % task
    end % cell type within experiment
end % number of experiments

end
