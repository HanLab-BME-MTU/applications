%% Generic accumulation of features from cell type + condition (over all locations)
%   Params:
%        featsInDname   % input directory
%        featsOutDname  % output directory
%        metaDataFname  % meta data (from excel DB)
%        fGetFeats      % rule for extracting the features (e.g.,different time points, combination of features)
%                                   also retruns the ID of cells (from the raw feature file) for cell explorer
%        featsStrIn     % feature string for the input files
%        featsStrOut    % feature string for the output files

function [] = pcAccumulateFeatsGeneric2018(params)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;
always = true;%false;

featsInDname = params.featsInDname;
featsOutDname = params.featsOutDname;
featsStrIn = params.featsStrIn; % prefix for the input of the features
featsStrOut = params.featsStrOut; % prefix for output of the features
featsStrID = params.featsStrID; % prefix for feature ID  (used in cell explorer)
metaDataFname = params.metaDataFname; % meta data (from excel DB)
fGetFeats = params.fGetFeats;

if ~exist(featsOutDname,'dir')
    unix(sprintf('mkdir %s',featsOutDname));
end

load(metaDataFname);%metaData

% features
accFeatsFnameAll = [featsOutDname filesep featsStrOut '_' num2str(iScale) '_all.mat'];

%% Accumulate
if ~exist(accFeatsFnameAll,'file') || always
    tic;
    out = accumulateFeatsGenericThirdGen(metaData,featsInDname,featsStrIn,featsStrID,iScale,fGetFeats);    
    allCells = out.allCells;
    clear out;    
    save(accFeatsFnameAll,'allCells','-v7.3');    
    tt = toc;    
end
end

%% Here we do the accumulation for every well!
function out = accumulateFeatsGenericThirdGen(metaData,featsInDname,featsStrIn,featsStrID,iScale,fGetFeats)

%% init
out.allCells.accFeats = {};
out.allCells.accLocation = {};
out.allCells.accCellID = {};

out.allCells.strs = {};
out.allCells.date = {};
out.allCells.task = {}; % needed to go back to the cell movie for cell explorer!
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

out.allCells.expStr = {};
out.allCells.featsStrID = [featsStrID '_' num2str(iScale)];

%% Accumulation
for iexp = 1 : metaData.experiments.N
    expFname = metaData.experiments.fnames{iexp};
    expDate = expFname(1:6);
    %     curSources = metaData.experiments.source{iexp};
    
    expCellTypes = metaData.experiments.cellTypes{iexp};
    uniqueExpCellTypes = unique(expCellTypes);
    nExpCellTypes = length(uniqueExpCellTypes);
    
    for iCellType = 1 : nExpCellTypes
        
        curCellType = uniqueExpCellTypes{iCellType};
        tasksItr = find(strcmp(curCellType,expCellTypes));
        
        %% count the next open spot for each condition
        nextPosition = length(out.allCells.accFeats) + 1;
        
        for itask = tasksItr
            
            % in exclude list
            if ismember(itask,metaData.experiments.exclude{iexp})
                continue;
            end
            
            cellTypeInd = find(strcmpi(curCellType,metaData.cellTypes.ids));
            curSource = metaData.cellTypes.source{cellTypeInd};
            curMetEff = metaData.cellTypes.metastaticEfficiency(cellTypeInd);
            
            featsInFname = [featsInDname num2str(iScale) filesep expFname sprintf('_s%02d_%s.mat',itask,featsStrIn)];
            if ~exist(featsInFname,'file')
                error('file %s does not exist',featsInFname);
            end
            
            [curFeats,curIDs,curTXY] = getLchFeats(fGetFeats,dataTask);
            nCurCells = size(curFeats,2);
            
            
            %% Actual accumulation
            
            if length(out.allCells.accFeats) < nextPosition
                out.allCells.accFeats{nextPosition} = [];
                out.allCells.accLocation{nextPosition} = [];
                out.allCells.accCellID{nextPosition} = [];
                
                % loaction
                out.allCells.locations{nextPosition}.locationFeats = {};
                out.allCells.locations{nextPosition}.locationStr = {};
                out.allCells.locations{nextPosition}.locationTXY = {};
            end
            
            % loaction
            iLocation = length(out.allCells.locations{nextPosition}.locationFeats) + 1;
            out.allCells.locations{nextPosition}.locationFeats{iLocation} = curFeats;
            out.allCells.locations{nextPosition}.locationStr{iLocation} = sprintf('%02d',itask);
            out.allCells.locations{nextPosition}.locationTXY{iLocation} = curTXY;
            
            % accumulation
            out.allCells.accFeats{nextPosition} = [out.allCells.accFeats{nextPosition},curFeats]; %% TODO??
            
            out.allCells.accLocation{nextPosition} = [out.allCells.accLocation{nextPosition},ones(1,nCurCells)*itask];
            out.allCells.accCellID{nextPosition} = [out.allCells.accCellID{nextPosition},curIDs];
            
            out.allCells.strs{nextPosition} = [curCellType '_' expDate]; % just repeats for every time
            out.allCells.date{nextPosition} = expDate;
            out.allCells.cellType{nextPosition} = curCellType;
            out.allCells.cellTypeInd{nextPosition} = cellTypeInd; % index in metaData.cellTypes.ids
            assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
            out.allCells.source{nextPosition} = curSource;
            out.allCells.metEff{nextPosition} = curMetEff;
            
            out.allCells.expStr{nextPosition} = expFname;
        end % task
    end % cell type within experiment
end % number of experiments

end
