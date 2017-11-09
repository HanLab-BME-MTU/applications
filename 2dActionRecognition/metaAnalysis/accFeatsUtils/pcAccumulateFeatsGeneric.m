%% Generic accumulation of features from cell type + condition (over all locations)
%   Params: 
%        featsInDname   % input directory
%        featsOutDname  % output directory
%        metaDataFname  % meta data (from excel DB)
%        fGetFeats      % rule for extracting the features (e.g.,different time points, combination of features)
%        featsStrIn     % feature string for the input files
%        featsStrOut    % feature string for the output files

function [] = pcAccumulateFeatsGeneric(params)

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
        out = accumulateFeatsGeneric(metaData,featsInDname,featsStrIn,iScale,fGetFeats);
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
function out = accumulateFeatsGeneric(metaData,featsInDname,featsStrIn,iScale,fGetFeats)

%% init
out.allCells.fov.accFeats = {};
out.allCells.strs = {};
out.allCells.date = {};
out.allCells.cellType = {};
out.allCells.source = {};
out.allCells.metEff = {};

%% Accumulation
for iexp = 1 : metaData.experiments.N
    curFname = metaData.experiments.fnames{iexp};
    % TODO: CHANGE BASED ON THE NEW FORMAT
    curDate = curFname(1:6); % 7 instead of 6!
    for in = 1 : 2
        %% count the next open spot for each condition
        curAll = length(out.allCells.fov.accFeats) + 1;        
               
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
            
            featsInFname = [featsInDname num2str(iScale) filesep curFname sprintf('_s%02d_%s.mat',itask,featsStrIn)];
            if ~exist(featsInFname,'file')
                error('file %s does not exist',featsInFname);
            end
            
            cellData = getCellData(featsInFname);
            curFeats = getLchFeats(fGetFeats,cellData);
            
            
            %% Actual accumulation
            
            if length(out.allCells.fov.accFeats) < curAll
                out.allCells.fov.accFeats{curAll} = [];
                % loaction
                out.allCells.fov.locations{curAll}.locationDeltaLbp = {};
                out.allCells.fov.locations{curAll}.locationStr = {};
            end
            
            % loaction
            iLocation = length(out.allCells.fov.locations{curAll}.locationDeltaLbp) + 1;
            out.allCells.fov.locations{curAll}.locationDeltaLbp{iLocation} = curFeats;
            out.allCells.fov.locations{curAll}.locationStr{iLocation} = sprintf('%02d',itask);
            
            % accumulation
            out.allCells.fov.accFeats{curAll} = [out.allCells.fov.accFeats{curAll},curFeats];
            out.allCells.strs{curAll} = [curCellType '_' curDate]; % just repeats for every time
            out.allCells.date{curAll} = curDate;
            out.allCells.cellType{curAll} = curCellType;
            out.allCells.cellTypeInd{curAll} = cellTypeInd; % index in metaData.cellTypes.ids
            assert(strcmp(curCellType,metaData.cellTypes.ids{cellTypeInd}));
            out.allCells.source{curAll} = curSource;
            out.allCells.metEff{curAll} = curMetEff;
        end % task
    end % first or second cell type in well
end % number of experiments

end
