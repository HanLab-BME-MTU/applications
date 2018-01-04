%% LCH_orginizeFeatures2018
% Orginize in conivnient data structures + normalize
function [] = LCH_orginizeFeatures2018(featsDname,featsStrID,always)
featsFname = [featsDname filesep featsStrID '_all.mat'];
load(featsFname); % allCells

outdir = [featsDname filesep featsStrID];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

outFname = [outdir filesep featsStrID '.mat'];

if exist(outFname,'file') && ~always
    fprintf(sprintf('%s pre-exists',outFname));
    return;
end


analysisDname = [featsDname filesep '..' filesep '..' filesep 'All' filesep];

% Number of cell type - date combinations
nCellTypeDate = length(allCells.strs);

allCellsMovieData = {};
ncells = 0;
for iCellTypeDate = 1 : nCellTypeDate 
    expStr = allCells.expStr{iCellTypeDate};
    expDate = allCells.date{iCellTypeDate};
    cellType = allCells.cellType{iCellTypeDate};
    metEff = allCells.metEff{iCellTypeDate};        
    
    locationData = allCells.locations{iCellTypeDate};
    nLocations = length(locationData.locationTXY);
    
    for ilocation = 1 : nLocations
        curLocationStr = locationData.locationStr{ilocation};
        curLocationFeats = locationData.locationFeats{ilocation};
        curLocationTXY = locationData.locationTXY{ilocation};        
                
        %% 
        nCellsLocation = length(curLocationTXY);
        
        for icell = 1 : nCellsLocation
            ncells = ncells + 1;
            curCellMovieData.expStr = expStr;
            curCellMovieData.date = expDate;
            curCellMovieData.cellType = cellType;
            curCellMovieData.metEff = metEff;
            
            curCellMovieData.ts = curLocationTXY{icell}.ts;
            curCellMovieData.xs = curLocationTXY{icell}.xs;
            curCellMovieData.ys = curLocationTXY{icell}.ys;
            
            curCellMovieData.locationStr = curLocationStr;                        
            
            % Unique cell ID (assumes the trajecotry length is constant)
            curCellMovieData.key = [...
                date '_' cellType ...
                '_s' curCellMovieData.locationStr ...
                '_t' num2str(curCellMovieData.ts(1)) ...
                '_x' num2str(round(curCellMovieData.xs(1))) ...
                '_y' num2str(round(curCellMovieData.ys(1)))];
            
            mdFname = [analysisDname curCellMovieData.expStr filesep...
                curCellMovieData.expStr '_s' curCellMovieData.locationStr filesep...
                curCellMovieData.expStr '_s' curCellMovieData.locationStr '.mat'];                        
            
            % For locations where the location do not have a trailing zero
            if ~exist(mdFname,'file')
                mdFname = [analysisDname curCellMovieData.expStr filesep...
                    curCellMovieData.expStr '_s' sprintf('%d',curCellMovieData.locationStr) filesep...
                    curCellMovieData.expStr '_s' sprintf('%d',curCellMovieData.locationStr) '.mat'];
            end
            
            assert(logical(exist(mdFname,'file')));
            
            % Movie data to access the cell movies
            curCellMovieData.MD = mdFname;
            
            % Empty annotation map
            curCellMovieData.annotations = containers.Map(...
                {'blebbing', 'moving', 'short extensions',...
                'dead','debris','out of focus','plastic',...
                'stretched','collagen deformation', 'long extensions',...
                'multiple cells', 'big', 'small', 'medium'},...
                {nan,nan,nan,nan,nan,...
                nan,nan,nan,nan,nan,...
                nan,nan,nan,nan});
            
            % features
            curFeats = curLocationFeats(:,icell);
            %             curCellMovieData.featureMap = containers.Map({featsStrID},{curFeats});
            curCellMovieData.feats = curFeats;
            
            % include cell
            allCellsMovieData{ncells} = curCellMovieData;
        end
    end    
end

%% To list
[nCells,cellFeats,cellTypes,metEffs,dates,sources] = cellMD2list(allCellsMovieData);

%% To matrix
[featsAll,indsAll] = LCH_feats2Matrix(cellFeats);

%% Normalize
[meanValsAll,stdValsAll] = LCH_getNormFeats(featsAll);
[featsNormAll] = LCH_setNormFeats(featsAll,meanValsAll,stdValsAll);

%% TODO:
% 1. PCA
% 2. Partition to: all, tumors, cell lines, melanocytes?

%% Save
save(outFname,...
    'allCellsMovieData',...
    'nCells','cellFeats','cellTypes','metEffs','dates','sources',...
    'featsAll','indsAll','meanValsAll','stdValsAll','featsNormAll'...
    );

end


%%
function [nCells,cellFeats,cellTypes,metEffs,dates,sources] = cellMD2list(allCellsMovieData)
nCells = length(allCellsMovieData);
cellFeats = cell(1,nCells);
cellTypes = cell(1,nCells);
metEffs = nan(1,nCells);
dates = cell(1,nCells);
sources = cell(1,nCells);

for icell = 1 : nCells
    curCell = allCellsMovieData{icell};
    cellFeats{icell} = curCell.feats;
    cellTypes{icell} = curCell.cellType;
    metEffs(icell) = curCell.metEff;
    dates{icell} = curCell.date;
    sources{icell} = LCH_getSource(cellTypes{icell});
end
end
