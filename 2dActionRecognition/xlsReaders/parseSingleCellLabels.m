% update single cell DB
function [] = parseSingleCellLabels(xlsFname,dbFname,metaDataFname,analysisDname)

%% File names
if nargin == 0
    xlsFname = 'SingleCells20160114_ESW_AZ.xlsx';
end

if nargin < 2
    dbFname = 'singleCellLabelDB.mat';
end

if nargin < 3
    metaDataFname = 'Experiments20151023.mat';
end

if nargin < 4
    analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
end

xlsFname = [analysisDname 'SingleCellAnalysis/' xlsFname];
dbFname = [analysisDname 'SingleCellAnalysis/' dbFname];
metaDataFname = [analysisDname 'MetaData/' metaDataFname];
dataDname = [analysisDname 'Data/'];

%% read xls and merge with database
load(metaDataFname); % metaData
[cellInfo,cellTXY] = readCellInfo(xlsFname,metaData,dataDname);

% validate(cellInfo); % TODO: IMPLEMENT (consistency with file names, cell types etc.) - need another parameter pointing to the metaData files (experiments and cell types in the MetaData folder)

if exist(dbFname,'file')
    load(dbFname); % cellLabelDB
else
    cellInfoDB = {};
    TYXDB = [];
end

[cellInfoDB,TYXDB,cellGroups] = merge(cellInfo,cellTXY,cellInfoDB,TYXDB);

assert(~hasDuplicates(TYXDB));

save(dbFname,'cellInfoDB','TYXDB','cellGroups');

end

%%
function [cellInfo,cellTXY] = readCellInfo(xlsFname,metaData,dataDname)
[numbers strings misc] = xlsread(xlsFname);
firstRow = 2;
lastRow = size(strings,1);

[numsAll,txtAll,rawAll] = xlsread(xlsFname,sprintf('A%d:E%d',firstRow,lastRow));

% n = 0;
n = lastRow - firstRow + 1;
cellInfo = cell(1,n);
cellTXY = nan(n,3);
for line = firstRow : lastRow
    i = line - firstRow + 1;
    
    curVideoName = txtAll{i,1};
    curCellType = txtAll{i,2};
    curLabels = strsplit(txtAll{i,4},', ');
    curCellID = numsAll(i);
    curDescr = txtAll{i,5};
    
    curCellSource = getCellSource(curCellType,metaData.cellTypes);   
    
    curDname = [dataDname curCellSource filesep curVideoName(1:end-4) filesep curVideoName filesep];
    
    if ~exist(curDname,'dir')
        error(['Directory '  curDname ' does not exists']);
    end
    
    load([curDname 'tracking/cellIdTYX.mat']);                
    
    cellInfo{i}.fname = curVideoName;
    cellInfo{i}.date = curVideoName(1:6);
    cellInfo{i}.ts = cellTYX{curCellID}.ts;
    cellInfo{i}.ys = cellTYX{curCellID}.ys;
    cellInfo{i}.xs = cellTYX{curCellID}.xs;
    cellInfo{i}.excelID = curCellID;
    cellInfo{i}.cellType = curCellType;
    cellInfo{i}.source = curCellSource;
    cellInfo{i}.dir =  curDname;
    cellInfo{i}.labels = curLabels;
    cellInfo{i}.desc = curDescr;
    cellTXY(i,:) = [cellInfo{i}.ts(1),cellInfo{i}.ys(1),cellInfo{i}.xs(1)];        
end

end


%% groupsByTreatments - cell array tha holds pairs <treatmentStr,[list of indices]>
function [groupsByTreatments] = clusterByTreatments(curCellInfo)

ncells = length(curCellInfo);

% accumulate labels and create cell array of unique labels
accLables = {};
nlabels = 0;
for icell = 1 : ncells
    curLabels = curCellInfo{icell}.labels;
    for ilabel = 1 : length(curLabels)
        nlabels = nlabels + 1;
        accLables{nlabels} =  curLabels{ilabel};
    end
end

uniqueLabels = unique(accLables);
nUniqueLabels = length(uniqueLabels);

% initiate groups
groupsByTreatments = cell(1,nUniqueLabels);
for ilabel = 1 : nUniqueLabels
    groupsByTreatments{ilabel}.str = uniqueLabels{ilabel};
    groupsByTreatments{ilabel}.inds = [];
end

% update based on labels
for icell = 1 : ncells
    curCellLabels = curCellInfo{icell}.labels;
    for ilabel = 1 : length(curCellLabels)
        curLabel = curCellLabels{ilabel};
        indLabel = find(strcmp(curLabel,uniqueLabels));
        groupsByTreatments{indLabel}.inds = [groupsByTreatments{indLabel}.inds, icell];
    end
end
end
    
function [groupsByDates] = clusterByDates(curCellInfo)

ncells = length(curCellInfo);

% accumulate labels and create cell array of unique labels
accDates = cell(1,ncells);
for icell = 1 : ncells
    curDate = curCellInfo{icell}.date;
    accDates{icell} = curDate;
end

uniqueDates = unique(accDates);
nUniqueDates = length(uniqueDates);

% initiate groups
groupsByDates = cell(1,nUniqueDates);
for idate = 1 : nUniqueDates
    groupsByDates{idate}.str = uniqueDates{idate};
    groupsByDates{idate}.inds = [];
end

% update based on dates
for icell = 1 : ncells
    curCellDate = curCellInfo{icell}.date;    
    indDate = find(strcmp(curCellDate,uniqueDates));
    groupsByDates{indDate}.inds = [groupsByDates{indDate}.inds, icell];    
end
end

function [groupsByCellTypes] = clusterByCellTypes(curCellInfo)

ncells = length(curCellInfo);

% accumulate labels and create cell array of unique labels
accCellTypes = cell(1,ncells);
for icell = 1 : ncells
    curCellType = curCellInfo{icell}.cellType;
    accCellTypes{icell} = curCellType;
end

uniqueCellTypes = unique(accCellTypes);
nUniqueCellTypes = length(uniqueCellTypes);

% initiate groups
groupsByCellTypes = cell(1,nUniqueCellTypes);
for icelltype = 1 : nUniqueCellTypes
    groupsByCellTypes{icelltype}.str = uniqueCellTypes{icelltype};
    groupsByCellTypes{icelltype}.inds = [];
end

% update based on dates
for icell = 1 : ncells
    curCellCellType = curCellInfo{icell}.cellType;    
    indCellType = find(strcmp(curCellCellType,uniqueCellTypes));
    groupsByCellTypes{indCellType}.inds = [groupsByCellTypes{indCellType}.inds, icell];    
end
end





%%
function [cellInfoDB,TYX,cellGroups] = merge(curCellLabels,cellTXY,cellInfoDB,TYX)

for i = (length(cellInfoDB) + 1) : (length(cellInfoDB) + length(curCellLabels))
    curCellLabels{i}.serialNum = i;
end

cellInfoDB = [cellInfoDB,curCellLabels];
TYX = [TYX,cellTXY];

cellGroups.byTreatments = clusterByTreatments(cellInfoDB);
cellGroups.byDates = clusterByDates(cellInfoDB);
cellGroups.byCellTypes = clusterByCellTypes(cellInfoDB);
end

%%
function [hasDuplicates] = hasDuplicates(TYX)
fprintf('hasDuplicates not implemented\n');
hasDuplicates = false;
return;
D = pdist2(TYX);


% nCur = length(curCellLabels);
% nDB = length(cellInfoDB);
% nReplicateIgnored = 0;
% 
% for icur = 1 : nCur
%     newCellInfo = curCellLabels{icur};
% for idb = 1 : nDB
%     dbCellInfo = cellInfoDB{idb};
%     
% end
% end
end
