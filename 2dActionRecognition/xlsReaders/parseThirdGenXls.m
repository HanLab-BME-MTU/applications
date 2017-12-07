
function [] = parseThirdGenXls(experimentsFname,cellTypeFname)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/xlsReaders'));

if nargin < 2
    %     experimentsFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/ThirdGen20170118_noControl.xlsx';
    %     cellTypeFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/CellTypeThirdGen20170118.xlsx';
    experimentsFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/ExperimentsAll20170509.xlsx';
    cellTypeFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/CellTypeAll20170509.xlsx';
end

[pathstr,name,ext] = fileparts(experimentsFname);

metaData.cellTypes = getCellTypes(cellTypeFname); % cellTypes links cell type to source and other information
metaData.experiments = getExperiments(experimentsFname,metaData.cellTypes); % experiments, N1/N2 for 2 cell types, exclude list
metaData.tasks = getTasks(metaData.experiments); % tasks include list of tuples (exp,task) to be dispatched to the cluster

outFname = [pathstr filesep name '.mat'];

save(outFname,'metaData');
end

%%
function [cellTypes] = getCellTypes(fname)
firstRow = 2;

[numbers strings misc] = xlsread(fname);
lastRow = size(strings,1);
N = lastRow - firstRow + 1;

cellTypes.ids = cell(1,N);
cellTypes.source = cell(1,N);
cellTypes.metastaticEfficiency = -1 * ones(1,N);


[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:C%d',firstRow,lastRow));  

for line = firstRow : lastRow 
    i = line - firstRow + 1;
    cellTypes.ids{i} = rawAll{i,1};
    cellTypes.source{i} = rawAll{i,2};
    cellTypes.metastaticEfficiency(i) = cell2mat(rawAll(i,3));    
end
end

%%
function [exps] = getExperiments(fname,cellTypes)
firstRow = 2;

[numbers strings misc] = xlsread(fname);
lastRow = size(strings,1);
exps.N = lastRow - firstRow + 1;

exps.fnames = cell(1,exps.N);

exps.cellTypes = cell(1,exps.N);
exps.sources = cell(1,exps.N);
exps.ns = cell(1,exps.N);
exps.strs = cell(1,exps.N);

exps.exclude = cell(1,exps.N);
exps.experimentalist = cell(1,exps.N);


[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:I%d',firstRow,lastRow));  

for line = firstRow : lastRow 
    i = line - firstRow + 1;
    fname = cell2mat(rawAll(i,1));
    
    curExp = strsplit(fname,'_');
    
    assert(length(curExp{1}) == 6 && length(curExp) >= 2); % we could add assertion to the cell type names
    
    date = curExp{1};

    curTask = 0;
    cellTypesCur = {};
    cellSourcesCur = {};
    strs = {};
    for iCellType = 2 : length(curExp)
        curCellType = lower(curExp{iCellType});
        curCellSource = cellTypes.source{ismember(lower(cellTypes.ids),curCellType)};
        curN = numsAll(i,iCellType-1);

        % do one by one...
        for icount = 1 : curN
            curTask = curTask + 1;
            cellTypesCur = [cellTypesCur curCellType];
            cellSourcesCur = [cellSourcesCur curCellSource];
            strs = [strs [date '_' curCellType]];
        end
    end
    
    exps.fnames{i} = fname;
    exps.date{i} = date;
    exps.ns{i} = curTask;
    exps.cellTypes{i} = cellTypesCur;
    exps.source{i} = cellSourcesCur;
    exps.strs{i} = strs;   
    
    excludeListStr = cell2mat(rawAll(i,6));
    exps.exclude{i} = getExcludeList(excludeListStr); 
   
    exps.experimentalist{i} = cell2mat(rawAll(i,7));
end
end

%% 
function [tasks] = getTasks(exps)

% each task is a tuple (experiment index, task index)
tasks.exps = [];
tasks.tasks = [];
tasks.N = 0;

for i = 1 : exps.N
    curExpN = exps.ns{i};
    tasksList = 1 : curExpN;
    tasksList = tasksList(~ismember(tasksList,exps.exclude{i}));
    nTasks = length(tasksList);
    
    tasks.exps = [tasks.exps, i*ones(1,nTasks)];
    tasks.tasks = [tasks.tasks, tasksList];    
end
tasks.N = length(tasks.tasks);
end

%% Get exclude list for pattern such as '1,2,3,4-8,9-11,12'
function [excludeList] = getExcludeList(str)

if isnan(str)
    excludeList = [];
    return;
end

if isnumeric(str)
    str = int2str(str);
end

splitstring = regexp(str,',','split');
n = length(splitstring);

excludeList = [];
for i = 1 : n
    curStr = splitstring{i};
    splitstring1 = regexp(curStr,'-','split');
    assert(length(splitstring1) == 1 || length(splitstring1) == 2)
    if length(splitstring1) == 1
        excludeList = [excludeList str2num(splitstring1{1})];    
    else
        from = str2num(splitstring1{1});
        to = str2num(splitstring1{2});
        assert(from < to);        
        for curTask = from : to
            excludeList = [excludeList curTask];
        end
    end
end
end

