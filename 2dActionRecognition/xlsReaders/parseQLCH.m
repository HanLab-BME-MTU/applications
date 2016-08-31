% experimentsFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/Experiments20150616_AN.xlsx';
% cellTypeFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/CellType20150616_AN.xlsx';
% experimentsFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/Experiments20150521.xlsx';
% cellTypeFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/CellType20150520.xlsx';
% experimentsFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/Experiments20150807.xlsx';
% cellTypeFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/CellType20150806.xlsx';
% experimentsFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/Experiments20150901.xlsx';
% cellTypeFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/CellType20150901.xlsx';
% experimentsFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/Experiments20151023.xlsx';
% cellTypeFname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/MetaData/CellType20151022.xlsx';


% parseQLCH(experimentsFname,cellTypeFname)
function [] = parseQLCH(experimentsFname,cellTypeFname)

if nargin < 2
    experimentsFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/Experiments20151023.xlsx';
    cellTypeFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/CellType20160111.xlsx';
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

exps.cellType1 = cell(1,exps.N);
exps.source1 = cell(1,exps.N); % CellLines, Melanocytes, Tumors
exps.n1 = cell(1,exps.N);
exps.day1 = cell(1,exps.N);
exps.str1 = cell(1,exps.N);

exps.cellType2 = cell(1,exps.N);
exps.source2 = cell(1,exps.N);
exps.n1 = cell(1,exps.N);
exps.day2 = cell(1,exps.N);
exps.str2 = cell(1,exps.N);

exps.exclude = cell(1,exps.N);


[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:F%d',firstRow,lastRow));  

for line = firstRow : lastRow 
    i = line - firstRow + 1;
    fname = cell2mat(rawAll(i,1));
    
    curExp = strsplit(fname,'_');
    
    assert(length(curExp{1}) == 6 && (length(curExp) == 2 || length(curExp) == 3));
    
    date = curExp{1};
    cellType1 = lower(curExp{2});
        
    exps.fnames{i} = fname;    
    exps.date{i} = date; 
    
    exps.cellType1{i} = cellType1;
    exps.source1{i} = cellTypes.source{ismember(lower(cellTypes.ids),cellType1)};
    exps.n1{i} = cell2mat(rawAll(i,2));
    exps.day1{i} = cell2mat(rawAll(i,3));
    exps.str1{i} = [date '_' cellType1];
    
    if length(curExp) == 3
        cellType2 = lower(curExp{3});
        
        exps.cellType2{i} = cellType2;
        exps.source2{i} = cellTypes.source{ismember(lower(cellTypes.ids),cellType2)};
        exps.n2{i} = cell2mat(rawAll(i,4));        
        exps.day2{i} = cell2mat(rawAll(i,5));
        exps.str2{i} = [date '_' cellType2];
    else
        exps.n2{i} = 0;
    end
    
    excludeListStr = cell2mat(rawAll(i,6));
    exps.exclude{i} = getExcludeList(excludeListStr);    
end
end

%% 
function [tasks] = getTasks(exps)

% each task is a tuple (experiment index, task index)
tasks.exps = [];
tasks.tasks = [];
tasks.N = 0;

for i = 1 : exps.N
    curExpN = exps.n1{i} + exps.n2{i};
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

