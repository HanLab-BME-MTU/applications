function [param] = UICreateCombinedMovieList(varargin)
%Creates CombinedMovieList object by prompting user to select MovieLists.
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%The ML must have TimePoints process.
%
%SYNOPSIS [] = UICreateCombinedMovieList()
%
%Tae H Kim, July 2015

ip = inputParser;
ip.StructExpand = true;
ip.addParameter('CML_FileName',[]);
ip.addParameter('CML_FilePath',[]);
ip.addParameter('nameCML',[]);
ip.addParameter('ML_FullPath',[]);
ip.addParameter('fileName',[]);
ip.addParameter('filePath',[]);
ip.addParameter('align_userChoice',[]);
ip.addParameter('CML_analysisPara',[]);
ip.parse(varargin{:});
param = ip.Results;


%% Prompt user
%prompt user to determine where to save CML
if(isempty(param.CML_FileName) || isempty(param.CML_FilePath))
    [param.CML_FileName, param.CML_FilePath] = uiputfile('*.mat', 'Find a place to save your combined movie list');
end
%prompt user to enter the name
if(isempty(param.nameCML))
    param.nameCML = inputdlg('Enter a name or the experimental conditions of combined movie list', 'Enter a name', 1, {'no VEGF'});
    param.nameCML = param.nameCML{1};
end
%prompts user to select ML files until the user presses cancel
if(isempty(param.ML_FullPath))
    param.ML_FullPath = {};
    [param.fileName, param.filePath] = uigetfile('*.mat', 'Select MovieLists', 'MultiSelect', 'on');
    while ~isnumeric(param.fileName)
        if iscell(param.fileName)
            param.ML_FullPath = [param.ML_FullPath cellfun(@(x) [param.filePath x], param.fileName, 'UniformOutput', false)]; %#ok<*AGROW>
        else
            param.ML_FullPath{end+1} = [param.filePath param.fileName];
        end
        [param.fileName, param.filePath] = uigetfile('*.mat', 'Select MovieLists', 'MultiSelect', 'on');
    end
end
%Parameter prompt
if(isempty(param.CML_analysisPara))
    align_StringList = {'start', 'VEGF_added'};
    param.align_userChoice = listdlg('PromptString','Select the align event', 'SelectionMode','single', 'ListString', align_StringList);
    if param.align_userChoice == 1
        param.CML_analysisPara.alignEvent = 'start';
    end
    if param.align_userChoice == 2
        param.CML_analysisPara.alignEvent = 'VEGF_added';
    end
end
%% Input check

%load MLs w/o sanity check or loading MDs
MLs = cellfun(@load, param.ML_FullPath);
%check if all are MovieList object
try
    assert(all(arrayfun(@(x) isa(x.ML, 'MovieList'), MLs)));
    %2 ways error could be thrown:
    %1: loaded ML is not a MovieList
    %2: loaded file doesn't have variable 'ML'
    %in second case, error will be thrown before assert
catch
    error('Selected files do not all contain MovieList object.');
end
%check if all MovieList object has TimePoints process
assert(~any(arrayfun(@(x) isempty(x.ML.getProcessIndex('TimePoints')), MLs)), 'Not all Movie Lists contain TimePoints process. Makes sure MovieLists have timePoints process.');

%% Create CombinedMovieList object and save
CML = CombinedMovieList(param.ML_FullPath, param.CML_FilePath(1:end-1), 'name', param.nameCML);
CML.fileName_ = param.CML_FileName;
CML.filePath_ = param.CML_FilePath;
CML.analysisPara_ = param.CML_analysisPara;
CML.save();
end

