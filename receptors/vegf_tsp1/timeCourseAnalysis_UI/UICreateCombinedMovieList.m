function [] = UICreateCombinedMovieList()
%Creates CombinedMovieList object by prompting user to select MovieLists.
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%The ML must have TimePoints process.
%
%SYNOPSIS [] = UICreateCombinedMovieList()
%
%Tae H Kim, July 2015

%% Prompt user
%prompt user to determine where to save CML
[CML_FileName, CML_FilePath] = uiputfile('*.mat', 'Find a place to save your combined movie list');
%prompt user to enter the name
% nameCML is obsolete because each channel may be a different condition
% nameCML = inputdlg('Enter a name or the experimental conditions of combined movie list', 'Enter a name', 1, {'no VEGF'});
% nameCML = nameCML{1};
%prompts user to select ML files until the user presses cancel
ML_FullPath = {};
[fileName, filePath] = uigetfile('*.mat', 'Select MovieLists', 'MultiSelect', 'on');
while ~isnumeric(fileName)
    if iscell(fileName)
        ML_FullPath = [ML_FullPath cellfun(@(x) [filePath x], fileName, 'UniformOutput', false)]; %#ok<*AGROW>
    else
        ML_FullPath{end+1} = [filePath fileName];
    end
    [fileName, filePath] = uigetfile('*.mat', 'Select MovieLists', 'MultiSelect', 'on');
end
%Parameter prompt
align_StringList = {'start', 'VEGF_added'};
align_userChoice = listdlg('PromptString','Select the align event', 'SelectionMode','single', 'ListString', align_StringList);
if align_userChoice == 1
    CML_analysisPara.alignEvent = 'start';
end
if align_userChoice == 2
    CML_analysisPara.alignEvent = 'VEGF_added';
end
%% Input check

%load MLs w/o sanity check or loading MDs
MLs = cellfun(@load, ML_FullPath);
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
CML = CombinedMovieList(ML_FullPath, CML_FilePath(1:end-1));
CML.fileName_ = CML_FileName;
CML.filePath_ = CML_FilePath;
CML.analysisPara_ = CML_analysisPara;
CML.save();
end

