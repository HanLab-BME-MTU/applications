function [anaDat,pathName,dataProperties]=calcAnaDatFromProjectData(pathName,fileName)
%function that calculates anaDat from projectData fileName located in pathName

if nargin == 0 | isempty(pathName) 
    %check if default biodata-dir exists and cd if exist
    mainDir = cdBiodata(2);
    
    %get project data file
    [fileName,pathName] = uigetfile({'*-data*-??-???-????-??-??-??.mat','project data files'},'select project data file');
    
    if fileName == 0
        error('No file loaded')
    end
    
end

if isdir(pathName)
    cd(pathName);
else
    error(['directory ''',pathName,''' does not exist']);
end

%load dataProperties
load(fileName,'dataProperties');
if ~exist('dataProperties','var')
    error('No dataProperties in project data: corrupt data file');
end

%get lastResult, load corresponding idlist (if possible)
load(fileName,'lastResult');
if isempty(findstr('idlist',lastResult))
    error('No idlist of any kind in project data: run spotID first');
end

%load either idlist, idlist_L, idlisttrack, or idlisttrack_L
load(fileName,lastResult);

%name loaded data idlist
eval(['idlist = ',lastResult,';']);


%-------------calculate anaDat
anaDat = adgui_calc_anaDat(idlist,dataProperties,lastResult);
%-----------------------------


