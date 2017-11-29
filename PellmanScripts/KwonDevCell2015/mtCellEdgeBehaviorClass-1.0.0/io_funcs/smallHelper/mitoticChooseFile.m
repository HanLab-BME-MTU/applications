function fileName=chooseFile(fileString,fileDir,chooseMode,excludeStr)
%chooseFile chooses a file in a directory containing fileString according to chooseMode
%
%SYNOPSIS fileName=chooseFile(fileString,fileDir,chooseMode,excludeStr)
%
%INPUT      fileString: string the files to be looked for contains (can be
%                   regular expression) 
%           opt fileDir: directory in which to look for files (if empty: current dir)
%           opt chooseMode: 'old' selects the oldest of the files if multiple
%                               files are found
%                       'new' selects the newest of the files if multiple
%                               files are found
%                       'GUI' launches a GUI if multiple files are found (default)
%                       'all' returns the names of all matching files as cell array
%
%           opt excludeStr: files containing this string are excluded from
%                   list (can be regular expression)
%
%OUTPUT     fileName: [] if no such file
%                     name of the file if file has been found. 
%                     if chooseMode is 'all', fileName is a n by-1 cell array of strings
%
%	                     to load the file:
%	                     if ~isempty(fileName)
%	                           data=imread([fileDir,fileName]);
%	                           OR
%	                           load([fileDir,fileName]);
%	                     end
%
%                       
%
%c: 2-03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

%test input
if ~isa(fileString,'char') | nargin==0
    error('fileString has to be string');
else
    fileString = lower(fileString);
end

if nargin==1 | isempty(fileDir)
    fileDir=pwd;
end

if exist(fileDir)~=7
    error('fileDir is no directory')
end

if nargin<3 | isempty(chooseMode)
    chooseMode='GUI';
end

if ~(strcmp(chooseMode,'new')+strcmp(chooseMode,'old')+strcmp(chooseMode,'GUI')+strcmp(chooseMode,'all'))
    error('wrong chooseMode')
end
if nargin<4|isempty(excludeStr) 
    excludeStr='';
end
excludeStr = lower(excludeStr);


%find all files in directory, create list of filenames and list of dates
dirList=dir(fileDir);

nameList = arrayfun(@(x) x.name,dirList,'uniformoutput',0);
dateList = arrayfun(@(x) x.date, dirList,'uniformoutput',0);
isDirList = arrayfun(@(x) x.isdir,dirList);

% %lookfor file. strfind only allows string-by-string comparison, hence loop
% % fileIdx = [];
% % for i = 1:length(nameList)
% %     if strfind(lower(nameList{i}),fileString)&~isdir(nameList{i})
% %         fileIdx = [fileIdx; i];
% %     end
% % end
% % fileIdx=strmatch(fileString,nameList);
% % %check if any of the fileIdx contains the excludeString
% % if ~isempty(excludeStr)
% %     for i=1:length(fileIdx)
% %         %use strfind, else files are excluded that are contained in excludeStr
% %         if strfind(lower(char(nameList(fileIdx(i)))),excludeStr)
% %             fileIdx(i)=0;
% %         end
% %     end
% %     fileIdx(find(fileIdx==0))=[];
% % end

% Improvement:
% look for file with regexp; don't call the slow isdir anymore
includeCell = regexp(lower(nameList),fileString);
excludeCell = regexp(lower(nameList),excludeStr);
fileIdx = zeros(length(nameList));
% we have to loop because isempty does not work on a cell in the way we'd
% need
for i=length(nameList):-1:1
    if ~isDirList(i) && ~isempty(includeCell{i}) && isempty(excludeCell{i})
        fileIdx(i) = 1;
    end
end
fileIdx = find(fileIdx);




%if no file returned: return empty
%if one file returned: return fname
%if more than one file returned: choose according to chooseMode
switch length(fileIdx)
    case 0 %no file found
        fileName=[];
    case 1 %exactly one file found
        fileName=char(nameList(fileIdx));
    otherwise %more than one file found: continue according to chooseMode
        switch chooseMode
            case 'old' %select oldest file
                dateListFiles=dateList(fileIdx);
                [dummy,ranked]=sort(datenum(dateListFiles));
                %take the first file of the date-sorted indices
                fileName=char(nameList(fileIdx(ranked(1))));
            case 'new' %select newest file
                dateListFiles=dateList(fileIdx);
                [dummy,ranked]=sort(datenum(dateListFiles));
                %take the last file of the date-sorted indices
                fileName=char(nameList(fileIdx(ranked(end))));
            case 'GUI' %launch GUI
                listData=nameList(fileIdx);
                chooseNum=chooseFileGUI(listData);
                fileName=char(listData(chooseNum));
            case 'all' %keep all files found
                fileName = nameList(fileIdx);
        end
end
