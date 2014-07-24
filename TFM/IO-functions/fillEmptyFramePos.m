function []=fillEmptyFramePos(inputFileList,target_dir)
% this program fills in missing frames by multiple copies of the previous
% frame in the filestack.

nargin=0;

%read in Stack of bead images:
if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image of the Stack to be analyzed');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

%get the target directory:
if nargin < 2 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

n = numel(inputFileList);
padZeros=floor(log10(2*n))+1;

for i=1:n-1
    sourceFileName=inputFileList{i};    
    [path body no_curr, ext]=getFilenameBody(inputFileList{i});
    [ ~,    ~, no_next,  ~ ]=getFilenameBody(inputFileList{i+1});
    
    %first copy the original image to the new position:
    targetFileName=[target_dir filesep body num2str(str2double(no_curr),['%0.',int2str(padZeros),'d']) ext];
    copyfile(sourceFileName,targetFileName);
    
    % In case there is a missing slot, copy the current image to this
    % slot
    no_diff=str2double(no_next)-str2double(no_curr);
    if no_diff>1
        for k=1:no_diff-1;
            no_fill=str2double(no_curr)+k;
            targetFileName=[target_dir filesep body num2str(no_fill,['%0.',int2str(padZeros),'d']) ext];
            copyfile(sourceFileName,targetFileName);
        end
    end
end

%copy the last image:
sourceFileName=inputFileList{n};    
[path body no_curr, ext]=getFilenameBody(inputFileList{n});
targetFileName=[target_dir filesep body num2str(str2double(no_curr),['%0.',int2str(padZeros),'d']) ext];
copyfile(sourceFileName,targetFileName);

