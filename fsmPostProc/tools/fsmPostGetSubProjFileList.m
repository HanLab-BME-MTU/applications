function [fileList,success]=fsmPostGetSubProjFileList(fromDir,root)

% INPUT
% fromDir       : directory to scan
% root          : string containing the root of the files to look for (e.g. 'cands' for all cands###.mat files)
%
% OUTPUT
% fileList      : file list

% Set flag for success
success=0;

% Initialize output
fileList=[];

% Check that the passed directory exists and is a directory
if ~isdir(fromDir)
    disp('The passed directory does not exist.');
    return
end

% Read the content of the directory
everything=dir(fromDir);

% Discard what is not a file
everything(find([everything.isdir]==1))=[];

if isempty(everything)
    disp('The directory is empty.');
    return
end

% Create empty fileList
fileList=cell(1,length(everything));

% Fill it
c=0;
for i=1:length(everything)
    if strncmp(everything(i).name,root,length(root))==1
        c=c+1;
        fileList(1,c)={[fromDir,filesep,everything(i).name]};
    end
end

% Make sure it is sorted
fileList=sort(fileList);

% Return success
success=1;