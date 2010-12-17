function [targetFileName]=collapseFileStack(inputFileList, target_dir)
% this program collapses the files stack, such that the first frame starts
% with one and the last frame ends with #frames.

%read in Stack of bead images:
if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.TIF';'*.tif';'*.jpg';'*.png';'*.*'}, ...
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

doOverwrite=0;
%get the target directory:
if nargin < 2 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
elseif nargin==2 && isnumeric(target_dir) && target_dir==-1
    % Then overwrite the old images.
    doOverwrite=1;
end

n = numel(inputFileList);

padZeros=floor(log10(n))+1;

mapfnos=zeros(length(inputFileList),2);
for i=1:n
    %create the new filenames for the data:
    [path, body, fno , ext]=getFilenameBody(inputFileList{i});
    sourceFileName=inputFileList{i};
    if doOverwrite==0
        targetFileName{i}=[target_dir filesep body num2str(i,['%0.',int2str(padZeros),'d']) ext];
        copyfile(sourceFileName,targetFileName{i});
    else
        % Old data will be overwritten. Append a new number to the complete old filename:
        targetFileName{i}=[path,filesep, body,'_newfno_',num2str(i,['%0.',int2str(padZeros),'d']) ext];
        movefile(sourceFileName,targetFileName{i},'f');
        mapfnos(i,:)= [str2num(fno), i];
    end
end
if doOverwrite==1
   % save the map-file
   fid = fopen([path,filesep,'mapFrameNos-',date,'.txt'], 'w');
   fprintf(fid, '%d = %d\n', mapfnos');
   fclose(fid);
end