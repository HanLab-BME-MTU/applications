function []=insertRefFrame(ref_filename_path,inputFileList,target_dir)
% this program re-orders the bead images such that:
% 1. newImage_1  = reference frame
% 2. newImage_2  = beadImage_1
% 3. newImage_3  = reference frame
% 4. newImage_4  = beadImage_2
% 5. newImage_5  = reference frame
% 6. newImage_20 = beadImage_10
% 7. newImage_21 = reference frame

%read in ref-frame:
if nargin <1 || isempty(ref_filename_path)
    [ref_filename, ref_pathname] = uigetfile({'*.TIF';'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select the Reference Frame');
   ref_filename_path=[ ref_pathname filesep ref_filename];
end

%read in Stack of bead images:
if nargin < 2 || isempty(inputFileList)
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

%get the target directory:
if nargin < 3 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

n = numel(inputFileList);
padZeros=floor(log10(2*n))+1;

for i=1:n
    %create the new filenames for the data:
    [path body no ext]=getFilenameBody(inputFileList{i});
    newNo=2*round(str2double(no));
    sourceFileName=inputFileList{i};
    targetFileName=[target_dir filesep body num2str(newNo,['%0.',int2str(padZeros),'d']) ext];
    copyfile(sourceFileName,targetFileName);
    
    %create the filenames for the reference frame:
    targetFileName=[target_dir filesep body num2str(newNo-1,['%0.',int2str(padZeros),'d']) ext];    
    copyfile(ref_filename_path,targetFileName);
end