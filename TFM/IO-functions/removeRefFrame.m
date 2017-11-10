function []=removeRefFrame(inputFileList,target_dir)
%this program finds the the ref. frames and deletes them:

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


I1 = double(imread(inputFileList{1}));
I2 = double(imread(inputFileList{2}));
if length(inputFileList)>3
    I3 = double(imread(inputFileList{3}));
    I4 = double(imread(inputFileList{4}));
else
    I3=I1;
    I4=I1;
end

if sum(sum(I3-I1))==0
    fprintf('\nIt is assumed, that images no: 1,3,5... are the reference frames\n');
    toSave=2:2:n;
elseif sum(sum(I4-I2))==0
    fprintf('\nIt is assumed, that images no: 2,4,6... are the reference frames\n');
    toSave=1:2:n;
else
    error('Cannot find reference frames!');
end

padZeros=floor(log10(n))+1;

for i=toSave
    %create the new filenames for the data:
    [~, body, no, ext]=getFilenameBody(inputFileList{i});
    if mod(i,2)==0
        newNo=round(str2double(no))/2;
    else
        newNo=(round(str2double(no))+1)/2;
    end
    sourceFileName=inputFileList{i};
    
    targetFileName=[target_dir filesep body num2str(newNo,['%0.',int2str(padZeros),'d']) ext];
    
    copyfile(sourceFileName,targetFileName);
end