function []=deleteBadImages(inputFileList)
%remove bad pictures
percentile=0.50;


if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image of the Stack to be cleaned');
   
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

n = numel(inputFileList);
totInt=zeros(min(100,n),1);

for i = 1:min(100,n)
    % Read image.
    I = double(imread(inputFileList{i}));
    totInt(i)=sum(sum(I));
end

refValue=mean(totInt);

totInt=zeros(n,1);
indBadImages=[];
for i = 1:n
    I = double(imread(inputFileList{i}));
    totInt(i)=sum(sum(I));
    if abs(totInt(i)-refValue)/refValue>percentile
        indBadImages(end+1)=i;
        %imagesc(I)
    end
end
display('The following images will be removed from the stack:');

indBadImages

stepSize=input('How many frames do you want to analyze: 1 out of: ');

targetDir = uigetdir('','Please select target directory');

listToAnalyze=[];
for i = 1:stepSize:n
    while ismember(i,indBadImages)
        i=i+1;
    end
    sourceFileName=inputFileList{i};
    targetFileName=[targetDir filesep sourceFileName(length(pathname)+1:end)];
    copyfile(sourceFileName,targetFileName);
    listToAnalyze(end+1)=i;
end
display('saved the file-stack index, not the actual image number, this should be changed')
save([targetDir filesep 'indexOfDeletedImages.mat'],'indBadImages')

moreStacks= input('are there more stacks to treated in the same way? Yes=1, No=0: ');

if moreStacks==1
    if nargin < 1 || isempty(inputFileList)
       [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
           'Select First Image of the Stack');

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
    
    targetDir = uigetdir('','Please select target directory');

    for i = listToAnalyze
        sourceFileName=inputFileList{i};
        targetFileName=[targetDir filesep sourceFileName(length(pathname)+1:end)];
        copyfile(sourceFileName,targetFileName);        
    end
    save([targetDir filesep 'indexOfDeletedImages.mat'],'indBadImages')
end