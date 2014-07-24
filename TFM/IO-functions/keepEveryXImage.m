function []=keepEveryXImage(fileList,x)

numFiles=length(fileList);
keepList=1:x:numFiles;
badList =setdiff(1:numFiles,keepList);
for iframe=badList
    delete(fileList{iframe});
end

for k=1:length(keepList)
    keepFileList{k}=fileList{keepList(k)};
end
collapseFileStack(keepFileList,-1);
