function []=deleteUnusedFolders()
A = dir('wellshIIA_02*');
sizeA = size(A,1);

% these are the Folders to be deleted:
F{1}     = '3BeadsCollapsed';
F{end+1} = '3CellsCollapsed';
F{end+1} = '3EcadCollapsed';
F{end+1} = '3MyosCollapsed';
F{end+1} = '4BeadsCropSDC';

for i=1:sizeA
    folderName = A(i,1).name;
    
    cd([folderName,filesep,'data']);
    
    for k=1:length(F)
        if isdir(F{k})
            rmdir(F{k},'s')
        end
    end
    cd ../..
end