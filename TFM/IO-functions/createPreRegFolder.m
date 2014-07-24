function []=createPreRegFolder(path_images,path_tmp,path_PreRegT)
[sortedImageFileList]=getFileListFromFolder(path_images);

if ~isempty(sortedImageFileList)
    % Check if the dir exists:
    if isdir(path_tmp)
        rmdir(path_tmp,'s');
        mkdir(path_tmp);
    else
        mkdir(path_tmp);
    end
    
    perfRegInPixStep('pixelwise', sortedImageFileList, path_tmp, path_PreRegT);
    delete([path_images,filesep,'*']);
    movefile([path_tmp,filesep,'*'],path_images);
    rmdir(path_tmp,'s');
end





    
