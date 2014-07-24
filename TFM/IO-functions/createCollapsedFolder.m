function []=createCollapsedFolder(path_images,target_dir_WithRefFrame,target_dir_Collapsed,path_RefFrame)

[sortedImageFileList]=getFileListFromFolder(path_images);

if ~isempty(sortedImageFileList)    
    % Check if the dir exists:
    if isdir(target_dir_WithRefFrame)
        rmdir(target_dir_WithRefFrame,'s');
        mkdir(target_dir_WithRefFrame);
    else
        mkdir(target_dir_WithRefFrame);
    end

    % Insert the reference frame:
    insertRefFrame(path_RefFrame,sortedImageFileList,target_dir_WithRefFrame);

    sortedImagesWithRefFrameFileList=getFileListFromFolder(target_dir_WithRefFrame);

    % Check if the dir exists:
    if isdir(target_dir_Collapsed)
        rmdir(target_dir_Collapsed,'s');
        mkdir(target_dir_Collapsed);
    else
        mkdir(target_dir_Collapsed);
    end

    collapseFileStack(sortedImagesWithRefFrameFileList, target_dir_Collapsed);
    rmdir(target_dir_WithRefFrame,'s');
end