function [folderExists]=createRegFolder(path_images,target_dir_Reg,target_dir_Final,T_path,reg_method)

[sortedImagesCollapsedFileList]=getFileListFromFolder(path_images);

folderExists=~isempty(sortedImagesCollapsedFileList);
if folderExists
    % Check if the dir exists:
    if isdir(target_dir_Reg)
        rmdir(target_dir_Reg,'s');
        mkdir(target_dir_Reg);
    else
        mkdir(target_dir_Reg);
    end

    % Sub-pixel registration of cell images:
    perfRegInPixStep(reg_method, sortedImagesCollapsedFileList, target_dir_Reg, T_path)

    if ~isempty(target_dir_Final)
        [sortedImagesRegFileList]=getFileListFromFolder(target_dir_Reg);

        %remove the reference frames:
        % Check if the dir exists:
        if  isdir(target_dir_Final)
            rmdir(target_dir_Final,'s');
            mkdir(target_dir_Final)
        else
            mkdir(target_dir_Final)
        end

        removeRefFrame(sortedImagesRegFileList,target_dir_Final)
        rmdir(target_dir_Reg,'s');
    end
end