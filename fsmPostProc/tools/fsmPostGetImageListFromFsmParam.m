function [imageFileList,imageFirstIndex,imageLastIndex,status]=fsmPostGetImageListFromFsmParam(fsmParam)

status=0;

% Check whether fsmParm has been produced by older version of FSM. In the
% new version, fsmParam.specific.fileList is a cell array. To ensure
% backward-compatibility, we update that field accordingly.

if ~iscell(fsmParam.specific.fileList)
    n = size(fsmParam.specific.fileList,1);
    
    fsmParam.specific.fileList = arrayfun(@(i) ...
        strtrim(fsmParam.specific.fileList(i,:)), 1:n,'UniformOutput',false);

    save([fsmParam.main.path,filesep,'fsmParam.mat'],'fsmParam');
end

% Extract info from fsmParam
imageFileList=fsmParam.specific.fileList;
imageFirstIndex=fsmParam.specific.firstIndex;
imageLastIndex=fsmParam.specific.lastIndex;

% Check that the images exist
if exist(imageFileList{1},'file')==0
    warnStr=['I cannot find the image files stored in fsmParam.mat. You will now be asked to manually locate the file ',imageFileList{1},'.'];
    uiwait(msgbox(warnStr,'Warning'));
 
    % The user must select the first image of the stack 
    [tmp0,currentName,currentIndex,currentExt]=getFilenameBody(imageFileList{1});
    openStr=['Please locate ',[currentName,currentIndex,currentExt],'...'];
    [newImageName,newImagePath] = uigetfile(...
        {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
        '*.tif','TIF files (*.tif)'
        '*.tiff','TIFF files (*.tiff)'
        '*.jpg;','JPG files (*.jpg)'
        '*.jpeg;','JPEG files (*.jpeg)'
        '*.*','All Files (*.*)'},...
        openStr);
    if(isa(newImagePath,'char') && isa(newImageName,'char'))
        
        % Check index
        [tmp1,newName,newIndex,newExt]=getFilenameBody([newImagePath,newImageName]);
        
        choice=questdlg(['The selected file does not appear to be the same'...
            ' as pointed to in fsmParam. However, this may happen if you run'...
            ' SpeckTackle and the post-processing in two different operating'...
            ' systems, such as Windows and Linux. Do you want to use the image' ...
            ' you selected?'], 'Warning', 'Yes','No','No');
                   
        switch choice
            case 'Yes', % Continue like this
            case 'No', return
            otherwise, error('Impossible choice in the dialog.');
        end
    
        % Re-create list of files
        imageFileList=getFileStackNames([newImagePath,newImageName]);

        % Cut it
        imageFileList=imageFileList(1:imageLastIndex-imageFirstIndex+1);
    else
        return 
    end
end

% Set status to 1 to indicate success
status=1;