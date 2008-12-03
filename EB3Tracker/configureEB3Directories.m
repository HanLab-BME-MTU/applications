function [newFolderStatus]=configureEB3Directories
% moves Claudio's images corresponding to one movie into its own directory
% assume filenames have following format:
% targetName_oligoNumber_date_movieNumber_miscField_miscField_frame.tif


topDir=uigetdir(pwd,'Please select top-level directory containing targets');
[listOfFiles] = searchFiles('.tif',[],topDir,1); % time consuming...

counter=0;  % for number of files dealt with
counter2=1; % for number of new folders
newFolderStatus=[];
while counter<size(listOfFiles,1) 

    fileName=listOfFiles{counter+1,1};
    x=strfind(fileName,'_');

    % assume filenames have following format:
    % targetName_oligoNumber_date_movieNumber_miscField_miscField_frame.tif
    if length(x)~=6
        disp(['problem with ' listOfFiles(1,2) filesep listOfFiles(1,1)])
        return
    end

    % 1 for all files corresponding to the same movie
    temp=cellfun(@(y) ~isempty(strfind(y,fileName(1:x(4)-1))),listOfFiles(:,1));
    
    newDirName=[topDir filesep 'newDir' filesep ... % topLevel/newDir
        fileName(x(2)+1:x(3)-1) filesep ...         % /date
        fileName(1:x(1)-1) filesep ...              % /target name
        fileName(x(1)+1:x(2)-1) filesep ...         % /oligo reference number
        fileName(x(3)+1:x(4)-1) filesep ...         % /movie number
        'images'];                                  % /images
    
    % if directory doesn't exist yet, make it
    if exist(newDirName)~=7 
        mkdir(newDirName)
    end

    % check to make sure destination directory is different than source and
    % that there actually are some files to move
    [temp2] = searchFiles(fileName(1:x(4)-1),[],listOfFiles{counter+1,2},0);
    if ~isequal(listOfFiles{counter+1,2},newDirName) && ~isempty(temp2)
        movefile([listOfFiles{counter+1,2} filesep fileName(1:x(4)-1) '*'],newDirName)
        newFolderStatus=[newFolderStatus; counter2 sum(temp)];
    end
    
    counter=counter+sum(temp);
    counter2=counter2+1;
    

end