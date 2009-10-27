function fsmUpdateProjDir(rootDirectory)
% FSMUPDATEPROJDIR(rootDirectory) finds every FSMCenter project in a
% rootDirectory and checks that all absolute paths stored in FSMCenter
% parameter files are consistent with the location of those files. If the
% absolute paths are inconsitent, the function replaces those paths with
% the current project location.
%
% For each FSMCenter project, the following files are updated:
%
% lastProjSettings.mat
% edge/parameters.dat (use read_parameters.m and save_parameters.m)
% tack/fsmParam.mat
% tack/fsmImages.mat
%
% This function returns the previous project location stored in the
% parameter files as a string array.
%
% Sylvain Berlemont, 2009

% find every fsmCenter directory in rootDirectory
[status, result] = system(['find ' rootDirectory ' -name lastProjSettings.mat']);

if status
    return;
end

newProjDirList = regexp(result, '\n', 'split')';
if ~numel(newProjDirList{end})
    newProjDirList = newProjDirList(1:end-1);
end
n = length('/lastProjSettings.mat');
newProjDirList = cellfun(@(x) x(1:end-n), newProjDirList,...
    'UniformOutput', false);

for i = 1:numel(newProjDirList)
    newProjDir = newProjDirList{i};

    disp(['Updating ' newProjDir '...']);
    
    % Update lastProjSettings.mat
    
    if newProjDir(end) == filesep && length(newProjDir) ~= 1
        newProjDir = newProjDir(1:end-1);
    end

    filename = [newProjDir filesep 'lastProjSettings.mat'];

    if ~exist(filename, 'file')
        error([filename ' does not exist.']);
    end

    load(filename);
    projDir = projSettings.projDir;
    
    if strcmp(projDir, newProjDir)
        continue;
    end
    
    imgDirList = projSettings.unix_imgDirList{1};

    str1 = projDir;
    str2 = imgDirList;
    n = max(length(str1), length(str2));
    str1(n+1) = ' ';
    str2(n+1) = ' ';
    ind = find(str1 ~= str2, 1);

    if ~isempty(ind)
        projSubDir = projDir(ind:end);
        imgSubDir = imgDirList(ind:end);
    else
        projSubDir = [];
        imgSubDir = [];
    end

    % First case: images are either directly inside projDir or in a sub
    % folder of projDir.
    if isempty(projSubDir)
        projSettings.unix_imgDirList = {[newProjDir imgSubDir]};
    else
        % Second case: fsmCenter project and images are in 2 different
        % directories but have a common parent folder.
        ind = strfind(newProjDir, projSubDir);
        projSettings.unix_imgDirList = {[newProjDir(1:ind-1) imgSubDir]};
    end

    ind = strfind(newProjDir, '/Volumes/');
    if ~isempty(ind) && ind(1) == 1
        projSettings.unixMntRoot = '/Volumes/';
    else
        ind = strfind(newProjDir, '/mnt/');
        if ~isempty(ind) && ind(1) == 1
            projSettings.unixMntRoot = '/mnt/';
        else
            projSettings.unixMntRoot = '/';
        end
    end

    projSettings.unix_imgDrive = {[ '/' strtok(projSettings.unix_imgDirList{1}, ...
        '/')]};

    projSettings.projDir = newProjDir;
    
    save(filename, 'projSettings');
    
    % Update edge/parameters.dat

    edgeSubDir = projSettings.subProjDir{4};
    filename = [newProjDir filesep edgeSubDir filesep 'parameters.dat'];

    if ~exist(filename, 'file')
        disp(['\t' filename ' does not exist (skipping).']);
    else
        parameters = read_parameters(filename, 0);
        parameters.file = [projSettings.unix_imgDirList{1} ...
            projSettings.firstImgList{1}];
        parameters.results = [newProjDir filesep edgeSubDir filesep];
        save_parameters(filename, parameters);
    end
    
    % Update tack/fsmParam.mat
    
    tackSubDir = projSettings.subProjDir{1};
    filename = [newProjDir filesep tackSubDir filesep 'fsmParam.mat'];
    
    if ~exist(filename, 'file')
        disp(['\t' filename 'does not exist (skipping).']);
    else
        load(filename);
        fsmParam.project.path = newProjDir;
        fsmParam.main.path = [newProjDir filesep tackSubDir filesep];
        fsmParam.main.imagePath = [projSettings.unix_imgDirList{1}];
        if fsmParam.track.init
            n = numel(projDir);
            fsmParam.track.initPath = [newProjDir fsmParam.track.initPath(n+1:end)];
        end
        n = numel(imgDirList);
        m = size(fsmParam.specific.fileList, 1);
        fsmParam.specific.fileList = horzcat(repmat(projSettings.unix_imgDirList{1}, ...
            m, 1), fsmParam.specific.fileList(:, n+1:end)); 
        save(filename, 'fsmParam');
    end
    
    % Update tack/fsmImages.mat
    
    filename = [newProjDir filesep tackSubDir filesep 'fsmImages.mat'];
    
    if ~exist(filename, 'file')
        disp(['\t' filename 'does not exist (skipping).']);
    else
        load(filename);
        n = numel(imgDirList);
        fsmImages.firstName = [projSettings.unix_imgDirList{1} ...
            fsmImages.firstName(n+1:end)]; 
        fsmImages.lastName = [projSettings.unix_imgDirList{1} ...
            fsmImages.lastName(n+1:end)];
        save(filename, 'fsmImages');
    end
end