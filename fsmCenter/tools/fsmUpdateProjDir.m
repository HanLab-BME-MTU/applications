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
% This function works either from or to windows and linux.
%
% Sylvain Berlemont, 2009

% find every fsmCenter directory in rootDirectory
if rootDirectory(end) == filesep
    rootDirectory = rootDirectory(1:end-1);
end
[status, result] = system(['find ' rootDirectory ' -name lastProjSettings.mat']);

if status
    return;
end

newProjDirList = regexp(result, '\n', 'split')';
if ~numel(newProjDirList{end})
    newProjDirList = newProjDirList(1:end-1);
end
n = length([filesep 'lastProjSettings.mat']);
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
        disp('   lastProjSettings.mat does not exist.');
        continue;
    end

    load(filename);
    projDir = projSettings.projDir;
    
    if strcmp(projDir, newProjDir)
        continue;
    end
    
    createdOnLinux = ~isempty(projSettings.unix_imgDirList);
    createdOnWindows = ~isempty(projSettings.win_imgDirList);
    
    if ~xor(createdOnLinux, createdOnWindows)
        disp('  Invalid projSettings.*_imgDirList (skipping).');
        continue;
    end
       
    if createdOnLinux
        imgDirList = projSettings.unix_imgDirList{1};
    else
        imgDirList = projSettings.win_imgDirList{1};
    end

    str1 = projDir;
    str2 = imgDirList;
    n = max(length(str1), length(str2));
    str1(n+1) = ' ';
    str2(n+1) = ' ';
    ind = find(str1 ~= str2, 1);

    if ~isempty(ind)
        projSubDir = projDir(ind:end);
        imgSubDir = imgDirList(ind:end);
        c = imgSubDir(end);
        if c == '/' || c == '\'
            imgSubDir(end) = filesep;
        end
    else
        projSubDir = [];
        imgSubDir = [];
    end

    % First case: images are either directly inside projDir or in a sub
    % folder of projDir.
    if isempty(projSubDir)
        newImgDirList = [newProjDir imgSubDir];
    else
        % Second case: fsmCenter project and images are in 2 different
        % directories but have a common parent folder.
        ind = strfind(newProjDir, projSubDir);
        newImgDirList = [newProjDir(1:ind-1) imgSubDir];
    end
    
    % update projSettings.*_imgDirList
    if ispc
        projSettings.win_imgDirList = {newImgDirList};
        projSettings.unix_imgDirList = {};
    else
        projSettings.win_imgDirList = {};
        projSettings.unix_imgDirList = {newImgDirList};
    end
    
    % update projSettings.unixMntRoot
    if ~ispc
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
    end

    % update projSettings.*_imgDrive
    if ispc
        projSettings.win_imgDrive = {strtok(newImgDirList, filesep)};
        projSettings.unix_imgDrive = {};
    else
        projSettings.win_imgDrive = {};
        projSettings.unix_imgDrive = {[ '/' strtok(newImgDirList, ...
            filesep)]};
    end

    % update projSettings.projDir
    projSettings.projDir = newProjDir;
    
    % update projSettings.firstImgList
    if numel(projSettings.firstImgList) > 1
        projSettings.firstImgList = {projSettings.firstImgList{1}};
    end
    
    save(filename, 'projSettings');
    
    % Update edge/parameters.dat

    edgeSubDir = projSettings.subProjDir{4};
    filename = [newProjDir filesep edgeSubDir filesep 'parameters.dat'];

    if ~exist(filename, 'file')
        disp('   parameters.dat does not exist (skipping).');
    else
        parameters = read_parameters(filename, 0);
        parameters.file = [newImgDirList projSettings.firstImgList{1}];
        parameters.results = [newProjDir filesep edgeSubDir filesep];
        save_parameters(filename, parameters);
    end
    
    % Update tack/fsmParam.mat
    
    tackSubDir = projSettings.subProjDir{1};
    filename = [newProjDir filesep tackSubDir filesep 'fsmParam.mat'];
    
    if ~exist(filename, 'file')
        disp('   fsmParam.mat does not exist (skipping).');
    else
        load(filename);
        fsmParam.project.path = newProjDir;
        fsmParam.main.path = [newProjDir filesep tackSubDir filesep];
        fsmParam.main.imagePath = newImgDirList;
        if fsmParam.track.init
            n = numel(projDir);
            fsmParam.track.initPath = [newProjDir fsmParam.track.initPath(n+1:end)];
        end
        n = numel(imgDirList);
        
        if iscell(fsmParam.specific.fileList)
            fsmParam.specific.fileList = cellfun(@(c) [newImgDirList ...
                c(n+1:end)], fsmParam.specific.fileList);
        elseif ischar(fsmParam.specific.fileList)
            m = size(fsmParam.specific.fileList, 1);
            
            fsmParam.specific.fileList = arrayfun(@(i) [newImgDirList ...
                strtrim(fsmParam.specific.fileList(i, n+1:end))], 1:m, ...
                'UniformOutput',false);
        else
            disp('   Invalid fsmParam.specific.fileList field (skipping).');
        end
        save(filename, 'fsmParam');
    end
    
    % Update tack/fsmImages.mat
    
    filename = [newProjDir filesep tackSubDir filesep 'fsmImages.mat'];
    
    if ~exist(filename, 'file')
        disp('   fsmImages.mat does not exist (skipping).');
    else
        load(filename);
        n = numel(imgDirList);
        fsmImages.firstName = [newImgDirList fsmImages.firstName(n+1:end)]; 
        fsmImages.lastName = [newImgDirList fsmImages.lastName(n+1:end)];
        save(filename, 'fsmImages');
    end
end