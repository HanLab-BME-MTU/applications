function backupData

%choose top directory to seach for analysis directories in
topDir = uigetdir([],'Choose top directory where data will be copied from');

%choose directory where analysis directories will be copied to
targetDir = uigetdir([],'Choose top directory where data will be saved');

%recursively list subdirectories in this directory
dirList = getSubDirs(topDir);
numDir = length(dirList);

%find the analysis directories
isAnalysisDir = zeros(numDir,1);
for iDir = 1 : numDir
    dirCurr = dirList{iDir};
    fileSepLoc = regexp(dirCurr,'\\');
    tmp1 = regexp(dirCurr(fileSepLoc(end):end),'analysis','once');
    tmp2 = regexp(dirCurr,'oldStuff','once');
    if ~isempty(tmp1) && isempty(tmp2)
        isAnalysisDir(iDir) = 1;
    end
end

%keep only analysis directories
dirList = dirList(logical(isAnalysisDir));
numDir = length(dirList);

%keep only directories one level above analysis directories
keepDir = ones(numDir,1);
for iDir = 1 : numDir - 1
    dirCurr1 = dirList{iDir};
    dirCurr2 = dirList{iDir+1};
    fileSepLoc1 = regexp(dirCurr1,'\\');
    fileSepLoc2 = regexp(dirCurr2,'\\');
    dirCurr1 = dirCurr1(1:fileSepLoc1(end)-1);
    dirCurr2 = dirCurr2(1:fileSepLoc2(end)-1);
    if strcmp(dirCurr1,dirCurr2);
        keepDir(iDir) = 0;
    end
end
dirList = dirList(logical(keepDir));
numDir = length(dirList);
for iDir = 1 : numDir
    dirCurr = dirList{iDir};
    fileSepLoc = regexp(dirCurr,'\\');
    dirCurr = dirCurr(1:fileSepLoc(end)-1);
    dirList{iDir} = dirCurr;
end

%copy the data
for iDir = 1 : numDir
    dirCurr = dirList{iDir};
    disp(dirCurr)
    fileSepLoc = regexp(dirCurr,'\\');
    mkdir([targetDir filesep dirCurr(fileSepLoc(4)+1:end)])
    copyfile([dirCurr filesep 'analysis*'],[targetDir filesep dirCurr(fileSepLoc(4)+1:end)]);
end
