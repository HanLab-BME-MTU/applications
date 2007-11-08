function analysisStruct = makiCollectMovies(jobType)
%MAKICOLLECTMOVIES selects a group of movies that get combined in data analysis
%
%SYNOPSIS analysisStruct = makiCollectMovies(jobType)
%
%INPUT jobType: string which can take the values:
%               'TEST','HERCULES','DANUSER','MERLADI',SWEDLOW' or MCAINSH
%
%OUTPUT analysisStruct: Structure with fields:
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%
%Khuloud Jaqaman, July 2007

%assign default job type if not input
%define top directory to start search
if nargin < 1
    topDir = makiPathDef('$TESTDATA',jobType);
else
    jobType = upper(jobType);
    if strcmp(jobType,'TEST')
        topDir = makiPathDef('$TESTDATA',jobType);
    elseif strcmp(jobType,'HERCULES')        
        topDir = makiPathDef('$HERCULES',jobType);
    elseif strcmp(jobType,'DANUSER')
        topDir = makiPathDef('$DANUSER',jobType);
    elseif strcmp(jobType,'MERALDI')
        topDir = makiPathDef('$MERALDI',jobType);
    elseif strcmp(jobType,'SWEDLOW')
        topDir = makiPathDef('$SWEDLOW',jobType);
    elseif strcmp(jobType,'MCAINSH')
        topDir = makiPathDef('$MCAINSH',jobType);
    end
end

%allow user to choose directory
basePath = uigetdir(topDir,'Please select directory');

%find all makiAnalysis files in chosen directory
fileList1 = searchFiles('makiAnalysis',[],basePath,1);
numList1 = size(fileList1,1);

%find all makiData files in chosen directory
fileList2 = searchFiles('makiData',[],basePath,1);

%put the 2 lists in one
fileList = [fileList1; fileList2];

%allow user to choose files
selectIdx = listSelectGUI(fileList(:,1),[],'move');

%get number of chosen files
numFiles = length(selectIdx);

%keep only chosen files
fileList = fileList(selectIdx,:);

%find which files are makiData and which are makiAnalysis
isDataFile = selectIdx > numList1;

%generate analysisStruct and save it based on number and type of files
if numFiles == 1 %if only 1 file was uploaded

    %determine directory where to save analysisStruct
    dir2SaveRes = fileList{1,2};

    %if selected file is a makiData file ...
    if isDataFile

        %get file name and name base
        fileName = fileList{1,1};
        indx = regexp(fileName,'-makiData');
        fileNameBase = fileName(1:indx-1);

        %give a name to the file where analysisStruct will be saved
        fileName = ['makiAnalysis_' fileNameBase '_1'];

        %generate analysisStruct
        analysisStruct.fileName = fileName;
        analysisStruct.filePath = dir2SaveRes;
        analysisStruct.movies = fileList;
        
        %make analysisStruct platform independent
        analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

        %save analysisStruct
        save(fullfile(dir2SaveRes,fileName),'analysisStruct');
        
    else %if it is not a data file
        
        %get file name from user selection
        fileName = fileList{1,1};
        
        %load analysisStruct
        load(fullfile(dir2SaveRes,fileName));

    end

else %if more than 1 file was uploaded

    %create cell arrays where movie names are stored
    movieListAll1 = [];
    movieListAll2 = [];

    %go over all selected files
    for iFile = 1 : numFiles

        if isDataFile(iFile) %if this file is a data file

            %add movie name to list
            movieListAll1 = [movieListAll1; fileList(iFile,:)];

        else %if this file is an analysis file

            %load analysisStruct
            load(fullfile(fileList{iFile,2},fileList{iFile,1}));
            
            %add movie names to list
            movieListAll2 = [movieListAll2; analysisStruct.movies];

            %delete the uploaded analysisStruct
            clear analysisStruct

        end %(if isDataFile(iFile))

    end %(for iFile = 1 : numFiles)
    
    %change file paths in movieListAll1 into identifiers
    for i = 1 : size(movieListAll1,1)
        movieListAll1{i,2} = makiPathDef(movieListAll1{i,2},jobType);
    end
    
    %put all selected files together
    movieListAll = [movieListAll1; movieListAll2];
    
    %make sure that there are no movie repetitions
    [dummy,indxUnique] = unique(movieListAll(:,1));
    movieListAll = movieListAll(indxUnique,:);
    
    %get directory and file name where to save analysisStruct
    [fileName,dir2SaveRes] = uiputfile('makiAnalysis_*');

    if fileName == 0 %if user hit cancel

        disp('no file saved');

        %generate analysisStruct
        analysisStruct.fileName = fileName;
        analysisStruct.filePath = makiPathDef(dir2SaveRes,jobType);
        analysisStruct.movies = movieListAll;

    else %if user specified a file and directory

        %make sure that file name begins with makiAnalysis
        if length(fileName) < 13 || ~strcmp(fileName(1:13),'makiAnalysis_')
            fileName = ['makiAnalysis_' fileName];
        end

        %check whether the user typed in .mat
        if strcmp(fileName(end-3:end),'.mat')
            fileNameTmp = fileName(1:end-4);
        else
            fileNameTmp = fileName;
            fileName = [fileName '.mat'];
        end
        
        %make sure that file name ends with a version number
        versionNum = makiGetVersion(fileName);
        if isempty(versionNum)
            fileName = [fileNameTmp '_1.mat'];
        end
        
        %generate analysisStruct
        analysisStruct.fileName = fileName;
        analysisStruct.filePath = makiPathDef(dir2SaveRes,jobType);
        analysisStruct.movies = movieListAll;

        %save analysisStruct
        save(fullfile(dir2SaveRes,fileName),'analysisStruct');

    end

end %(if numFiles == 1 ... else ...)

