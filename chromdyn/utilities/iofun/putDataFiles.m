function putDataFiles
%copies movies and data from one biodata-directory to another
%
%SYNOPSIS  putDataFiles
%
%INPUT     the program asks for a source and a destination directory. The
%          source directory should contain all the project directories 
%          from which you want to copy files, the destination directory
%          should be the biodata top directory.
%
%OUTPUT    the progam will copy the data file and the corresponding
%          logfile. If necessary, new subdirectories are added and movies
%          and filtered movies are copied
%
%03/03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remember old path
oldPath=pwd;

%change to biodata
mainDir=cdBiodata(0);

mainDir=uigetdir(mainDir,'Select biodata directory from which to copy data');
mainDirLength=length(mainDir);

pathNameGet=uigetdir(mainDir,'Select directory containing project directories from which to copy data');

pathNamePut=uigetdir('','Select biodata directory to which to copy data');

%build dirList
dirList=dir(pathNameGet);

nDir=size(dirList,1);
copyList=cell(1,3);

%loop through all project directories and build list of datafiles (and
%filtered movies) to be copied
waitbarHandle=mywaitbar(0,[],nDir,'reading file structure');
for i=3:nDir
    cd(pathNameGet);
    if dirList(i).isdir
        cd(dirList(i).name);
        %get relative path
        absPath=pwd;
        relPath=absPath(mainDirLength+2:end);
        %choose dataFile
        datafileName=chooseFile([dirList(i).name,'-data'],[],'new','log');
        %find filtered movie name
        movieName=chooseFile('filtered_movie',[],'new');
        if isempty(movieName)
            movieName=chooseFile('moviedat',[],'new');
        end
        if ~isempty(movieName)
            copyList{i-2,1}=relPath;
            copyList{i-2,2}=datafileName;
            copyList{i-2,3}=movieName;
            copyList{i-2,4}=dirList(i).name;
        end
    end
    mywaitbar(i/(nDir),waitbarHandle,nDir);
end
close(waitbarHandle);


%loop through the copyList and copy all files (and if necessary, mkdir)
cd(mainDir);
for i=1:size(copyList,1)
    if ~isempty(copyList(i,1))
        pNPi=[pathNamePut,filesep,copyList{i,1}];
        %create dir if necessary
        if ~isdir(pNPi)
            cd(pathNamePut);
            mkdir(copyList{i,1});
            cd(mainDir);
        end
        %copy data file if there is any
        if ~isempty(copyList{i,2})
            disp([nowString,' copying ',copyList{i,1},filesep,copyList{i,2}]);
            copyfile([copyList{i,1},filesep,copyList{i,2}],[pNPi,filesep,copyList{i,2}]);
            %copy log file if it exists and if the data file exists
            if ~isempty(dir([copyList{i,1},filesep,copyList{i,2},'.log']))
                disp([nowString,' copying ',copyList{i,1},filesep,copyList{i,2},'.log']);
                copyfile([copyList{i,1},filesep,copyList{i,2},'.log'],[pNPi,filesep,copyList{i,2},'.log']);
            end
        end
        %copy filtered movie if necessary
        if isempty(dir([pNPi,filesep,copyList{i,3}]))
            disp([nowString,' copying ',copyList{i,1},filesep,copyList{i,3}]);
            copyfile([copyList{i,1},filesep,copyList{i,3}],[pNPi,filesep,copyList{i,3}]);
        end
        %copy movie if necessary
        if isempty(dir([pNPi,filesep,copyList{i,4},'.r3*']))
            disp([nowString,' copying ',copyList{i,1},filesep,copyList{i,4},'.r3x']);
            copyfile([copyList{i,1},filesep,copyList{i,4},'.r3*'],[pNPi]);
            if ~isempty(dir([copyList{i,1},filesep,copyList{i,4},'.r3d.log']))
                copyfile([copyList{i,1},filesep,copyList{i,4},'.r3d.log'],[pNPi,filesep,copyList{i,4},'.r3d.log']);
            end    
        end
    end
end

disp([nowString,' done']);
cd(oldPath);