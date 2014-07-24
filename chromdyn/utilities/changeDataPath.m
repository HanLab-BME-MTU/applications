function changeDataPath
%utility to change the projectProperties.dataPath information

%get mainDir
mainDir = cdBiodata(0);
mainLength = length(mainDir);
%get list of all dataFiles that should be changed
listOfDataFiles = searchFiles('-data-','log','ask');

%make sure user selected something
if isempty(listOfDataFiles)
    error('there are no dataFiles in the selected directory (incl. subdirectories)')
end

%ask user what to do (update path or change fileseps)
whichTask = questdlg('What do you want to change?','Select Task','updatePath','windows2Linux','linux2Windows','updatePath');

switch whichTask
    case 'updatePath'
    update = 1;
case 'windows2Linux'
    oldFileSep = '\';
    newFileSep = '/';
    update = 0;
case 'linux2Windows'
    oldFileSep = '/';
    newFileSep = '\';
    update = 0;
otherwise %user cancelled
    return
end

%loop through the dataFiles and perform task
for i = 1:size(listOfDataFiles,1)
    %go to the right directory
    cd(listOfDataFiles{i,2});
    %load projProperties
    load(listOfDataFiles{i,1},'projProperties');
    
    switch update
        case 1
            %update file location - check whether to use relative or absolute path
            newDir = listOfDataFiles{i,2};
            if strmatch(lower(mainDir),lower(newDir))
                projProperties.dataPath = newDir(mainLength+2:end);
            else
                projProperties.dataPath = newDir;
            end
        case 0
            %replace all oldFileSeps with newFileSeps
            fileSepIdx = strfind(projProperties.dataPath,oldFileSep);
            projProperties.dataPath(fileSepIdx)=newFileSep;
    end
    %save data
    save(listOfDataFiles{i,1},'projProperties','-append');
end %for i = 1:size(listOfDataFiles,1)
