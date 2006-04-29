function LG_navi_menuSaveIdlist(hObject,eventdata,handles)
% LG_navi_menuSaveIdlist saves idlists and dataProperties from labelgui2

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% if loadedFromOutside, there is no dataFileName - the outside function
% has to figure out what and where to save
% check if runCtBatch is active
if naviHandles.launchedFromOutside
    % set dataHasChanged to saved
    % dataHasChanged:
    % 0 : same as original
    % 1 : different from original, different from last saved
    % 2 : same as original, different from saved
    % 3 : different from original, saved
    if movieWindowHandles.dataHasChanged
        movieWindowHandles.dataHasChanged = 3;
    end
    % let runCtBatch resume: it will read idlist and then close the
    % movieWindow. Don't forget to store dataHasChanged
    guidata(movieWindowHandles.LG_movieWindow, movieWindowHandles);
    uiresume(naviHandles.LG_navigator)
    return
end

% to save, we need to do the following:
% 1) get idname
% 2) find dataFile (figure out whether it's data or data2)
% 3) save into dataFile
% 4) save into movieDir
% 5) adjust dataHasChanged

% get new idname. IdSaveName is the name of the idlist on disk
idlist = movieWindowHandles.idlist;
idname = movieWindowHandles.idname;
if strfind(idname,'track')
    idname = 'idlisttrack_L2';
    idSaveName = ['idlisttrack-L2-',nowString];
else
    idname = 'idlist_L2';
    idSaveName = ['idlist-L2-',nowString];
end
% name the variable idlist correctly
eval(sprintf('%s = idlist;',idname));

% find dataFile
dataFileName = movieWindowHandles.dataFileName;
if ~isempty(dataFileName)
    % add path to dataFileName. Check for correct idlistDir b/c of bug
    % (04/06)
    if isdir(movieWindowHandles.idlistDir)
        dataFileName = fullfile(movieWindowHandles.idlistDir,dataFileName);
    else
        dataFileName = fullfile(movieWindowHandles.movieDir,dataFileName);
    end
    save(dataFileName,idname,'-append');

    %update project properties
    load(dataFileName,'projProperties');
    prevStatus = bsum2bvec(projProperties.status);
    if strcmp(idname,'idlist_L2')
        projProperties.status = sum(prevStatus(find(prevStatus<8)))+8;
    else %it's an idlisttrack
        projProperties.status = sum(prevStatus(find(prevStatus<32)))+32;
    end
    %save project properties
    save(dataFileName,'projProperties','-append');

    %update lastResult
    load(dataFileName,'lastResult');
    lastResult = idname;
    save(dataFileName,'lastResult','-append');

    %save outside dataFile
    save([movieWindowHandles.idlistDir,idSaveName],idname);
    
    % save dataProperties
    if movieWindowHandles.dataPropertiesHasChanged
        dataProperties = movieWindowHandles.dataProperties;
        save(dataFileName,'dataProperties','-append');
    end

else
    % if user chose idlist, let user choose saveName/saveDir. Goto idlist
    % directory first
    oldDir = cd(movieWindowHandles.idlistDir);
    [fname,pathName] = uiputfile(idSaveName,'save idlist');
    if fname == 0
        % user cancelled
        return
    end
    save([pathName,fname],idname);
    
    % save dataProperties
    if movieWindowHandles.dataPropertiesHasChanged
        dataProperties = movieWindowHandles.dataProperties;
        save([pathName,'tmpDataProperties'],'dataProperties');
    end
    
    cd(oldDir);
end

% set dataHasChanged to saved
% dataHasChanged:
% 0 : same as original
% 1 : different from original, different from last saved
% 2 : same as original, different from saved
% 3 : different from original, saved
if movieWindowHandles.dataHasChanged
    movieWindowHandles.dataHasChanged = 3;
end
% set movieWindowHandles.dataPropertiesHasChanged to 0 in any case
movieWindowHandles.dataPropertiesHasChanged = 0;

guidata(movieWindowHandles.LG_movieWindow, movieWindowHandles);