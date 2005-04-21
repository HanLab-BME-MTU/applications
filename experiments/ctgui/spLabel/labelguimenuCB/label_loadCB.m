function label_loadCB
%loads both moviefile and idlist; tries to load from data file

%check if default biodata-dir exists and cd if exist
mainDir = cdBiodata(2);

%get project data file
[fileName,pathName] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files'},'select project data file');

if fileName==0;
    %user has to load data manually
    label_loadmovieCB;
    label_loadslistCB;
    return
else
    cd(pathName);
    data = load(fileName); %loads everything into the structure data
    if ~isfield(data,'dataProperties')
        h = warndlg('No dataProperties in project data: corrupt data file','Warning!');
        uiwait(h);
        return %end evaluation here
    else
        dataProperties = data.dataProperties;
    end
    
    %load slist
    if ~isfield(data,'slist')
        h = warndlg('No slist in project data!','Warning!');
        uiwait(h);
        return %end evaluation here
    else
        slist = data.slist;
    end
    
    %load projectProperties
    if ~isfield(data,'projProperties')
        h = warndlg('No projProperties in project data!','Warning!');
        uiwait(h);
        return %end evaluation here
    else
        projProperties = data.projProperties;
    end
    
    
    
    %--------------try to load filtered movie 
    % use cdLoadMovie
        [filteredMovie,dummy,loadStruct] = cdLoadMovie('latest');
        if ~strcmp(loadStruct.movieType,'filtered')
            warning(sprintf('filtered movie not found. Loading %s movie instead',loadStruct.movieType));
        end
        movieName = loadStruct.movieName;
        
%     %try to find filenames in the path from which projectData has been loaded
%     filteredMovieName = chooseFile('filtered_movie',[],'new');
%     altFilteredMovieName = chooseFile('moviedat',[],'new');
%     if isempty(filteredMovieName)
%         if isempty(altFilteredMovieName) %to ensure compatibility with earlier versions
%             disp('no filtered movie found. load unfiltered movie instead')
%             if findstr(projProperties.dataPath(end-10:end),'crop')|findstr(projProperties.dataPath(end-10:end),'corr')
%                 %cropped movie
%                 moviename = chooseFile('.r3c');
%                 filteredMovie  =  readmat(moviename);
%             else
%                 %normal movie
%                 moviename = chooseFile('.r3d');
%                 filteredMovie  =  r3dread(moviename);
%             end
%         else
%             filteredMovie = readmat(altFilteredMovieName);
%         end
%     else 
%         filteredMovie = readmat(filteredMovieName);
%     end;
%     
%     %test if everything correctly loaded
%     if ~exist('filteredMovie','var')
%         error('no movie found')
%         return
%     end
    %---let the user choose which idlist to load
    
    %find which idlists there are
    dataFieldNames = fieldnames(data);
    idnameListIdx = strmatch('idlist',dataFieldNames);
    idnameList = dataFieldNames(idnameListIdx);
    
    %have the user choose, if there is more than one entry left
    switch length(idnameList)
        case 0 %no idlist loaded. continue w/o loading
            
             idname = '[]';
            h = warndlg('No idlist found in project data','Warning!');
            uiwait(h);
            
        case 1 %only one idlist loaded. Continue
            
            idname = char(idnameList);
            
        otherwise %let the user choose
            idSelect = chooseFileGUI(idnameList);
            if isempty(idSelect)
                idname = '';
            else
                idname = idnameList{idSelect};
            end
    end
end

if isempty(idname)
    ans = myQuestdlg('Do you want to load an idlist?','No idlist loaded from project file','Yes','No','Don''t load anything','Yes');
    switch ans
        case 'Yes'
            label_loadslistCB;
        case 'No'
            % just no idlist, so continue
        otherwise
            % quit here, no load
            return
    end
end

label_loadmovieCB(filteredMovie,projProperties.projName,pathName);

%set dataProperties (do before idlist is loaded to prevent using refresh/plotlabels with wrong pixelsize)
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
SetUserData(imgFigureH,dataProperties,1);

%set save file name
dataFile.name = fileName;
dataFile.path = pathName;
SetUserData(imgFigureH,dataFile,1);

if ~isempty(idname) % do not load idlist if there isn't any!
    %load idlist
    eval(['label_loadslistCB(data.',idname,',idname,slist);']);
end