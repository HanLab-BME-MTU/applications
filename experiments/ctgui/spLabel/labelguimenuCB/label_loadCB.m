function label_loadCB
%loads both moviefile and idlist; tries to load from data file

%check if default biodata-dir exists and cd if exist
mainDir = cdBiodata(2);

%get project data file
[fileName,pathName,loadIdx] = uigetfile(...
    {'*-data-??-???-????-??-??-??.mat','project data files';...
    '*-data2-??-???-????-??-??-??.mat','project data files'},...
    'select project data file');

if fileName==0;
    %user has to load data manually
    label_loadmovieCB;
    label_loadslistCB;
    return
else
    cd(pathName);
    switch loadIdx
        case 1
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

            %find which idlists there are
            dataFieldNames = fieldnames(data);
            idnameListIdx = strmatch('idlist',dataFieldNames);
            idnameList = dataFieldNames(idnameListIdx);

            %have the user choose, if there is more than one entry left
            switch length(idnameList)
                case 0 %no idlist loaded. continue w/o loading

                    idname = '';
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


        case 2
            load(fileName)

            load(data2File.dataProperties)
            if ~isempty(data2File.slist)
                load(data2File.slist)
            else
                load(data2File.synthSlist)
            end

            % have user choose idlists, if there are any
            dataFieldNames = fieldnames(data2File);
            idnameListIdx = [strmatch('idlist',dataFieldNames);...
                strmatch('synthIdlist',dataFieldNames)];
            idnameList = dataFieldNames(idnameListIdx);
            % make sure the fields aren't empty
            for i = size(idnameList,1):-1:1
                if isempty(data2File.(idnameList{i}))
                    idnameList(i) = [];
                    idnameListIdx(i) = [];
                end
            end

 %have the user choose, if there is more than one entry left
            switch length(idnameList)
                case 0 %no idlist loaded. continue w/o loading

                    idname = '';
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
            
            % load idlist
            data = load(data2File.(idname));
            idname = char(fieldnames(data));

        otherwise
            h=warndlg('sorry, can''t just load any file')
            uiwait(h)
            return % end evaluation here
    end


    %--------------try to load filtered movie
    % use cdLoadMovie
    [filteredMovie,dummy,loadStruct] = cdLoadMovie('latest');
    if ~strcmp(loadStruct.movieType,'filtered')
        warning(sprintf('filtered movie not found. Loading %s movie instead',loadStruct.movieType));
    end
    movieName = loadStruct.movieName;



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

label_loadmovieCB(filteredMovie,movieName,pathName);

%set dataProperties (do before idlist is loaded to prevent using refresh/plotlabels with wrong pixelsize)
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
SetUserData(imgFigureH,dataProperties,1);

%set save file name
dataFile.name = fileName;
dataFile.path = pathName;
SetUserData(imgFigureH,dataFile,1);

if ~isempty(idname) % do not load idlist if there isn't any!
    %load idlist
    label_loadslistCB(data.(idname),idname,slist);
end