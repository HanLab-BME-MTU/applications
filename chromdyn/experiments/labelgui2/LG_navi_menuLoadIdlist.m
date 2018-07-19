function success = LG_navi_menuLoadIdlist(hObject,eventData,naviHandles)
%LG_navi_menuLoadIdlist loads idlists into labelgui2

success = 0;

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get movieDir
movieDir = movieWindowHandles.movieDir;

% if this function is called from loadAll, we simply want to read the
% idlists from the data file. If it is called directly from anywhere, we
% want to allow the user to choose an idlist from file. 
if hObject == 0
    % we're coming from loadAll. Find data (or data2) file
    [dataFile,token] = searchFiles('-data(\d?)-','log',movieDir);
    
    % die if no file
    if isempty(dataFile)
        h = errordlg('no dataFile found!','file not found');
        uiwait(h)
        return
    end
    
    % remember dataFile Name, directory
    dataFileName = dataFile{1};
    idlistDir = dataFile{2};
    if isempty(token{1})
        [idlist,dataProperties,projProperties,slist] = ...
            loadProjectData(dataFileName,movieDir,[],1,1);
        if isempty(idlist)
            success = 0;
            return
        end
        % if idlist is -1: offer to recalculate
        if isequal(idlist,-1)
            answer = myQuestdlg('No new idlist found. Do you want to recalculate?',...
                'Warning','Yes','No','Yes');
            if strcmp(answer,'Yes')
                
                %% calculate idlist from slist
                try
                idlist = linker(slist,dataProperties);
                catch
                    err=lasterror;
                    h = errordlg(err.message,'ERROR');
                    uiwait(h);
                    success = 0;
                    return
                end
                idlist(1).stats.idname = 'idlist2';
                
                % save idlist2 in -data-
                projProperties.status = 7; % I hope that's correct!
                save(fullfile(movieDir,dataFileName),'projProperties','-append');
                idlist2=idlist;
                save(fullfile(movieDir,dataFileName),'idlist2','-append');
                lastResult = 'idlist2';
                save(fullfile(movieDir,dataFileName),'lastResult','-append');
                
                % later: write to log
            else
                success = 0;
                return
            end
        end
        
        % remember idname
        idname = idlist(1).stats.idname;
        idlist(1).stats = rmfield(idlist(1).stats,'idname');
        
    else
        % we'll need to update loadData2 to have the choice of idlists
        h = errordlg('need to work on loadData2.m first','unfinished code');
        uiwait(h)
        return
    end
    
else
    % load idlist via selection dialog
    
    % uigetfile doesn't allow passing a directory. Therefore cd
    oldDir = cd(movieDir);
    [fname,idlistDir] = uigetfile('id*','select idlist file');
    cd(oldDir);
    
    idlistDir = idlistDir(1:end-1); % remove filesep
    
    % if no selection, die
    if(fname(1)==0)
        h = errordlg('no idlist loaded!');
        uiwait(h)
        return
    end
    %idFile has a field idlist/idlist_L/idlisttrack/idlisttrack_L. Read it
    %and assign idlist
    if ~strcmp(idlistDir(end),filesep)
        idlistDir = [idlistDir, filesep];
    end
    idFile = load([idlistDir,fname]); 
    idFileName = char(fieldnames(idFile));
    idlist = idFile.(idFileName);
    idname = idFileName;
    
    % assign empty dataFileName
    dataFileName = [];
end

% load idlist and set all the rest. Replace whatever idlist was there
success = LG_loadIdlist(idlist, 1, idname, idlistDir, dataFileName);


