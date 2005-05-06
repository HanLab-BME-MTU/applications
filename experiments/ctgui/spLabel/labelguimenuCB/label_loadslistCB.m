function label_loadslistCB(varargin);
%loads idlist (and slist) into LabelPanel
%changed to load idlist instead of slist - 10/02 JD (in 04/03, slist is loaded again...)



labelguiH = openfig('labelgui','reuse');
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    errordlg('no movie loaded.');
    return;
end;

if nargin==0|isempty(varargin{1}) %no (load-all&exist projectData)
    %check if default biodata-dir exists and cd if exist
    mainDir = cdBiodata(2);
    %change to path from which movie has been loaded
    loadPath = GetUserData(imgFigureH,'moviePath');
    if ~isempty(loadPath)
        cdBiodata(0); %change to biodata-main so that the relative path is valid
        cd(loadPath);
    end
    
    %try to load dataProperties: look for data file with the same name as directoryName
    fSeps = findstr(loadPath,filesep);
    projName = loadPath(fSeps(end-1)+1:fSeps(end)-1);
    dataFile = chooseFile([projName,'-data'],[],'GUI','log');
    if ~isempty(dataFile)
        load(dataFile,'dataProperties');
        SetUserData(imgFigureH,dataProperties,1);
        set(imgFigureH,'Name',dataProperties.name);
    else
        %try to load tmpDataProperties
        altDataFile = chooseFile('tmpDataProperties',[],'new');
        if ~isempty(altDataFile)
            load(altDataFile);
            SetUserData(imgFigureH,dataProperties,1);
            set(imgFigureH,'Name',dataProperties.name);
        else
            error('no dataProperties found!')
        end
    end
    
    %load idlist
    [fname,pathName] = uigetfile('id*','select idlist file');
    if(fname(1)==0)
        image = [];
        fname = [];
        return;
    end;
    cd(pathName);
    filename = [pathName fname];
    idFile = load(fname); %idFile has a field idlist/idlist_L/idlisttrack/idlisttrack_L
    idFileName = char(fieldnames(idFile));
    
    eval(['idlist = idFile.',idFileName,';']);
    idname = idFileName;
    
    
    
    %try to load slist
    slistName = chooseFile('slist',[],'new');
    if ~isempty(slistName)
        load(slistName);
    end
    
else %idlist is being passed down
    idlist = varargin{1};
    idname = varargin{2};
    slist = varargin{3};
end    


%save user data
if exist('idlist','var')|exist('idlisttrack','var')
    
    %check version
    latestVersionDate = GetUserData(labelguiH,'latestChange');
    latestVersionNum = datenum(latestVersionDate);
    if ~isfield(idlist(1).stats,'created')
        error('sorry, you loaded too old an idlist. run the newest version of spotID first')
        return
    elseif latestVersionNum>datenum(idlist(1).stats.created)
        error('sorry, you loaded too old an idlist. run the newest version of spotID first')
        return
    end
    SetUserData(imgFigureH,idlist,1);
    SetUserData(imgFigureH,idlist,1,'idlist_old'); %store 2nd copy of idlist for reverting
    SetUserData(imgFigureH,idname,1);%write name of idlist
    if exist('slist','var')
        SetUserData(imgFigureH,slist,1);
    end
else
    error('incorrect file loaded');
end;

labelgui('refresh');

%clear old 3D data if new idlist is loaded
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    close(view3DH);
end
view3DfigH = GetUserData(imgFigureH,'view3Dfig');
if ishandle(view3DfigH)
    close(view3DfigH);
end

%recall positions
positions = GetUserData(labelguiH,'positions');

%if option is selected, load xzyz-view
labelguiHandles = guidata(labelguiH);
XYXZisChecked = get(labelguiHandles.label_showXZYZ,'Checked');

if strcmp(XYXZisChecked, 'on')
    set(labelguiHandles.label_showXZYZ,'Checked','off'); %else the program thinks the menu's being unchecked
    label_showXZYZ_CB(labelguiHandles.label_showXZYZ,[],labelguiHandles);
    %set position if applicable
    if isfield(positions,'zFigPos')
        zFigH = GetUserData(imgFigureH, 'XZYZFigureH');
        set(zFigH,'Position',positions.zFigPos)
    end
end



%if intFigure is selected: open & update
intFigIsChecked = get(labelguiHandles.label_showIntFig,'Checked');
if strcmp(intFigIsChecked,'on')
    set(labelguiHandles.label_showIntFig,'Checked','off');
    label_showIntFig(labelguiHandles.label_showIntFig,[],labelguiHandles);
    figH = GetUserData(imgFigureH,'intFigH');
    if ishandle(figH)
        figHandles = guidata(figH(end));
        buttonH = findall(0,'Tag','updateIntFig_PB');
        if ishandle(buttonH)
            label_updateIntFigCB(buttonH(end),[],figHandles);
        end
    end
    %set position if applicable
    if isfield(positions,'intFigPos')
        set(figH,'Position',positions.intFigPos)
    end
end

%if disFigure is selected: open & update
disFigIsChecked = get(labelguiHandles.label_showDisFig,'Checked');
if strcmp(disFigIsChecked,'on')
    set(labelguiHandles.label_showDisFig,'Checked','off');
    label_showDisplacementFigure(labelguiHandles.label_showDisFig,[],labelguiHandles);
    %set position if applicable
    if isfield(positions,'disFigPos')
        disFigH = GetUserData(imgFigureH,'disFigH');
        set(disFigH,'Position',positions.disFigPos)
    end
    
end


%if movieDataFigure is selected: open
dataIsChecked = get(labelguiHandles.label_showMovieData,'Checked');
if strcmp(dataIsChecked,'on')
    label_showMovieData(labelguiHandles.label_showMovieData,[],labelguiHandles,1);
    %set position if applicable
    if isfield(positions,'dataFigPos')
        dataFigH = GetUserData(imgFigureH, 'CurrentMovieData');
        set(dataFigH,'Position',positions.dataFigPos)
    end
end