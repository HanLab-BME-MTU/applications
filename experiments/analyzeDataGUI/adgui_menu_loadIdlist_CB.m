function adgui_menu_loadIdlist_CB(hObject, eventdata, handles)
% loads different idlist for the same project
% hObject    handle to adgui_menu_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check that a valid movieFile has been selected
currentSelection = get(handles.adgui_filename_PD,'Value');

%the PD can not remain on 'add Data'
if currentSelection==1
    h = errordlg('You first have to load a regular data file before loading alternate idlists!')
    uiwait(h)
    return
end

%load data for current movie and switch to directory
currentData = handles.data(currentSelection-1);
cd(currentData.dataPath);

%ask for input
fileName=chooseFile('idlist',[],'GUI');

%load file (idlist or idlisttrack)
if ~isempty(fileName)
load(fileName);
else 
    return
end

%turn loaded data into idlist
if strmatch('idlisttrack',fileName);
    if exist('idlisttrack','var')
        idlist = idlisttrack;
    elseif exist('idlisttrack_L','var')
        idlist = idlisttrack_L;
    else
        h = errordlg('bad idlisttrack');
        uiwait(h)
        return
    end
end

%read dataProperties
dataProperties = currentData.dataProperties;
%label idlist
lastResult = fileName;

%-------------calculate anaDat
anaDat = adgui_calc_anaDat(idlist,dataProperties,lastResult);
%-----------------------------

%store anaDat into adguiHandles
numData = size(handles.data,2);
handles.data(numData+1).anaDat = anaDat;

%store dataProperties
handles.data(numData+1).dataProperties = dataProperties;

%store path
handles.data(numData+1).dataPath = currentData.dataPath;

%write moviename into "current file text", add 'alt'
oldText = get(handles.adgui_filename_PD,'String');
set(handles.adgui_filename_PD,'String',[oldText;{['alt ',anaDat(1).info.name]}],...
    'Value',size(oldText,1)+1,'ForegroundColor',extendedColors(handles.colorList(mod(numData+1-1,size(handles.colorList,1))+1,:)),'BackgroundColor','w');

%set tag text
set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles;handles.cenHandles(2);handles.alignHandles(2:3)],...
    'Value',1,'String',[{'set tag'};anaDat(1).info.labelColor]);

%set tag text
set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles;handles.cenHandles(2);handles.alignHandles(2:3)],...
    'Value',1,'String',[{'set tag'};anaDat(1).info.labelColor]);

guidata(hObject,handles);