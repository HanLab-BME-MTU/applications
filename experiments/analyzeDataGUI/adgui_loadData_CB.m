function adgui_loadData_CB(hObject, eventdata, handles)
% loads idlist and dataProperties into analyzeDataGUI

%hObject is the handle to the pulldown!

%get anaDat
try
    [anaDat,pathName,dataProperties] = calcAnaDatFromProjectData;
catch
    errormsg = lasterr;
    if strcmp(lasterr,'No file loaded')
        %reset only
        set(hObject,'Value',handles.lastFile);
        return
    else
        %return error
        h = errordlg(errormsg,'Warning'),
        uiwait(h)
        return
    end
end

%store anaDat into adguiHandles
numData = size(handles.data,2);
handles.data(numData+1).anaDat = anaDat;

%store dataProperties
handles.data(numData+1).dataProperties = dataProperties;

%store pathName
handles.data(numData+1).dataPath = pathName;

%write moviename into "current file text"; delete first entry if necessary
oldText = get(hObject,'String');
if strcmp(oldText{1},'no file loaded')
    oldText(1)=[];
end
set(hObject,'String',[oldText;{anaDat(1).info.name}],...
    'Value',size(oldText,1)+1,'ForegroundColor',extendedColors(handles.colorList(mod(numData+1-1,size(handles.colorList,1))+1,:)),'BackgroundColor','w');
handles.lastFile = size(oldText,1)+1;
%set tag text
set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles;handles.cenHandles(2);handles.alignHandles(2:3)],...
    'Value',1,'String',[{'set tag'};anaDat(1).info.labelColor]);

guidata(hObject,handles);