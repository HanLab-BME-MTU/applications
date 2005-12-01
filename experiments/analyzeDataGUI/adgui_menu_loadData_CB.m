function adgui_menu_loadData_CB(hObject, eventdata, handles)
% loads idlist and dataProperties into analyzeDataGUI
% hObject    handle to adgui_menu_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%check if default biodata-dir exists and cd if exist
mainDir = cdBiodata(2);

%get project data file
[fileName,pathName] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files'},'select project data file');

if fileName==0;
    %user has to load data manually
    
    %not implemented yet
    return
else
    cd(pathName);
    %load dataProperties
    load(fileName,'dataProperties');
    if ~exist('dataProperties','var')
        h = errordlg('No dataProperties in project data: corrupt data file','Warning!');
        return %end evaluation here
    end
    
    %get lastResult, load corresponding idlist (if possible)
    load(fileName,'lastResult');
    if isempty(findstr('idlist',lastResult))
        errordlg('No idlist of any kind in project data: run spotID first','Warning!');
    end
    
    %load either idlist, idlist_L, idlisttrack, or idlisttrack_L
    load(fileName,lastResult);
    
    %name loaded data idlist
    eval(['idlist = ',lastResult,';']);
    
end

%-------------calculate anaDat
anaDat = adgui_calc_anaDat(idlist,dataProperties,lastResult);
%-----------------------------

%store anaDat into adguiHandles
numData = size(handles.data,2);
handles.data(numData+1).anaDat = anaDat;

%store dataProperties
handles.data(numData+1).dataProperties = dataProperties;


%write moviename into "current file text"; make pulldown-menu if necessary
switch numData 
case 0
	set(handles.adgui_filename_txt,'String',anaDat(1).info.name);
case 1 %create PD
	oldText = get(handles.adgui_filename_txt,'String');
	set(handles.adgui_filename_txt,'Style','popupmenu','String',{oldText;anaDat(1).info.name},...
	'Value',2,'ForegroundColor',extendedColors(handles.colorList(2,:)),'BackgroundColor','w');
otherwise %it's already a PD
	oldText = get(handles.adgui_filename_txt,'String');
	set(handles.adgui_filename_txt,'String',[oldText;{anaDat(1).info.name}],...
	'Value',size(oldText,1)+1,'ForegroundColor',extendedColors(handles.colorList(numData+1,:)),'BackgroundColor','w');
end

%set tag text
set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles;handles.cenHandles(2);handles.alignHandles(2:3)],...
    'Value',1,'String',[{'set tag'};anaDat(1).info.labelColor]);

guidata(hObject,handles);