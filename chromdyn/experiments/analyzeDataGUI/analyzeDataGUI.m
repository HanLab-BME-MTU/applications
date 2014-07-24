function varargout = analyzeDataGUI(varargin)
% ANALYZEDATAGUI M-file for analyzeDataGUI.fig
%      ANALYZEDATAGUI, by itself, creates a new ANALYZEDATAGUI or raises the existing
%      singleton*.
%
%      H = ANALYZEDATAGUI returns the handle to a new ANALYZEDATAGUI or the handle to
%      the existing singleton*.
%
%      ANALYZEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZEDATAGUI.M with the given input arguments.
%
%      ANALYZEDATAGUI('Property','Value',...) creates a new ANALYZEDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyzeDataGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyzeDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyzeDataGUI

% Last Modified by GUIDE v2.5 16-Jul-2003 16:44:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzeDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzeDataGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before analyzeDataGUI is made visible.
function analyzeDataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyzeDataGUI (see VARARGIN)

% Choose default command line output for analyzeDataGUI
handles.output = hObject;

%-----------------------initialize GUI

%collect handles
handles.labelHandles = [handles.adgui_xLabel_PD;...
        handles.adgui_yLabel_PD;...
        handles.adgui_zLabel_PD];

handles.xTagHandles = [handles.adgui_xTag1_PD;...
        handles.adgui_xTag2_PD;...
        handles.adgui_xTag3_PD;...
        handles.adgui_xTag4_PD];

handles.yTagHandles = [handles.adgui_yTag1_PD;...
        handles.adgui_yTag2_PD;...
        handles.adgui_yTag3_PD;...
        handles.adgui_yTag4_PD];

handles.zTagHandles = [handles.adgui_zTag1_PD;...
        handles.adgui_zTag2_PD;...
        handles.adgui_zTag3_PD;...
        handles.adgui_zTag4_PD];

handles.cenHandles = [handles.adgui_cen_txt;...
        handles.adgui_cenTag1_PD];

handles.alignHandles = [handles.adgui_align_txt;...
        handles.adgui_alignTag1_PD;...
        handles.adgui_alignTag2_PD];


%get data for pulldown menus
handles.PD_data = adgui_PD_data;

%set labelMenus and hide tag-PD (none -> no )
set(handles.labelHandles,'String',handles.PD_data(:,1));
set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles],'Visible','off');

%hide cen/align options
set([handles.cenHandles; handles.alignHandles],'Visible','off');

%init fields
handles.plotFigureH = [];
handles.data=[];

%remember last selection of data PD
handles.lastSelected = ones(3,1);

%remember last selection of file
handles.lastFile = 1;

%set possible linestyles
handles.lineStyles = {'-.','x';'-','d';':','+';'--','*'};

%set the color order (list for extendedColors)
handles.colorList = ...
    [1,0;...
        2,0;...
        3,0;...
        11,-0.5;...
        5,-0.5;...
        13,-0.5;...
        6,0;...
        17,0;...
        18,-0.5;...
        23,0;...
        13,0;...
        14,0;...
        15,0;...
        12,0;...
        10,-0.5;...
        22,0];

%for playing games
handles.isFirst = 1;
handles.clicks = 0;

% Update handles structure
guidata(hObject, handles);

%---------------------------------------------------------------------
%------------------------xyz pulldown callbacks------------------------

% --- Executes on selection change in adgui_xLabel_PD.
function adgui_xLabel_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_xLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%activate necessary tag-PD's
selectedItem = get(hObject,'Value');

necessaryTags=handles.PD_data{selectedItem, 2};
cenOrAlign = handles.PD_data{selectedItem, 5};

if necessaryTags == -1
    %user selected a separator
    set(hObject,'Value',handles.lastSelected(1));
else
    set(handles.xTagHandles,'Visible','off');
    set(handles.xTagHandles(1:necessaryTags),'Visible','on');
    handles.lastSelected(1) = selectedItem;
    
    %turn on cen/align
    switch cenOrAlign
        case 0 
            set([handles.cenHandles;handles.alignHandles],'Visible','off');
        case 1
            set(handles.cenHandles,'Visible','on');
            set(handles.alignHandles,'Visible','off');
        case 2
            set([handles.cenHandles;handles.alignHandles],'Visible','on');
    end
    
    guidata(hObject,handles);
end


% --- Executes on selection change in adgui_yLabel_PD.
function adgui_yLabel_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_yLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%activate necessary tag-PD's
selectedItem = get(hObject,'Value');

necessaryTags=handles.PD_data{selectedItem,2};
cenOrAlign = handles.PD_data{selectedItem, 5};

if necessaryTags == -1
    %user selected a separator
    set(hObject,'Value',handles.lastSelected(2));
else
    set(handles.yTagHandles,'Visible','off');
    set(handles.yTagHandles(1:necessaryTags),'Visible','on');
    handles.lastSelected(2) = selectedItem;
    
    %turn on cen/align
    switch cenOrAlign
        case 0 
            set([handles.cenHandles;handles.alignHandles],'Visible','off');
        case 1
            set(handles.cenHandles,'Visible','on');
            set(handles.alignHandles,'Visible','off');
        case 2
            set([handles.cenHandles;handles.alignHandles],'Visible','on');
    end

    guidata(hObject,handles);
end

% --- Executes on selection change in adgui_zLabel_PD.
function adgui_zLabel_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_zLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%activate necessary tag-PD's
selectedItem = get(hObject,'Value');

necessaryTags=handles.PD_data{selectedItem,2};
cenOrAlign = handles.PD_data{selectedItem, 5};

if necessaryTags == -1
    %user selected a separator
    set(hObject,'Value',handles.lastSelected(3));
else
    set(handles.zTagHandles,'Visible','off');
    set(handles.zTagHandles(1:necessaryTags),'Visible','on');
    handles.lastSelected(3) = selectedItem;
    
    %turn on cen/align
    switch cenOrAlign
        case 0 
            set([handles.cenHandles;handles.alignHandles],'Visible','off');
        case 1
            set(handles.cenHandles,'Visible','on');
            set(handles.alignHandles,'Visible','off');
        case 2
            set([handles.cenHandles;handles.alignHandles],'Visible','on');
    end

    guidata(hObject,handles);
end

% --- Executes on selection change in adgui_filename_PD.
function adgui_filename_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_filename_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pdString = get(hObject,'String');
currentSelection = get(hObject,'Value');

switch pdString{currentSelection}
    case 'no file loaded'
        %user clicked on wrong item - stop execution
        handles.lastFile = 1;
        guidata(hObject,handles);
        return
    case 'add file...'
        %load new movie
        adgui_loadData_CB(hObject,[],handles);
    otherwise
        %change PD menu selection (use currentSelection-1 since add file... is #1!)
        set(hObject,'ForegroundColor',extendedColors(handles.colorList(mod(currentSelection-1-1,size(handles.colorList,1))+1,:)),...
            'BackgroundColor','w');
        anaDat = handles.data(currentSelection-1).anaDat;
        set([handles.xTagHandles;handles.yTagHandles;handles.zTagHandles;handles.cenHandles(2);handles.alignHandles(2:3)],...
            'Value',1,'String',[{'set tag'};anaDat(1).info.labelColor]);
        handles.lastFile = currentSelection;
        guidata(hObject,handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = analyzeDataGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over figure background.
function adgui_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to adgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function adgui_xLabel_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_xLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function adgui_xTag1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_xTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_xTag1_PD.
function adgui_xTag1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_xTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_xTag1_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_xTag1_PD


% --- Executes during object creation, after setting all properties.
function adgui_xTag4_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_xTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_xTag4_PD.
function adgui_xTag4_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_xTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_xTag4_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_xTag4_PD


% --- Executes during object creation, after setting all properties.
function adgui_xTag3_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_xTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_xTag3_PD.
function adgui_xTag3_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_xTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_xTag3_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_xTag3_PD


% --- Executes during object creation, after setting all properties.
function adgui_xTag2_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_xTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_xTag2_PD.
function adgui_xTag2_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_xTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_xTag2_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_xTag2_PD


% --- Executes during object creation, after setting all properties.
function adgui_yLabel_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_yLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes during object creation, after setting all properties.
function adgui_yTag1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_yTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_yTag1_PD.
function adgui_yTag1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_yTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_yTag1_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_yTag1_PD


% --- Executes during object creation, after setting all properties.
function adgui_yTag4_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_yTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_yTag4_PD.
function adgui_yTag4_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_yTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_yTag4_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_yTag4_PD


% --- Executes during object creation, after setting all properties.
function adgui_yTag3_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_yTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_yTag3_PD.
function adgui_yTag3_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_yTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_yTag3_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_yTag3_PD


% --- Executes during object creation, after setting all properties.
function adgui_yTag2_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_yTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_yTag2_PD.
function adgui_yTag2_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_yTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_yTag2_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_yTag2_PD


% --- Executes during object creation, after setting all properties.
function adgui_zLabel_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_zLabel_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes during object creation, after setting all properties.
function adgui_zTag1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_zTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_zTag1_PD.
function adgui_zTag1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_zTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_zTag1_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_zTag1_PD


% --- Executes during object creation, after setting all properties.
function adgui_zTag4_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_zTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_zTag4_PD.
function adgui_zTag4_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_zTag4_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_zTag4_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_zTag4_PD


% --- Executes during object creation, after setting all properties.
function adgui_zTag3_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_zTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_zTag3_PD.
function adgui_zTag3_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_zTag3_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_zTag3_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_zTag3_PD


% --- Executes during object creation, after setting all properties.
function adgui_zTag2_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_zTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_zTag2_PD.
function adgui_zTag2_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_zTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_zTag2_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_zTag2_PD





% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.isFirst
    handles.timer = timer('StartDelay',5,'TimerFcn','clickFcn');
    start(handles.timer);
    handles.isFirst = 0;
    handles.clicks = 1;
    guidata(hObject,handles);
else
    handles.clicks = handles.clicks +1;
    guidata(hObject,handles);
end

    



% --- Executes on button press in adgui_showLegend_TB.
function adgui_showLegend_TB_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_showLegend_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Value')
    case 1 %button is toggled
        adguiLegend; %show legend
    case 0
        legendH = findall(0,'Tag','adguiLegend');
        close(legendH); %close legend
end

% --- Executes during object creation, after setting all properties.
function adgui_filename_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_filename_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes during object creation, after setting all properties.
function adgui_cenTag1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_cenTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_cenTag1_PD.
function adgui_cenTag1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_cenTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_cenTag1_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_cenTag1_PD


% --- Executes during object creation, after setting all properties.
function adgui_alignTag2_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_alignTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_alignTag2_PD.
function adgui_alignTag2_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_alignTag2_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_alignTag2_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_alignTag2_PD


% --- Executes during object creation, after setting all properties.
function adgui_alignTag1_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adgui_alignTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adgui_alignTag1_PD.
function adgui_alignTag1_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_alignTag1_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns adgui_alignTag1_PD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from adgui_alignTag1_PD


% --- Executes on button press in adgui_errorbar_check.
function adgui_errorbar_check_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_errorbar_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adgui_errorbar_check


% --------------------------------------------------------------------
function adgui_menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function adgui_menu_loadIdlist_Callback(hObject, eventdata, handles)
% hObject    handle to adgui_menu_loadIdlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


