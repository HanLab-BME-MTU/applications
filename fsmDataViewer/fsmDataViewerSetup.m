function varargout = fsmDataViewerSetup(varargin)
% FSMDATAVIEWERSETUP M-file for fsmDataViewerSetup.fig
%      FSMDATAVIEWERSETUP, by itself, creates a new FSMDATAVIEWERSETUP or raises the existing
%      singleton*.
%
%      H = FSMDATAVIEWERSETUP returns the handle to a new FSMDATAVIEWERSETUP or the handle to
%      the existing singleton*.
%
%      FSMDATAVIEWERSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMDATAVIEWERSETUP.M with the given input arguments.
%
%      FSMDATAVIEWERSETUP('Property','Value',...) creates a new
%      FSMDATAVIEWERSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmDataViewerSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmDataViewerSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmDataViewerSetup

% Last Modified by GUIDE v2.5 14-Apr-2009 15:47:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmDataViewerSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmDataViewerSetup_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fsmDataViewerSetup is made visible.
function fsmDataViewerSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmDataViewerSetup (see VARARGIN)

% Choose default command line output for fsmDataViewerSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fsmDataViewer wait for user response (see UIRESUME)
uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = fsmDataViewerSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargout >= 1
    if isempty(handles)
        settings = '';
    else
        settings = get(hObject, 'UserData');

        delete(hObject);
    end
    varargout{1} = settings;
end

function editRootDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to editRootDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRootDirectory as text
%        str2double(get(hObject,'String')) returns contents of editRootDirectory as a double


% --- Executes during object creation, after setting all properties.
function editRootDirectory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRootDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', pwd);

% --- Executes on button press in pushbuttonRootDirectory.
function pushbuttonRootDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRootDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = get(hObject, 'Parent');

rootDirectory = uigetdir('', 'Select a root directory:');

if rootDirectory
    h = findobj(hFig, 'Tag', 'editRootDirectory');
    set(h, 'String', rootDirectory);
end

% --- Executes on selection change in listboxChannelType.
function listboxChannelType_Callback(hObject, eventdata, handles)
% hObject    handle to listboxChannelType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxChannelType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxChannelType

channelType = get(hObject, 'Value');

status = 'on';
if channelType == 1
    status = 'off';
end
    
hFig = get(hObject, 'Parent');
hPushButtonChannel = findobj(hFig, 'Tag', 'pushbuttonChannel');
set(hPushButtonChannel, 'Enable', status);
 
% --- Executes during object creation, after setting all properties.
function listboxChannelType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxChannelType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

channelPlugins = getPlugins();

set(hObject, 'String', {'Choose Channel Type...', channelPlugins(:).desc});

% --- Executes on button press in pushbuttonChannel.
function pushbuttonChannel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = get(hObject, 'Parent');
hFig = get(hPanel, 'Parent');

% Get the plugins list
channelPlugins = getPlugins();

% Get the root directory
h = findobj(hFig, 'Tag', 'editRootDirectory');
rootDirectory = get(h, 'String');
currentDirectory = pwd;

% Change directory
if ~isempty(rootDirectory)
    cd(rootDirectory);
end

% Get the channel type ID
h = findobj(hPanel, 'Tag', 'listboxChannelType');
channelTypeID = get(h, 'Value') - 1;

filterSpec = channelPlugins(channelTypeID).filterSpec;

% Get the image file
[fileName, directoryName] = uigetfile(filterSpec, 'Select an image file');

if ischar(fileName) && ischar(directoryName)
    h = findobj(hPanel, 'Tag', 'uitableChannels');
    data = get(h, 'Data');
    columnFormat = get(h, 'ColumnFormat');
    colorNames = columnFormat{3};
    newData = {true,...
        channelPlugins(channelTypeID).desc,...
        colorNames{1},...
        [directoryName fileName]};
    data = vertcat(data, newData);
    set(h, 'Data', data);
    set(h, 'ColumnWidth', {20, 110, 65, 280});
end

% Go back to the current directory
if ~isempty(rootDirectory)
    cd(currentDirectory);
end

% --- Executes when selected cell(s) is changed in uitableChannels.
function uitableChannels_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitableChannels (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uitableChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitableChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.

set(hObject, 'Data', cell(0, 4));

% --- Executes on button press in checkboxMask.
function checkboxMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMask
enable = get(hObject, 'Value');

hPanel = get(hObject, 'Parent');
hEditMask = findobj(hPanel, 'Tag', 'editMask');
hPushbuttonMask = findobj(hPanel, 'Tag', 'pushbuttonMask');

status = 'off';

if enable
    status = 'on';
end

set(hEditMask, 'Enable', status);
set(hPushbuttonMask, 'Enable', status);

function editMask_Callback(hObject, eventdata, handles)
% hObject    handle to editMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMask as text
%        str2double(get(hObject,'String')) returns contents of editMask as a double


% --- Executes during object creation, after setting all properties.
function editMask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonMask.
function pushbuttonMask_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = get(hObject, 'Parent');
hFig = get(hPanel, 'Parent');

% Get the root directory
h = findobj(hFig, 'Tag', 'editRootDirectory');
rootDirectory = get(h, 'String');
currentDirectory = pwd;

% Change directory
if ~isempty(rootDirectory)
    cd(rootDirectory);
end

[fileName, directoryName] = uigetfile({'*.tif'}, 'Select a cell mask');

if ischar(fileName) && ischar(directoryName)
    h = findobj(hPanel, 'Tag', 'editMask');
    set(h, 'String', [directoryName fileName]);
end

% Go back to the current directory
if ~isempty(rootDirectory)
    cd(currentDirectory);
end


% --- Executes on selection change in listboxLayerType.
function listboxLayerType_Callback(hObject, eventdata, handles)
% hObject    handle to listboxLayerType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxLayerType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxLayerType

channelType = get(hObject, 'Value');

status = 'on';
if channelType == 1
    status = 'off';
end
    
hFig = get(hObject, 'Parent');
hPushButtonLayer = findobj(hFig, 'Tag', 'pushbuttonLayer');
set(hPushButtonLayer, 'Enable', status);

% --- Executes during object creation, after setting all properties.
function listboxLayerType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxLayerType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

[channelPlugins, layerPlugins] = getPlugins();

set(hObject, 'String', {'Choose Layer Type...', layerPlugins(:).desc});

% --- Executes on button press in pushbuttonLayer.
function pushbuttonLayer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = get(hObject, 'Parent');
hFig = get(hPanel, 'Parent');

% Get the plugins list
[channelPlugins layerPlugins] = getPlugins();

% Get the root directory
h = findobj(hFig, 'Tag', 'editRootDirectory');
rootDirectory = get(h, 'String');
currentDirectory = pwd;

% Change directory
if ~isempty(rootDirectory)
    cd(rootDirectory);
end

% Get the layer type ID
h = findobj(hPanel, 'Tag', 'listboxLayerType');
layerTypeID = get(h, 'Value') - 1;

filterSpec = layerPlugins(layerTypeID).filterSpec;

% Get the image file
[fileName, directoryName] = uigetfile(filterSpec, 'Select an image file');

if ischar(fileName) && ischar(directoryName)
    h = findobj(hPanel, 'Tag', 'uitableLayers');
    data = get(h, 'Data');
    columnFormat = get(h, 'ColumnFormat');
    colorNames = columnFormat{3};
    newData = {true,...
        layerPlugins(layerTypeID).desc,...
        colorNames{1},...
        [directoryName fileName]};
    data = vertcat(data, newData);
    set(h, 'Data', data);
    set(h, 'ColumnWidth', {20, 110, 65, 280});
end

% Go back to the current directory
if ~isempty(rootDirectory)
    cd(currentDirectory);
end

% --- Executes when selected cell(s) is changed in uitableLayers.
function uitableLayers_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitableLayers (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function uitableLayers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitableLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Data', cell(0, 4));

% --- Executes on button press in pushbuttonOK.
function pushbuttonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = get(hObject, 'Parent');

[settings status] = getSettings(hFig);

if (status)
    set(hFig, 'UserData', settings);
    
    uiresume(hFig);
end

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = get(hObject, 'Parent');

set(hFig, 'UserData', []);

uiresume(hFig);

