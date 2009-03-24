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
%      FSMDATAVIEWERSETUP('Property','Value',...) creates a new FSMDATAVIEWERSETUP or raises the
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

% Last Modified by GUIDE v2.5 23-Mar-2009 12:10:40

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


% --- Executes on selection change in listboxBackground.
function listboxBackground_Callback(hObject, eventdata, handles)
% hObject    handle to listboxBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxBackground contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxBackground

hFsmDataViewerSetup = get(hObject, 'Parent');

status = 'on';

if get(hObject, 'Value') == 2
    status = 'off';
end

h = findobj(hFsmDataViewerSetup, 'Tag', 'textGreenChannel');
set(h, 'Enable', status);
h = findobj(hFsmDataViewerSetup, 'Tag', 'textBlueChannel');
set(h, 'Enable', status);
h = findobj(hFsmDataViewerSetup, 'Tag', 'editRedChannel');
set(h, 'String', '');
h = findobj(hFsmDataViewerSetup, 'Tag', 'editGreenChannel');
set(h, 'Enable', status);
set(h, 'String', '');
h = findobj(hFsmDataViewerSetup, 'Tag', 'editBlueChannel');
set(h, 'Enable', status);
set(h, 'String', '');
h = findobj(hFsmDataViewerSetup, 'Tag', 'pushButtonGreenChannel');
set(h, 'Enable', status);
h = findobj(hFsmDataViewerSetup, 'Tag', 'pushButtonBlueChannel');
set(h, 'Enable', status);


% --- Executes during object creation, after setting all properties.
function listboxBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editRedChannel_Callback(hObject, eventdata, handles)
% hObject    handle to editRedChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRedChannel as text
%        str2double(get(hObject,'String')) returns contents of editRedChannel as a double


% --- Executes during object creation, after setting all properties.
function editRedChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRedChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonRedChannel.
function pushButtonRedChannel_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonRedChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if (SpeedMap | KineticsMap) is selected, we expect '.mat' files.
hFsmDataViewerSetup = get(hObject, 'Parent');
h = findobj(hFsmDataViewerSetup, 'Tag', 'listboxBackground');

if (get(h, 'Value') == 1)
    [fileName, directoryName] = uigetfile({'*.tif';'*jpg';'*.png'}, 'Select first image');
else
    [fileName, directoryName] = uigetfile('*.mat', 'Select first image data');
end

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');
    h = findobj(hFsmDataViewerSetup, 'Tag', 'editRedChannel');
    set(h, 'String', [directoryName fileName]);
end

function editGreenChannel_Callback(hObject, eventdata, handles)
% hObject    handle to editGreenChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGreenChannel as text
%        str2double(get(hObject,'String')) returns contents of editGreenChannel as a double


% --- Executes during object creation, after setting all properties.
function editGreenChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGreenChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonGreenChannel.
function pushButtonGreenChannel_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonGreenChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if (SpeedMap | KineticsMap) is selected, we expect '.mat' files.
hFsmDataViewerSetup = get(hObject, 'Parent');
h = findobj(hFsmDataViewerSetup, 'Tag', 'listboxBackground');

if (get(h, 'Value') == 1)
    [fileName, directoryName] = uigetfile({'*.tif';'*jpg';'*.png'}, 'Select first image');
else
    [fileName, directoryName] = uigetfile('*.mat', 'Select first image data');
end

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editGreenChannel');
    set(h, 'String', [directoryName fileName]);
end


function editBlueChannel_Callback(hObject, eventdata, handles)
% hObject    handle to editBlueChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBlueChannel as text
%        str2double(get(hObject,'String')) returns contents of editBlueChannel as a double


% --- Executes during object creation, after setting all properties.
function editBlueChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBlueChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonBlueChannel.
function pushButtonBlueChannel_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonBlueChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if (SpeedMap | KineticsMap) is selected, we expect '.mat' files.
hFsmDataViewerSetup = get(hObject, 'Parent');
h = findobj(hFsmDataViewerSetup, 'Tag', 'listboxBackground');

if (get(h, 'Value') == 1)
    [fileName, directoryName] = uigetfile({'*.tif';'*jpg';'*.png'}, 'Select first image');
else
    [fileName, directoryName] = uigetfile('*.mat', 'Select first image data');
end

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editBlueChannel');
    set(h, 'String', [directoryName fileName]);
end


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


% --- Executes on button press in pushButtonMask.
function pushButtonMask_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.tif'}, 'Select first image');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editMask');
    set(h, 'String', [directoryName fileName]);
end

function editLayer1_Callback(hObject, eventdata, handles)
% hObject    handle to editLayer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLayer1 as text
%        str2double(get(hObject,'String')) returns contents of editLayer1 as a double


% --- Executes during object creation, after setting all properties.
function editLayer1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLayer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonLayer1.
function pushButtonLayer1_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonLayer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.mat'}, 'Select first image data');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editLayer1');
    set(h, 'String', [directoryName fileName]);
end


function editLayer2_Callback(hObject, eventdata, handles)
% hObject    handle to editLayer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLayer2 as text
%        str2double(get(hObject,'String')) returns contents of editLayer2 as a double


% --- Executes during object creation, after setting all properties.
function editLayer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLayer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonLayer2.
function pushButtonLayer2_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonLayer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.mat'}, 'Select first image data');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editLayer2');
    set(h, 'String', [directoryName fileName]);
end


function editLayer3_Callback(hObject, eventdata, handles)
% hObject    handle to editLayer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLayer3 as text
%        str2double(get(hObject,'String')) returns contents of editLayer3 as a double


% --- Executes during object creation, after setting all properties.
function editLayer3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLayer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonLayer3.
function pushButtonLayer3_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonLayer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.mat'}, 'Select first image data');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editLayer3');
    set(h, 'String', [directoryName fileName]);
end


function editLayer4_Callback(hObject, eventdata, handles)
% hObject    handle to editLayer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLayer4 as text
%        str2double(get(hObject,'String')) returns contents of editLayer4 as a double


% --- Executes during object creation, after setting all properties.
function editLayer4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLayer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonLayer4.
function pushButtonLayer4_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonLayer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.mat'}, 'Select first image data');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editLayer4');
    set(h, 'String', [directoryName fileName]);
end


function editLayer5_Callback(hObject, eventdata, handles)
% hObject    handle to editLayer5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLayer5 as text
%        str2double(get(hObject,'String')) returns contents of editLayer5 as a double


% --- Executes during object creation, after setting all properties.
function editLayer5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLayer5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushButtonLayer5.
function pushButtonLayer5_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonLayer5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, directoryName] = uigetfile({'*.mat'}, 'Select first image data');

if ischar(fileName) && ischar(directoryName)
    hFsmDataViewerSetup = get(hObject, 'Parent');

    h = findobj(hFsmDataViewerSetup, 'Tag', 'editLayer5');
    set(h, 'String', [directoryName fileName]);
end

% --- Executes on button press in pushButtonColorLayer1.
function pushButtonColorLayer1_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonColorLayer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentColor = get(hObject, 'backgroundColor');

color = uisetcolor(currentColor);

if numel(color) == 3
    set(hObject, 'backgroundColor', color);
end

% --- Executes on button press in pushButtonColorLayer2.
function pushButtonColorLayer2_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonColorLayer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentColor = get(hObject, 'backgroundColor');

color = uisetcolor(currentColor);

if numel(color) == 3
    set(hObject, 'backgroundColor', color);
end

% --- Executes on button press in pushButtonColorLayer3.
function pushButtonColorLayer3_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonColorLayer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentColor = get(hObject, 'backgroundColor');

color = uisetcolor(currentColor);

if numel(color) == 3
    set(hObject, 'backgroundColor', color);
end

% --- Executes on button press in pushButtonColorLayer4.
function pushButtonColorLayer4_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonColorLayer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentColor = get(hObject, 'backgroundColor');

color = uisetcolor(currentColor);

if numel(color) == 3
    set(hObject, 'backgroundColor', color);
end

% --- Executes on button press in pushButtonColorLayer5.
function pushButtonColorLayer5_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonColorLayer5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentColor = get(hObject, 'backgroundColor');

color = uisetcolor(currentColor);

if numel(color) == 3
    set(hObject, 'backgroundColor', color);
end


% --- Executes on button press in pushButtonOK.
function pushButtonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFsmDataViewerSetup = get(hObject, 'Parent');

[settings status] = getFsmDataViewerSettings(hFsmDataViewerSetup);

if (status)
    set(hFsmDataViewerSetup, 'UserData', settings);
    
    uiresume(hFsmDataViewerSetup);
end

% --- Executes on button press in pushButtonCancel.
function pushButtonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFsmDataViewerSetup = get(hObject, 'Parent');

set(hFsmDataViewerSetup, 'UserData', []);

uiresume(hFsmDataViewerSetup);