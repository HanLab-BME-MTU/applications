function varargout = detectionPreviewGUI(varargin)
% DETECTIONPREVIEWGUI M-file for detectionPreviewGUI.fig
%      DETECTIONPREVIEWGUI, by itself, creates a new DETECTIONPREVIEWGUI or raises the existing
%      singleton*.
%
%      H = DETECTIONPREVIEWGUI returns the handle to a new DETECTIONPREVIEWGUI or the handle to
%      the existing singleton*.
%
%      DETECTIONPREVIEWGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTIONPREVIEWGUI.M with the given input arguments.
%
%      DETECTIONPREVIEWGUI('Property','Value',...) creates a new DETECTIONPREVIEWGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detectionPreviewGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detectionPreviewGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detectionPreviewGUI

% Last Modified by GUIDE v2.5 09-Jun-2011 12:09:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detectionPreviewGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @detectionPreviewGUI_OutputFcn, ...
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


% --- Executes just before detectionPreviewGUI is made visible.
function detectionPreviewGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Useful tools:
% 
% userData.mainFig - handle of main figure
% userData.handles_main - the "handles" of main figure
%
% userData.useLastImg - 1: last frame 0: first frame
%
%



userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

userData.default_scales = [1 4];
userData.default_mulfac = 3; 

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.handles_main = guidata(userData.mainFig);
userData.useLastImg = 0;

userData.projData = userData.handles_main.projList;

set(handles.slider_2_1, 'Value', get(userData.handles_main.slider_1, 'Value'))
set(handles.slider_2_2, 'Value', get(userData.handles_main.slider_2, 'Value'))
set(handles.slider_2_3, 'Value', get(userData.handles_main.slider_3, 'Value'))

set(handles.edit_2_1, 'String', get(userData.handles_main.edit_detect_1, 'String'))
set(handles.edit_2_2, 'String', get(userData.handles_main.edit_detect_2, 'String'))
set(handles.edit_2_3, 'String', get(userData.handles_main.edit_detect_3, 'String'))

set(handles.uipanel_3, 'SelectionChangeFcn', @uipanel_3_SelectionChangeFcn);

userData.bitDepth = str2double(get(userData.handles_main.bitDepth,'String'));
if mod(userData.bitDepth,2)~=0, 
    errordlg('Invalid bit depth first');
    delete(handles.figure1);
end
plusTipTroubleShootDetect(userData.projData, userData.default_scales, ...
    userData.bitDepth, userData.useLastImg, userData.default_mulfac, 0, 1, hObject, handles.axes_1_1);

plusTipTroubleShootDetect(userData.projData, userData.default_scales, ...
    userData.bitDepth, userData.useLastImg, userData.default_mulfac, 0, 0, hObject, handles.axes_1_2);

scales = [get(handles.slider_2_1, 'Value') get(handles.slider_2_2, 'Value')];
mulfac = get(handles.slider_2_3, 'Value');

plusTipTroubleShootDetect(userData.projData, scales,...
    userData.bitDepth, userData.useLastImg, mulfac, 0, 1, hObject, handles.axes_2_1);
plusTipTroubleShootDetect(userData.projData, scales, ...
    userData.bitDepth, userData.useLastImg, mulfac, 0, 0, hObject, handles.axes_2_2);

% Set zoom callbacks
axesList={'axes_1_1' 'axes_1_2' 'axes_2_1' 'axes_2_2'};
axesHandles = cellfun(@(x) handles.(x),axesList);
linkaxes(axesHandles);

% Update handles structure
set(hObject, 'UserData', userData)
guidata(hObject, handles);

% UIWAIT makes detectionPreviewGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detectionPreviewGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit_1_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_1_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_1_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_1_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_1_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_1_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_1_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_1_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_1_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_1_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_1_3_Callback(hObject, eventdata, handles)


value = str2double(get(hObject, 'String'));

if isnan(value) || value < get(handles.slider_1_3, 'Min') || value > get(handles.slider_1_3, 'Max')
    
    set(hObject, 'String', num2str(get(handles.slider_1_3, 'Value')))
    return
    
end

if value == get(handles.slider_1_3, 'Value')
   return 
end

set(handles.slider_1_3, 'Value', value);

% Recalculate
userData = get(handles.figure1, 'UserData');
scales = userData.default_scales;
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, value, 0, 0, handles.figure1, handles.axes_1_2);



% --- Executes during object creation, after setting all properties.
function edit_1_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_1_3_Callback(hObject, eventdata, handles)
value = get(hObject, 'Value');

if value == str2double( get(handles.edit_1_3, 'String') )
    return
end

set(handles.edit_1_3, 'String', num2str(value))

% Recalculate
userData = get(handles.figure1, 'UserData');
scales = userData.default_scales;
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, value, 0, 0, handles.figure1, handles.axes_1_2);


% --- Executes during object creation, after setting all properties.
function slider_1_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_1_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');


set(userData.handles_main.slider_1, 'Value', get(handles.slider_2_1, 'Value'))
set(userData.handles_main.slider_2, 'Value', get(handles.slider_2_2, 'Value'))
set(userData.handles_main.slider_3, 'Value', get(handles.slider_2_3, 'Value'))

set(userData.handles_main.edit_detect_1, 'String', num2str(get(handles.slider_2_1, 'Value')))
set(userData.handles_main.edit_detect_2, 'String', num2str(get(handles.slider_2_2, 'Value')))
set(userData.handles_main.edit_detect_3, 'String', num2str(get(handles.slider_2_3, 'Value')))

delete(handles.figure1)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1)



function edit_2_1_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = str2double(get(hObject, 'String'));

if isnan(value) || value < get(handles.slider_2_1, 'Min') || value > get(handles.slider_2_1, 'Max')
    
    set(hObject, 'String', num2str(get(handles.slider_2_1, 'Value')))
    return
end

if value == get(handles.slider_2_1, 'Value')
   return 
end

if get(handles.slider_2_2, 'Value') <= value

    set(hObject, 'String', num2str(get(handles.slider_2_1, 'Value')))
    return
end

set(handles.slider_2_1, 'Value', value);

% Recalculate
scales = [value get(handles.slider_2_2, 'Value')];
mulfac = get( handles.slider_2_3, 'Value');


plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 1, handles.figure1, handles.axes_2_1);
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 0, handles.figure1, handles.axes_2_2);

% --- Executes during object creation, after setting all properties.
function edit_2_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_2_2_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = str2double(get(hObject, 'String'));

if isnan(value) || value < get(handles.slider_2_2, 'Min') || value > get(handles.slider_2_2, 'Max')
    
    set(hObject, 'String', num2str(get(handles.slider_2_2, 'Value')))
    return
end

if value == get(handles.slider_2_2, 'Value')
   return 
end

if value <= get(handles.slider_2_1, 'Value') 
    
    set(hObject, 'String', num2str(get(handles.slider_2_2, 'Value')))
    return
end

set(handles.slider_2_2, 'Value', value)

% Recalculate
scales = [get(handles.slider_2_1, 'Value') value];
mulfac = get( handles.slider_2_3, 'Value');

plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 1, handles.figure1, handles.axes_2_1);
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 0, handles.figure1, handles.axes_2_2);


% --- Executes during object creation, after setting all properties.
function edit_2_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_2_3_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = str2double(get(hObject, 'String'));

if isnan(value) || value < get(handles.slider_2_3, 'Min') || value > get(handles.slider_2_3, 'Max')
    
    set(hObject, 'String', num2str(get(handles.slider_2_3, 'Value')))
    return
    
end

if value == get(handles.slider_2_3, 'Value')
   return 
end

set(handles.slider_2_3, 'Value', value);

% Recalculate
scales = [get(handles.slider_2_1, 'Value') get(handles.slider_2_2, 'Value')];
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, value, 0, 0, handles.figure1, handles.axes_2_2);


% --- Executes during object creation, after setting all properties.
function edit_2_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_2_1_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = get(hObject, 'Value');

if value == str2double( get(handles.edit_2_1, 'String') )
    return
end

if get(handles.slider_2_2, 'Value') <= value

    set(hObject, 'Value', str2double(get(handles.edit_2_1, 'String')))
    return
end

set(handles.edit_2_1, 'String', num2str(value))

% Recalculate
scales = [value get(handles.slider_2_2, 'Value')];
mulfac = get( handles.slider_2_3, 'Value');


plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 1, handles.figure1, handles.axes_2_1);
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 0, handles.figure1, handles.axes_2_2);


% --- Executes during object creation, after setting all properties.
function slider_2_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_2_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_2_2_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = get(hObject, 'Value');

if value == str2double( get(handles.edit_2_2, 'String') )
    return
end

if value <= get(handles.slider_2_1, 'Value') 
    
    set(handles.slider_2_2, 'Value', str2double(get(handles.edit_2_2, 'String')))
    return
end

set(handles.edit_2_2, 'String', num2str(value))

% Recalculate
scales = [get(handles.slider_2_1, 'Value') value];
mulfac = get( handles.slider_2_3, 'Value');

plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 1, handles.figure1, handles.axes_2_1);
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 0, handles.figure1, handles.axes_2_2);


% --- Executes during object creation, after setting all properties.
function slider_2_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_2_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_2_3_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
value = get(hObject, 'Value');

if value == str2double( get(handles.edit_2_3, 'String') )
    return
end

set(handles.edit_2_3, 'String', num2str(value))

% Recalculate

scales = [get(handles.slider_2_1, 'Value') get(handles.slider_2_2, 'Value')];
plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, value, 0, 0, handles.figure1, handles.axes_2_2);


% --- Executes during object creation, after setting all properties.
function slider_2_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_2_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_both.
function checkbox_both_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_both


function uipanel_3_SelectionChangeFcn(hObject, eventdata)
% Call back function of ration button group uipanel_1
handles = guidata(hObject);
userData = get(handles.figure1, 'UserData');

if isequal(userData.useLastImg, get(handles.radiobutton_1, 'Value'))
    
    if get(handles.radiobutton_1, 'Value')
    
        userData.useLastImg = 0;
        
    else
        userData.useLastImg = 1;
    end
    

    plusTipTroubleShootDetect(userData.projData, userData.default_scales, userData.bitDepth, userData.useLastImg, userData.default_mulfac, 0, 1, handles.figure1, handles.axes_1_1);

    plusTipTroubleShootDetect(userData.projData, userData.default_scales, userData.bitDepth, userData.useLastImg, userData.default_mulfac, 0, 0, handles.figure1, handles.axes_1_2);

    scales = [get(handles.slider_2_1, 'Value') get(handles.slider_2_2, 'Value')];
    mulfac = get(handles.slider_2_3, 'Value');

    plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 1, handles.figure1, handles.axes_2_1);

    plusTipTroubleShootDetect(userData.projData, scales, userData.bitDepth, userData.useLastImg, mulfac, 0, 0, handles.figure1, handles.axes_2_2);
    
end

set(handles.figure1, 'UserData', userData)
