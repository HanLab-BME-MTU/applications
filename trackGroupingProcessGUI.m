function varargout = trackGroupingProcessGUI(varargin)
% trackGroupingProcessGUI M-file for trackGroupingProcessGUI.fig
%      trackGroupingProcessGUI, by itself, creates a new trackGroupingProcessGUI or raises the existing
%      singleton*.
%
%      H = trackGroupingProcessGUI returns the handle to a new trackGroupingProcessGUI or the handle to
%      the existing singleton*.
%
%      trackGroupingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in trackGroupingProcessGUI.M with the given input arguments.
%
%      trackGroupingProcessGUI('Property','Value',...) creates a new trackGroupingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackGroupingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackGroupingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackGroupingProcessGUI

% Last Modified by GUIDE v2.5 21-Mar-2013 17:45:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackGroupingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackGroupingProcessGUI_OutputFcn, ...
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


% --- Executes just before trackGroupingProcessGUI is made visible.
function trackGroupingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set default parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set-up parameters
userData.numParams = {'minLifetime', 'maxDistance', 'minOverlap', 'bandWidth',...
                      'minDistance','alpha'};
for i =1 : numel(userData.numParams)
    paramName = userData.numParams{i};
    set(handles.(['edit_' paramName]), 'String', funParams.(paramName));
end

% Update GUI user data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = trackGroupingProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

delete(userData.helpFig(ishandle(userData.helpFig))); 

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
funParams.ChannelIndex = get(handles.listbox_selectedChannels, 'Userdata');

% Get numerical parameters
userData = get(handles.figure1, 'UserData');

nParam = numel(userData.numParams);

for i = 1:nParam    
    paramName = userData.numParams{i};
    value = str2double(get(handles.(['edit_' paramName]),{'String'}));
    if isnan(value) || value < 0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' paramName]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(paramName)=value; 
end



processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



function edit_minLifetime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minLifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minLifetime as text
%        str2double(get(hObject,'String')) returns contents of edit_minLifetime as a double


% --- Executes during object creation, after setting all properties.
function edit_minLifetime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minLifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxDistance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxDistance as text
%        str2double(get(hObject,'String')) returns contents of edit_maxDistance as a double


% --- Executes during object creation, after setting all properties.
function edit_maxDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minDistance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minDistance as text
%        str2double(get(hObject,'String')) returns contents of edit_minDistance as a double


% --- Executes during object creation, after setting all properties.
function edit_minDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minOverlap as text
%        str2double(get(hObject,'String')) returns contents of edit_minOverlap as a double


% --- Executes during object creation, after setting all properties.
function edit_minOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bandWidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bandWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bandWidth as text
%        str2double(get(hObject,'String')) returns contents of edit_bandWidth as a double


% --- Executes during object creation, after setting all properties.
function edit_bandWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bandWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
