function varargout = focalAdhesionSegmentationProcessGUI(varargin)
% focalAdhesionSegmentationProcessGUI M-file for focalAdhesionSegmentationProcessGUI.fig
%      focalAdhesionSegmentationProcessGUI, by itself, creates a new focalAdhesionSegmentationProcessGUI or raises the existing
%      singleton*.
%
%      H = focalAdhesionSegmentationProcessGUI returns the handle to a new focalAdhesionSegmentationProcessGUI or the handle to
%      the existing singleton*.
%
%      focalAdhesionSegmentationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in focalAdhesionSegmentationProcessGUI.M with the given input arguments.
%
%      focalAdhesionSegmentationProcessGUI('Property','Value',...) creates a new focalAdhesionSegmentationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before focalAdhesionSegmentationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to focalAdhesionSegmentationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help focalAdhesionSegmentationProcessGUI

% Last Modified by GUIDE v2.5 26-Apr-2013 08:08:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @focalAdhesionSegmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @focalAdhesionSegmentationProcessGUI_OutputFcn, ...
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


% --- Executes just before focalAdhesionSegmentationProcessGUI is made visible.
function focalAdhesionSegmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set default parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set-up parameters
userData.numParams = {'SteerableFilterSigma', 'OpeningRadiusXY', ...
                      'OpeningHeightT', 'MinVolTime'};
for i =1 : numel(userData.numParams)
    paramName = userData.numParams{i};
    set(handles.(['edit_' paramName]), 'String', funParams.(paramName));
end

% Update GUI user data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = focalAdhesionSegmentationProcessGUI_OutputFcn(~, ~, handles) 
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



function edit_SteerableFilterSigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SteerableFilterSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SteerableFilterSigma as text
%        str2double(get(hObject,'String')) returns contents of edit_SteerableFilterSigma as a double


% --- Executes during object creation, after setting all properties.
function edit_SteerableFilterSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SteerableFilterSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MinVolTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinVolTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinVolTime as text
%        str2double(get(hObject,'String')) returns contents of edit_MinVolTime as a double


% --- Executes during object creation, after setting all properties.
function edit_MinVolTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinVolTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_OpeningHeightT_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OpeningHeightT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OpeningHeightT as text
%        str2double(get(hObject,'String')) returns contents of edit_OpeningHeightT as a double


% --- Executes during object creation, after setting all properties.
function edit_OpeningHeightT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OpeningHeightT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_OpeningRadiusXY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OpeningRadiusXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OpeningRadiusXY as text
%        str2double(get(hObject,'String')) returns contents of edit_OpeningRadiusXY as a double


% --- Executes during object creation, after setting all properties.
function edit_OpeningRadiusXY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OpeningRadiusXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
