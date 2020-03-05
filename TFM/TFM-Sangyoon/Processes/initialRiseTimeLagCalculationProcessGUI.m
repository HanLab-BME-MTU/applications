function varargout = initialRiseTimeLagCalculationProcessGUI(varargin)
% initialRiseTimeLagCalculationProcessGUI M-file for initialRiseTimeLagCalculationProcessGUI.fig
%      initialRiseTimeLagCalculationProcessGUI, by itself, creates a new initialRiseTimeLagCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = initialRiseTimeLagCalculationProcessGUI returns the handle to a new initialRiseTimeLagCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      initialRiseTimeLagCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in initialRiseTimeLagCalculationProcessGUI.M with the given input arguments.
%
%      initialRiseTimeLagCalculationProcessGUI('Property','Value',...) creates a new initialRiseTimeLagCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before initialRiseTimeLagCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to initialRiseTimeLagCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help initialRiseTimeLagCalculationProcessGUI

% Last Modified by GUIDE v2.5 02-Nov-2017 14:19:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @initialRiseTimeLagCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @initialRiseTimeLagCalculationProcessGUI_OutputFcn, ...
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


% --- Executes just before initialRiseTimeLagCalculationProcessGUI is made visible.
function initialRiseTimeLagCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

set(handles.edit_mainSlave,'String',num2str(funParams.mainSlave));

% Choose default command line output for initialRiseTimeLagCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = initialRiseTimeLagCalculationProcessGUI_OutputFcn(~, ~, handles) 
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

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

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
userData = get(handles.figure1, 'UserData');

funParams.mainSlave=str2double(get(handles.edit_mainSlave,'String'));

% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_saveTractionField.
function checkbox_saveTractionField_Callback(hObject, eventdata, handles)
%
%
%


% --- Executes on button press in checkbox_fill.
function checkbox_fill_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fill



function edit_mainSlave_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mainSlave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mainSlave as text
%        str2double(get(hObject,'String')) returns contents of edit_mainSlave as a double


% --- Executes during object creation, after setting all properties.
function edit_mainSlave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mainSlave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
