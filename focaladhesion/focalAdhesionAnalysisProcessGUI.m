function varargout = focalAdhesionAnalysisProcessGUI(varargin)
% focalAdhesionAnalysisProcessGUI M-file for focalAdhesionAnalysisProcessGUI.fig
%      focalAdhesionAnalysisProcessGUI, by itself, creates a new focalAdhesionAnalysisProcessGUI or raises the existing
%      singleton*.
%
%      H = focalAdhesionAnalysisProcessGUI returns the handle to a new focalAdhesionAnalysisProcessGUI or the handle to
%      the existing singleton*.
%
%      focalAdhesionAnalysisProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in focalAdhesionAnalysisProcessGUI.M with the given input arguments.
%
%      focalAdhesionAnalysisProcessGUI('Property','Value',...) creates a new focalAdhesionAnalysisProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before focalAdhesionAnalysisProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to focalAdhesionAnalysisProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help focalAdhesionAnalysisProcessGUI

% Last Modified by GUIDE v2.5 28-Feb-2017 13:41:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @focalAdhesionAnalysisProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @focalAdhesionAnalysisProcessGUI_OutputFcn, ...
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


% --- Executes just before focalAdhesionAnalysisProcessGUI is made visible.
function focalAdhesionAnalysisProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set default parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set-up parameters
userData.numParams = {'bandwidthNA', 'minLifetime', 'minFALengthMicron'};
userData.checkBoxes = {'reTrack','onlyEdge','matchWithFA',...
                       'getEdgeRelatedFeatures'};
                        
% Set edit strings/numbers
for paramName = userData.numParams
    set(handles.(['edit_' paramName{1}]), 'String', funParams.(paramName{1}));
end

% Set edit strings/numbers
for paramName = userData.checkBoxes
    set(handles.(['checkbox_' paramName{1}]), 'Value', funParams.(paramName{1}));
end

% Update GUI user data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = focalAdhesionAnalysisProcessGUI_OutputFcn(~, ~, handles) 
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

if ~isempty(userData)
    delete(userData.helpFig(ishandle(userData.helpFig))); 
    set(handles.figure1, 'UserData', userData);
    guidata(hObject,handles);
end


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

% Get number params
for paramName = userData.numParams    
    value = str2double(get(handles.(['edit_' paramName{1}]),{'String'}));
    if isnan(value) || value < 0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' paramName{1}]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(paramName{1})=value; 
end

% Get checkbox params
for paramName = userData.checkBoxes
    funParams.(paramName{1}) = logical(get(handles.(['checkbox_' paramName{1}]),'Value')); 
end

processGUI_ApplyFcn(hObject, eventdata, handles, funParams);



% --- Executes on button press in checkbox_reTrack.
function checkbox_reTrack_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_reTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_reTrack



% --- Executes on button press in checkbox_onlyEdge.
function checkbox_onlyEdge_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_onlyEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_onlyEdge


% --- Executes on button press in checkbox_getEdgeRelatedFeatures.
function checkbox_getEdgeRelatedFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_getEdgeRelatedFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_getEdgeRelatedFeatures


% --- Executes on button press in checkbox_matchWithFA.
function checkbox_matchWithFA_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_matchWithFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_matchWithFA
