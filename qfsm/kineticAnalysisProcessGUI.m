function varargout = kineticAnalysisProcessGUI(varargin)
% kineticAnalysisProcessGUI M-file for kineticAnalysisProcessGUI.fig
%      kineticAnalysisProcessGUI, by itself, creates a new kineticAnalysisProcessGUI or raises the existing
%      singleton*.
%
%      H = kineticAnalysisProcessGUI returns the handle to a new kineticAnalysisProcessGUI or the handle to
%      the existing singleton*.
%
%      kineticAnalysisProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in kineticAnalysisProcessGUI.M with the given input arguments.
%
%      kineticAnalysisProcessGUI('Property','Value',...) creates a new kineticAnalysisProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kineticAnalysisProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kineticAnalysisProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kineticAnalysisProcessGUI

% Last Modified by GUIDE v2.5 19-Jul-2011 09:02:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kineticAnalysisProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @kineticAnalysisProcessGUI_OutputFcn, ...
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


% --- Executes just before kineticAnalysisProcessGUI is made visible.
function kineticAnalysisProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% ---------------------- Parameter Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

try 
    set(handles.edit_alpha,'String',userData.crtPackage.processes_{4}.funParams_.alpha);
catch ME
    set(handles.edit_alpha,'String','NA');
end

set(handles.popupmenu_bleachRed,'Value',...
    floor(funParams.bleachRed/7.2500e-05)+1);
set(handles.edit_sigma,'String',funParams.sigma);
set(handles.edit_timeWindow,'String',funParams.timeWindow);

% ----------------------------------------------------------------

% Choose default command line output for kineticAnalysisProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = kineticAnalysisProcessGUI_OutputFcn(~, ~, handles) 
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

if isfield(userData, 'previewFig') && ishandle(userData.previewFig)
   delete(userData.previewFig) 
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
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% Check user input

if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
timeWindow = str2double(get(handles.edit_timeWindow,'String'));
if isnan(timeWindow) || timeWindow<=0 || timeWindow > userData.MD.nFrames_ || mod(timeWindow,2)==0
    errordlg(sprintf('The time window should be an odd number between 0 and %g.',...
        userData.MD.nFrames_),'Setting Error','modal')
    return;
end


% Process Sanity check ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message 'Please double check your data.'],'Setting Error','modal');
    return;
end

% Set parameter 
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
funParams.timeWindow = str2double(get(handles.edit_timeWindow,'String'));
funParams.sigma = str2double(get(handles.edit_sigma,'String'));
funParams.bleachRed = 7.2500e-05*(get(handles.popupmenu_bleachRed,'Value')-1);


% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
