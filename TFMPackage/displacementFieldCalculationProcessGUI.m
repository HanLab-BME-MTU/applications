function varargout = displacementFieldCalculationProcessGUI(varargin)
% displacementFieldCalculationProcessGUI M-file for displacementFieldCalculationProcessGUI.fig
%      displacementFieldCalculationProcessGUI, by itself, creates a new displacementFieldCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = displacementFieldCalculationProcessGUI returns the handle to a new displacementFieldCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      displacementFieldCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in displacementFieldCalculationProcessGUI.M with the given input arguments.
%
%      displacementFieldCalculationProcessGUI('Property','Value',...) creates a new displacementFieldCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before displacementFieldCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to displacementFieldCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help displacementFieldCalculationProcessGUI

% Last Modified by GUIDE v2.5 17-Mar-2013 14:46:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @displacementFieldCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @displacementFieldCalculationProcessGUI_OutputFcn, ...
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


% --- Executes just before displacementFieldCalculationProcessGUI is made visible.
function displacementFieldCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)


processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
set(handles.edit_referenceFramePath,'String',funParams.referenceFramePath);
userData.numParams ={'alpha','minCorLength', 'maxFlowSpeed'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams);
set(handles.edit_maxFlowSpeedNmMin,'String',...
    funParams.maxFlowSpeed*userData.MD.pixelSize_/userData.MD.timeInterval_*60);
set(handles.checkbox_highRes, 'Value', funParams.highRes);
set(handles.checkbox_mode, 'Value', strcmp(funParams.mode, 'accurate'));

% Propagate stage drift correction parameters if no process and stage drift
% correction parameters has been set up
stageDriftCorrProc = userData.crtPackage.processes_{1};
if ~isempty(stageDriftCorrProc)
    set(handles.edit_referenceFramePath,'Enable','off');
    set(handles.pushbutton_selectReferenceFrame,'Enable','off');
    if isempty(userData.crtPackage.processes_{userData.procID}) 
        set(handles.edit_alpha,'String', stageDriftCorrProc.funParams_.alpha);
    end
end

% Override default channels callback function
set(handles.checkbox_all,'Callback',@(hObject,eventdata)...
    checkallChannels(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_select,'Callback',@(hObject,eventdata)...
    selectChannel(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_delete,'Callback',@(hObject,eventdata)...
    deleteChannel(hObject,eventdata,guidata(hObject)));

% Choose default command line output for displacementFieldCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% Update value of psf sigma
update_psfSigma(handles);


% --- Outputs from this function are returned to the command line.
function varargout = displacementFieldCalculationProcessGUI_OutputFcn(~, ~, handles) 
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
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
else
    channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
    funParams.ChannelIndex = channelIndex;
end

% Retrieve reference frame path
if strcmpi(get(handles.edit_referenceFramePath,'Enable'),'on')
    funParams.referenceFramePath=get(handles.edit_referenceFramePath,'String');
    if isempty(funParams.referenceFramePath)
        errordlg('Please select a reference frame.','Setting Error','modal')
        return;
    end
end

% Read numeric information
for i = 1:numel(userData.numParams),
    value = get(handles.(['edit_' userData.numParams{i}]),'String');
    if isempty(value)
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=str2double(value); 
end

funParams.highRes = get(handles.checkbox_highRes, 'Value');
if get(handles.checkbox_mode, 'Value'),
    funParams.mode = 'accurate';
else
    funParams.mode = 'fast';
end


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


% --- Executes on button press in pushbutton_selectReferenceFrame.
function pushbutton_selectReferenceFrame_Callback(hObject, eventdata, handles)
[file path]=uigetfile({'*.tif;*.TIF;*.stk;*.STK;*.bmp;*.BMP;*.jpg;*.JPG',...
    'Image files (*.tif,*.stk,*.bmp,*.jpg)'},...
    'Select the reference frame');
if ~isequal(file,0) && ~isequal(path,0)
    set(handles.edit_referenceFramePath,'String',[path file]);
end


function edit_maxFlowSpeed_Callback(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
value=str2double(get(handles.edit_maxFlowSpeed,'String'));
set(handles.edit_maxFlowSpeedNmMin,'String',...
    value*userData.MD.pixelSize_/userData.MD.timeInterval_*60);

function selectChannel(hObject, eventdata, handles)

selectChannel_Callback(hObject, eventdata, handles);
update_psfSigma(handles);


function deleteChannel(hObject, eventdata, handles)

deleteChannel_Callback(hObject, eventdata, handles);
update_psfSigma(handles);

function checkallChannels(hObject, eventdata, handles)

checkallChannels_Callback(hObject, eventdata, handles);
update_psfSigma(handles);

function update_psfSigma(handles)

userData = get(handles.figure1, 'UserData');
selectedChannels = get(handles.listbox_selectedChannels,'UserData');
if ~isempty(selectedChannels)
    psfSigma = userData.MD.channels_(selectedChannels(1)).psfSigma_;
else
    psfSigma = '';
end
set(handles.edit_psfSigma, 'String', psfSigma);
