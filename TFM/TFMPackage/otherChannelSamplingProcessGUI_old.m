function varargout = otherChannelSamplingProcessGUI(varargin)
% otherChannelSamplingProcessGUI M-file for otherChannelSamplingProcessGUI.fig
%      otherChannelSamplingProcessGUI, by itself, creates a new otherChannelSamplingProcessGUI or raises the existing
%      singleton*.
%
%      H = otherChannelSamplingProcessGUI returns the handle to a new otherChannelSamplingProcessGUI or the handle to
%      the existing singleton*.
%
%      otherChannelSamplingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in otherChannelSamplingProcessGUI.M with the given input arguments.
%
%      otherChannelSamplingProcessGUI('Property','Value',...) creates a new otherChannelSamplingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before otherChannelSamplingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to otherChannelSamplingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help otherChannelSamplingProcessGUI

% Last Modified by GUIDE v2.5 01-Mar-2026 11:28:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @otherChannelSamplingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @otherChannelSamplingProcessGUI_OutputFcn, ...
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


% --- Executes just before otherChannelSamplingProcessGUI is made visible.
function otherChannelSamplingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

userData = get(handles.figure1, 'UserData');

% Propagate stage drift correction parameters if no process AND stage drift
% correction parameters has been set up
stageDriftCorrProc = userData.crtPackage.processes_{1};

if ~isempty(stageDriftCorrProc) 
    set(handles.edit_referenceFramePath, 'String',stageDriftCorrProc.funParams_.referenceFramePath)
    set(handles.edit_referenceFramePath,'Enable','off');
    set(handles.pushbutton_selectReferenceFrame,'Enable','off');
end
    
if ~isempty(stageDriftCorrProc) && strcmp(stageDriftCorrProc.name_, 'Bead Tracking Drift Correction')
    if isempty(userData.crtPackage.processes_{userData.procID}) 
        set(handles.edit_alpha,'String', stageDriftCorrProc.funParams_.alpha);
    end
end


if ~isempty(stageDriftCorrProc) 
    channelString = userData.MD.getChannelPaths(stageDriftCorrProc.funParams_.iBeadChannel);

    set(handles.listbox_selectedChannels,'String',channelString,...
        'UserData',stageDriftCorrProc.funParams_.iBeadChannel);
end

% Set process parameters
funParams = userData.crtProc.funParams_;
set(handles.edit_referenceFramePath,'String',funParams.referenceFramePath);
userData.numParams ={'alpha','minCorLength', 'maxFlowSpeed'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams);
if ~isempty(userData.MD.timeInterval_)
    set(handles.edit_maxFlowSpeedNmMin,'String',...
        funParams.maxFlowSpeed*userData.MD.pixelSize_/userData.MD.timeInterval_*60);
else
    set(handles.edit_maxFlowSpeedNmMin,'String',...
        funParams.maxFlowSpeed*userData.MD.pixelSize_);
end

% Override default channels callback function
set(handles.checkbox_all,'Callback',@(hObject,eventdata)...
    checkallChannels(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_select,'Callback',@(hObject,eventdata)...
    selectChannel(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_delete,'Callback',@(hObject,eventdata)...
    deleteChannel(hObject,eventdata,guidata(hObject)));

% Choose default command line output for otherChannelSamplingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% Update value of psf sigma
update_psfSigma(handles);


% --- Outputs from this function are returned to the command line.
function varargout = otherChannelSamplingProcessGUI_OutputFcn(~, ~, handles) 
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

if get(handles.checkbox_mode, 'Value'),
    funParams.mode = 'accurate';
else
    funParams.mode = 'fast';
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end



function selectChannel(hObject, eventdata, handles)

selectChannel_Callback(hObject, eventdata, handles);
update_psfSigma(handles);


function deleteChannel(hObject, eventdata, handles)

deleteChannel_Callback(hObject, eventdata, handles);
update_psfSigma(handles);

function checkallChannels(hObject, eventdata, handles)

checkallChannels_Callback(hObject, eventdata, handles);
update_psfSigma(handles);


%--- Executes during object creation, after setting all properties.
function edit_sigCrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

