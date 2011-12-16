function varargout = flowTrackingProcessGUI(varargin)
% flowTrackingProcessGUI M-file for flowTrackingProcessGUI.fig
%      flowTrackingProcessGUI, by itself, creates a new flowTrackingProcessGUI or raises the existing
%      singleton*.
%
%      H = flowTrackingProcessGUI returns the handle to a new flowTrackingProcessGUI or the handle to
%      the existing singleton*.
%
%      flowTrackingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in flowTrackingProcessGUI.M with the given input arguments.
%
%      flowTrackingProcessGUI('Property','Value',...) creates a new flowTrackingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flowTrackingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flowTrackingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flowTrackingProcessGUI

% Last Modified by GUIDE v2.5 29-Nov-2011 10:45:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flowTrackingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @flowTrackingProcessGUI_OutputFcn, ...
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


% --- Executes just before flowTrackingProcessGUI is made visible.
function flowTrackingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Choose default command line output for noiseEstimationProcessGUI
handles.output = hObject;

% Parameter Setup 
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

userData.numParams = {'firstImage','lastImage','timeWindow','timeStepSize',...
  'minCorLength','maxCorLength','minFeatureSize',...
  'edgeErodeWidth','maxFlowSpeed'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),userData.numParams);
set(handles.edit_maxFlowSpeedNmMin,'String',...
    funParams.maxFlowSpeed*userData.MD.pixelSize_/userData.MD.timeInterval_*60);

% Stationary background parameters
substractBackground = (funParams.numStBgForAvg~=0);
set(handles.checkbox_substractBackground,'Value',substractBackground)
if substractBackground
    if funParams.numStBgForAvg==-1
        set(handles.checkbox_useAllStBgImages,'Value',1);
        set(handles.edit_numStBgForAvg,'String','','Enable','off');
    else
        set(handles.checkbox_useAllStBgImages,'Value',0);
        set(handles.edit_numStBgForAvg,'String',funParams.numStBgForAvg,...
            'Enable','on');
    end
else
    set(handles.checkbox_useAllStBgImages,'Value',0,'Enable','off');
    set(handles.edit_numStBgForAvg,'String','','Enable','off');
end

% Outlier parameters
detectOutliers = ~isempty(funParams.outlierThreshold);
set(handles.checkbox_filterOutliers,'Value',detectOutliers);
if detectOutliers
    set(handles.edit_outlierThreshold,'String',funParams.outlierThreshold,...
        'Enable','on');
else
    set(handles.edit_outlierThreshold,'String','','Enable','off');
end


% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = flowTrackingProcessGUI_OutputFcn(~, ~, handles) 
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

function pushbutton_done_Callback(hObject, eventdata, handles)

% Input check
userData = get(handles.figure1, 'UserData');
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

% Process Sanity check  ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Retrieve GUI-defined parameters
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

% Retrieve numeric parameters
for i = 1:numel(userData.numParams),
    value = str2double(get(handles.(['edit_' userData.numParams{i}]),'String'));
    if isnan(value) || value<0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=value;
end

% Retrieve stationary background parameters
if ~get(handles.checkbox_substractBackground,'Value')
    funParams.numStBgForAvg=0;
else
    if get(handles.checkbox_useAllStBgImages,'Value')
        funParams.numStBgForAvg=0-1;
    else
        numStBgForAvg= str2double(get(handles.edit_numStBgForAvg,'String'));
        if isnan(numStBgForAvg) || numStBgForAvg<=0
            errordlg(['Please enter a valid value for the '...
                get(handles.text_numStBgForAvg,'String') '.'],'Setting Error','modal');
            return;
        end
        funParams.numStBgForAvg=numStBgForAvg;
    end
end

% Retrieve outlier parameters
if get(handles.checkbox_filterOutliers,'Value')
    funParams.outlierThreshold=str2double(get(handles.edit_outlierThreshold,...
        'String'));
else
    funParams.outlierThreshold=[];
end

% Set parameters and update main window
setLastImage =@(x) parseProcessParams(x,struct('lastImage',...
    min(x.funParams_.lastImageNum,x.owner_.nFrames_)));
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on button press in checkbox_filterOutliers.
function checkbox_filterOutliers_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.edit_outlierThreshold,'Enable','on')
else
    set(handles.edit_outlierThreshold,'Enable','off')
end


function edit_maxFlowSpeed_Callback(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
value=str2double(get(handles.edit_maxFlowSpeed,'String'));
set(handles.edit_maxFlowSpeedNmMin,'String',...
    value*userData.MD.pixelSize_/userData.MD.timeInterval_*60);


% --- Executes on button press in checkbox_substractBackground.
function checkbox_substractBackground_Callback(hObject, eventdata, handles)

if get(hObject,'Value'), 
    set(handles.checkbox_useAllStBgImages,'Enable','on');
    set(handles.edit_numStBgForAvg,'Enable','on');
else
    set(handles.checkbox_useAllStBgImages,'Value',0,'Enable','off');
    set(handles.edit_numStBgForAvg,'String','','Enable','off');
end

% --- Executes on button press in checkbox_useAllStBgImages.
function checkbox_useAllStBgImages_Callback(hObject, eventdata, handles)

if get(hObject,'Value'), 
    set(handles.edit_numStBgForAvg,'String','','Enable','off');
else
    set(handles.edit_numStBgForAvg,'Enable','on');
end
