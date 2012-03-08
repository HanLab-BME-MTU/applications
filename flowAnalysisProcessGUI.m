function varargout = flowAnalysisProcessGUI(varargin)
% flowAnalysisProcessGUI M-file for flowAnalysisProcessGUI.fig
%      flowAnalysisProcessGUI, by itself, creates a new flowAnalysisProcessGUI or raises the existing
%      singleton*.
%
%      H = flowAnalysisProcessGUI returns the handle to a new flowAnalysisProcessGUI or the handle to
%      the existing singleton*.
%
%      flowAnalysisProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in flowAnalysisProcessGUI.M with the given input arguments.
%
%      flowAnalysisProcessGUI('Property','Value',...) creates a new flowAnalysisProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flowAnalysisProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flowAnalysisProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flowAnalysisProcessGUI

% Last Modified by GUIDE v2.5 30-Sep-2011 21:04:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flowAnalysisProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @flowAnalysisProcessGUI_OutputFcn, ...
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


% --- Executes just before flowAnalysisProcessGUI is made visible.
function flowAnalysisProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
userData.numericParams ={'timeWindow','corrLength','gridSize'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numericParams)

set(handles.edit_corrLengthMicrons,'String',...
    funParams.corrLength*userData.MD.pixelSize_/1000);
set(handles.edit_gridSizeMicrons,'String',...
    funParams.gridSize*userData.MD.pixelSize_/1000);

% Set popup-menu
flowProc = userData.crtProc.getFlowProcesses;
validProc =  cellfun(@(x) ~isempty(userData.MD.getProcessIndex(x,1)),flowProc);
flowProc=flowProc(validProc);

flowString = cellfun(@(x) eval([x '.getName']),flowProc,'Unif',false);
flowValue = find(strcmp(funParams.FlowProcess,flowProc));
if isempty(flowValue), flowValue=1; end
set(handles.popupmenu_FlowProcess,'String',flowString,'Value',flowValue,...
    'UserData',flowProc);

% Choose default command line output for speedMapsProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = flowAnalysisProcessGUI_OutputFcn(~, ~, handles) 
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

% Check user input
userData = get(handles.figure1, 'UserData');
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

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

for i=1:numel(userData.numericParams)  
   funParams.(userData.numericParams{i})=...
       str2double(get(handles.(['edit_' userData.numericParams{i}]),'String'));
end

flowProps = get(handles.popupmenu_FlowProcess,{'UserData','Value'});
funParams.FlowProcess = flowProps{1}{flowProps{2}};

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

function edit_corrLength_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
value = str2double(get(handles.edit_corrLength,'String'));
set(handles.edit_corrLengthMicrons,'String',value*userData.MD.pixelSize_/1000);

function edit_gridSize_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
value = str2double(get(handles.edit_gridSize,'String'));
set(handles.edit_gridSizeMicrons,'String',value*userData.MD.pixelSize_/1000);
