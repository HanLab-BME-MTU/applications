function varargout = protrusionProcessGUI(varargin)
% protrusionProcessGUI M-file for protrusionProcessGUI.fig
%      protrusionProcessGUI, by itself, creates a new protrusionProcessGUI or raises the existing
%      singleton*.
%
%      H = protrusionProcessGUI returns the handle to a new protrusionProcessGUI or the handle to
%      the existing singleton*.
%
%      protrusionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in protrusionProcessGUI.M with the given input arguments.
%
%      protrusionProcessGUI('Property','Value',...) creates a new protrusionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before protrusionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to protrusionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help protrusionProcessGUI

% Last Modified by GUIDE v2.5 01-Mar-2012 20:30:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @protrusionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @protrusionProcessGUI_OutputFcn, ...
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


% --- Executes just before protrusionProcessGUI is made visible.
function protrusionProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

segProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
segProcID=find(segProc);
segProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(segProc),'Unif',false);
segProcString = vertcat('Choose later',segProcNames(:));
segProcData=horzcat({[]},num2cell(segProcID));
segProcValue = find(cellfun(@(x) isequal(x,funParams.SegProcessIndex),segProcData));
if isempty(segProcValue), segProcValue = 1; end
set(handles.popupmenu_SegProcessIndex,'String',segProcString,...
    'UserData',segProcData,'Value',segProcValue);

userData.numParams ={'DownSample','SplineTolerance'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams)

% Update channels listboxes depending on the selected process
popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

% Choose default command line output for protrusionProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = protrusionProcessGUI_OutputFcn(~, ~, handles) 
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

% Retrieve mask process index and class (for propagation)
props=get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
funParams.SegProcessIndex = props{1}{props{2}};
if ~isempty(funParams.SegProcessIndex)
    segProcessClass=class(userData.MD.processes_{funParams.SegProcessIndex});
else
    segProcessClass = '';
end

for i=1:numel(userData.numParams)
    value = get(handles.(['edit_' userData.numParams{i}]),'String');
    if isempty(value)
        errordlg(['Please enter a valid value for the '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],'Setting Error','modal');
        return;
    end
    funParams.(userData.numParams{i})=str2double(value);
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
setMaskProcess = @(x) parseProcessParams(x, struct('SegProcessIndex',...
    x.owner_.getProcessIndex(segProcessClass,1,false)));
processGUI_ApplyFcn(hObject, eventdata, handles,funParams,'settingFcn',{setMaskProcess});


% --- Executes on selection change in popupmenu_SegProcessIndex.
function popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

% Retrieve selected process ID
props= get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
procID = props{1}{props{2}};

% Read process and check available channels
userData = get(handles.figure1, 'UserData');
if isempty(procID)
    allChannelIndex=1:numel(userData.MD.channels_);
else
    allChannelIndex = find(userData.MD.processes_{procID}.checkChannelOutput);
end

% Set up available channels listbox
if ~isempty(allChannelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(allChannelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,allChannelIndex);
    end
else
    channelString = {};
end
set(handles.listbox_availableChannels,'String',channelString,'UserData',allChannelIndex);

% Set up selected channels listbox
channelIndex = get(handles.listbox_selectedChannels, 'UserData');
channelIndex = intersect(channelIndex,allChannelIndex);
if ~isempty(channelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(channelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,channelIndex);
    end
else
    channelString = {};
end
set(handles.listbox_selectedChannels,'String',channelString,'UserData',channelIndex);