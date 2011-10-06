function varargout = windowSamplingProcessGUI(varargin)
% windowSamplingProcessGUI M-file for windowSamplingProcessGUI.fig
%      windowSamplingProcessGUI, by itself, creates a new windowSamplingProcessGUI or raises the existing
%      singleton*.
%
%      H = windowSamplingProcessGUI returns the handle to a new windowSamplingProcessGUI or the handle to
%      the existing singleton*.
%
%      windowSamplingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in windowSamplingProcessGUI.M with the given input arguments.
%
%      windowSamplingProcessGUI('Property','Value',...) creates a new windowSamplingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before windowSamplingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to windowSamplingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help windowSamplingProcessGUI

% Last Modified by GUIDE v2.5 28-Sep-2011 11:43:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @windowSamplingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @windowSamplingProcessGUI_OutputFcn, ...
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


% --- Executes just before windowSamplingProcessGUI is made visible.
function windowSamplingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);
% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;


% Read available segmentaation processes
imProc =  cellfun(@(x) isa(x,'ImageProcessingProcess'),userData.MD.processes_);
imProcID=find(imProc);
imProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(imProc),'Unif',false);
imProcString = vertcat('Raw images',imProcNames(:));
imProcData=horzcat({[]},num2cell(imProcID));


% Read the default segmentation process index
% If empty, try to propagate segmentation process from protrusion process
initProcessIndex = funParams.ProcessIndex;
initProcessValue = find(cellfun(@(x) isequal(x,initProcessIndex),imProcData));
set(handles.popupmenu_ProcessIndex,'String',imProcString,...
    'UserData',imProcData,'Value',initProcessValue);

% Update channels listboxes depending on the selected process
popupmenu_ProcessIndex_Callback(hObject,eventdata,handles);

% Choose default command line output for windowSamplingProcessGUI
handles.output = hObject;

% Update user data and GUI data
% set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = windowSamplingProcessGUI_OutputFcn(~, ~, handles) 
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

% Retrieve selected process ID
props= get(handles.popupmenu_ProcessIndex,{'UserData','Value'});
procID = props{1}{props{2}};
funParams.ProcessIndex=procID;

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


% --- Executes on selection change in popupmenu_ProcessIndex.
function popupmenu_ProcessIndex_Callback(hObject, eventdata, handles)

% Retrieve selected process ID
props= get(handles.popupmenu_ProcessIndex,{'UserData','Value'});
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
