function varargout = speckleTrackingProcessGUI(varargin)
% speckleTrackingProcessGUI M-file for speckleTrackingProcessGUI.fig
%      speckleTrackingProcessGUI, by itself, creates a new speckleTrackingProcessGUI or raises the existing
%      singleton*.
%
%      H = speckleTrackingProcessGUI returns the handle to a new speckleTrackingProcessGUI or the handle to
%      the existing singleton*.
%
%      speckleTrackingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in speckleTrackingProcessGUI.M with the given input arguments.
%
%      speckleTrackingProcessGUI('Property','Value',...) creates a new speckleTrackingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speckleTrackingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speckleTrackingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speckleTrackingProcessGUI

% Last Modified by GUIDE v2.5 30-Sep-2011 21:01:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speckleTrackingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @speckleTrackingProcessGUI_OutputFcn, ...
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


% --- Executes just before speckleTrackingProcessGUI is made visible.
function speckleTrackingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% ---------------------- Parameter Setup -------------------------

userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
set(handles.edit_corrLength,'String',funParams.corrLength);
set(handles.edit_corrLengthMicrons,'String',...
    funParams.corrLength*userData.MD.pixelSize_/1000);
set(handles.edit_threshold,'String',funParams.threshold);
set(handles.checkbox_enhanced,'Value',funParams.enhanced);

% Enable interpolation methods only if flow tracking is set up
if ~isempty(userData.crtPackage.processes_{5});
    set(handles.popupmenu_interpolationMethod,'Enable','on');
else
    set(handles.popupmenu_interpolationMethod,'Enable','off');
    % If neither flow tracking nor speckle tracking is set up
    if isempty(userData.crtPackage.processes_{6});
        set(handles.checkbox_enhanced,'Value',1);
    end
end

% Load interpolation methods
intMethods=SpeckleTrackingProcess.getInterpolationMethods;
methodDescriptions = {intMethods(:).description};
methodNames = {intMethods(:).name};
try
    methodValue = find(strcmp(funParams.interpolationMethod,methodNames));
catch ME
    methodValue=1;
end
if isempty(methodValue), methodValue=1; end
set(handles.popupmenu_interpolationMethod,'String',methodDescriptions,...
    'UserData',methodNames,'Value',methodValue);

% Choose default command line output for speckleTrackingProcessGUI
handles.output = hObject;
% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = speckleTrackingProcessGUI_OutputFcn(~, ~, handles) 
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

% -------- Check user input --------

if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

% -------- Process Sanity check --------
% ( only check underlying data )
userData = get(handles.figure1, 'UserData');
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% -------- Set parameter --------

channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
funParams.corrLength =str2double(get(handles.edit_corrLength,'String'));
funParams.threshold =str2double(get(handles.edit_threshold,'String'));
funParams.enhanced =get(handles.checkbox_enhanced,'Value');

props = get(handles.popupmenu_interpolationMethod,{'UserData','Value'});
funParams.interpolationMethod = props{1}{props{2}};

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



function edit_corrLength_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
value = str2double(get(handles.edit_corrLength,'String'));
set(handles.edit_corrLengthMicrons,'String',value*userData.MD.pixelSize_/1000);
