function varargout = speckleDetectionProcessGUI(varargin)
% speckleDetectionProcessGUI M-file for speckleDetectionProcessGUI.fig
%      speckleDetectionProcessGUI, by itself, creates a new speckleDetectionProcessGUI or raises the existing
%      singleton*.
%
%      H = speckleDetectionProcessGUI returns the handle to a new speckleDetectionProcessGUI or the handle to
%      the existing singleton*.
%
%      speckleDetectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in speckleDetectionProcessGUI.M with the given input arguments.
%
%      speckleDetectionProcessGUI('Property','Value',...) creates a new speckleDetectionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speckleDetectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speckleDetectionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speckleDetectionProcessGUI

% Last Modified by GUIDE v2.5 03-Oct-2011 16:14:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speckleDetectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @speckleDetectionProcessGUI_OutputFcn, ...
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


% --- Executes just before speckleDetectionProcessGUI is made visible.
function speckleDetectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set-up parameters
userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;

% Set-up parameters
set(handles.listbox_maskChannels,'Value',funParams.MaskChannelIndex);
set(handles.edit_alpha,'String',funParams.alpha);
set(handles.popupmenu_speckleOrder,'Value',funParams.paramSpeckles(1));
set(handles.edit_percEdit,'String',funParams.paramSpeckles(2));

% Set the same channels for the masks as for the detectable channels
props=get(handles.listbox_selectedChannels,{'String','UserData'});
set(handles.listbox_maskChannels,'String',props{1},'UserData',props{1});  

% Store the image directories and filterSigma values (for multi-channel)
userData.filterSigma = funParams.filterSigma;

% Update status of noise parameter fields
userData.noiseParams={'I0','sDN','GaussRatio'};
if ~isempty(userData.crtPackage.processes_{1})
    cellfun(@(x)set(handles.(['edit_' x]),'Enable','off'),userData.noiseParams);
    set(handles.pushbutton_loadNoiseParameters,'Enable','off');
else
    cellfun(@(x)set(handles.(['edit_' x]),'String',funParams.(x)),userData.noiseParams);
end
    
set(handles.figure1, 'UserData', userData);
% Update values of filter sigma and noise parameters
listbox_selectedChannels_Callback(hObject,eventdata,handles);

% Update GUI user data
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = speckleDetectionProcessGUI_OutputFcn(~, ~, handles) 
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

% -------- Check user input --------

if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

if isempty(get(handles.listbox_maskChannels, 'Value'))
    errordlg('Please select at least one mask channel from ''Mask Channels''.','Setting Error','modal')
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

% Retrieve GUI-defined parameters
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
maskChannelIndex = get(handles.listbox_maskChannels, 'Value');
funParams.MaskChannelIndex = maskChannelIndex;
funParams.alpha=str2double(get(handles.edit_alpha, 'String'));
funParams.paramSpeckles=[get(handles.popupmenu_speckleOrder,'Value')...
    str2double(get(handles.edit_percEdit, 'String'))];
% Save the filterSigma if different from psfSigma
% In order not to override filterSigma in batch movie set up
if ~isequal(userData.filterSigma,[userData.MD.channels_.psfSigma_])
    funParams.filterSigma = userData.filterSigma;
end


if isempty(userData.crtPackage.processes_{1})
    for i=1:numel(userData.noiseParams)
        funParams.(userData.noiseParams{i}) = ...
            str2num(get(handles.(['edit_' userData.noiseParams{i}]),'String')); %#ok<ST2NM>
    end
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on selection change in popupmenu_speckleOrder.
function popupmenu_speckleOrder_Callback(hObject, eventdata, handles)

if get(hObject,'Value')==1
    set(handles.edit_percEdit,'Enable','off');
else
    set(handles.edit_percEdit,'Enable','on');
end

% --- Executes on selection change in listbox_selectedChannels.
function listbox_selectedChannels_Callback(hObject, eventdata, handles)

% Read channel index
userData = get(handles.figure1, 'UserData');
props = get(handles.listbox_selectedChannels, {'UserData','Value'});
chanIndx = props{1}(props{2});

set(handles.edit_filterSigma,'String',userData.filterSigma(chanIndx));
if ~isempty(userData.crtPackage.processes_{1})
    if userData.crtPackage.processes_{1}.checkChannelOutput() 
        [I0,sDN,GaussRatio] =userData.crtPackage.processes_{1}.loadChannelOutput(chanIndx);
        set(handles.edit_I0,'String',I0);
        set(handles.edit_sDN,'String',sDN);
        set(handles.edit_GaussRatio,'String',GaussRatio);
    else
        cellfun(@(x)set(handles.(['edit_' x]),'String','NA'),userData.noiseParams);
    end
end

function edit_filterSigma_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
chanIndx = props{1}(props{2});

value = str2double(get(hObject, 'String'));
if ~(value>0), 
    % Reset old value
    set(hObject,'String',userData.filterSigma(chanIndx))
else
    % Update the sigma value in the stored array-
    userData.filterSigma(chanIndx) = value;
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_loadNoiseParameters.
function pushbutton_loadNoiseParameters_Callback(hObject, eventdata, handles)

[file path]=uigetfile({'*.mat;*.MAT','Mat files (*.mat)'},...
    'Select the file containing the noise model parameters');
if ~isequal(file,0) && ~isequal(path,0)
    noiseParams={'I0','sDN','GaussRatio'};
    vars=whos(noiseParams{:},'-file',[path file]);
    if numel(vars)~=numel(noiseParams),
        errordlg('Please select a file containing valid noise model parameters');
        return 
    end
    
    s=load([path file],noiseParams{:});
    for i=1:numel(noiseParams)
        set(handles.(['edit_' noiseParams{i}]),'String',s.(noiseParams{i}));
        
    end
end
