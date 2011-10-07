function varargout = windowingProcessGUI(varargin)
% windowingProcessGUI M-file for windowingProcessGUI.fig
%      windowingProcessGUI, by itself, creates a new windowingProcessGUI or raises the existing
%      singleton*.
%
%      H = windowingProcessGUI returns the handle to a new windowingProcessGUI or the handle to
%      the existing singleton*.
%
%      windowingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in windowingProcessGUI.M with the given input arguments.
%
%      windowingProcessGUI('Property','Value',...) creates a new windowingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before windowingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to windowingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help windowingProcessGUI

% Last Modified by GUIDE v2.5 10-Aug-2011 09:25:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @windowingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @windowingProcessGUI_OutputFcn, ...
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


% --- Executes just before windowingProcessGUI is made visible.
function windowingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
userData.numParams ={'ParaSize','PerpSize','MinSize'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams)

% Read available segmentaation processes
segProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
segProcID=find(segProc);
segProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(segProc),'Unif',false);
segProcString = vertcat('Choose later',segProcNames(:));
segProcData=horzcat({[]},num2cell(segProcID));

% Read the default segmentation process index
% If empty, try to propagate segmentation process from protrusion process
initSegProcIndex = funParams.SegProcessIndex;
if isempty(initSegProcIndex) && ~isempty(userData.crtPackage.processes_{1});
    if ~isempty(userData.crtPackage.processes_{1}.funParams_.SegProcessIndex)
        initSegProcIndex = userData.crtPackage.processes_{1}.funParams_.SegProcessIndex;  
    end
end
segProcValue = find(cellfun(@(x) isequal(x,initSegProcIndex),segProcData));
set(handles.popupmenu_SegProcessIndex,'String',segProcString,...
    'UserData',segProcData,'Value',segProcValue);

% Create pop-up menu for windowing methods
methodString ={'Constant number';'Constant width';'Protrusion based';'PDE based'};
methodData ={'ConstantNumber';'ConstantWidth';'ProtrusionBased';'PDEBased'};
methodValue = find(strcmpi(funParams.MethodName,methodData));
set(handles.popupmenu_MethodName,'String',methodString,...
    'UserData',methodData,'Value',methodValue);

% Set-up pde parameters
pdeString ={'Select a model';'Viscous';'Viscous-convective';'Visco-elastic'};
pdeData ={'';'Viscous';'ViscousConvective';'ViscoElastic'};
if ~isempty(funParams.PDEPar)
    PDEValue = find(strcmpi(funParams.PDEPar,pdeData));
else
    PDEValue =1;
end
set(handles.popupmenu_PDEPar,'Value',PDEValue);
set(handles.edit_MeshQuality,'String',funParams.MeshQuality);
value=strcmpi(funParams.NonLinearSolver,'on');
set(handles.checkbox_NonLinearSolver,'Value',value);
set(handles.popupmenu_PDEPar,'String',pdeString,'UserData',pdeData);

popupmenu_MethodName_Callback(hObject, eventdata, handles)

% Choose default command line output for windowingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = windowingProcessGUI_OutputFcn(~, ~, handles) 
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

for i=1:numel(userData.numParams)  
    value = get(handles.(['edit_' userData.numParams{i}]),'String');
    if isempty(value)
        errordlg(['Please enter a valid value for the '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],'Setting Error','modal');
        return;
    end
    funParams.(userData.numParams{i})=str2double(value); 
end

%Retrieve segmentation process
props=get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
funParams.SegProcessIndex=props{1}{props{2}};

% Retrieve windowing method
props=get(handles.popupmenu_MethodName,{'UserData','Value'});
funParams.MethodName=props{1}{props{2}};

if strcmpi(funParams.MethodName,'pdebased');
    meshQuality = get(handles.edit_MeshQuality,'String');
    if isempty(meshQuality)
        errordlg('Please enter a valid value for the Quality of the triangular mesh.','Setting Error','modal');
        return;
    end
    funParams.MeshQuality=str2double(meshQuality);       
    props = get(handles.popupmenu_PDEPar,{'UserData','Value'});
    if props{2}==1
        errordlg('Select a valid window propagation method.','Setting Error','modal');
        return;
    end
    funParams.PDEPar=props{1}{props{2}};
    nonLinearSolver=get(handles.checkbox_NonLinearSolver,'Value');
    if nonLinearSolver
        funParams.NonLinearSolver='on';
    else
        funParams.NonLinearSolver='off';
    end
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


% --- Executes on selection change in popupmenu_MethodName.
function popupmenu_MethodName_Callback(hObject, eventdata, handles)

props=get(handles.popupmenu_MethodName,{'UserData','Value'});
methodName=props{1}{props{2}};

if strcmpi(methodName,'PDEBased')
    enableState='on';
else
    enableState='off';
end
set(get(handles.uipanel_pdeParameters,'Children'),'Enable',enableState);
