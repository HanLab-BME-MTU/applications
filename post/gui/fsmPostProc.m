function varargout = fsmPostProc(varargin)
% FSMPOSTPROC M-file for fsmPostProc.fig
%      FSMPOSTPROC, by itself, creates a new FSMPOSTPROC or raises the existing
%      singleton*.
%
%      H = FSMPOSTPROC returns the handle to a new FSMPOSTPROC or the handle to
%      the existing singleton*.
%
%      FSMPOSTPROC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMPOSTPROC.M with the given input arguments.
%
%      FSMPOSTPROC('Property','Value',...) creates a new FSMPOSTPROC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmPostProc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmPostProc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmPostProc

% Last Modified by GUIDE v2.5 30-Aug-2004 19:28:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmPostProc_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmPostProc_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fsmPostProc is made visible.
function fsmPostProc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmPostProc (see VARARGIN)

% Choose default command line output for fsmPostProc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fsmPostProc wait for user response (see UIRESUME)
% uiwait(handles.fsmPostProcessing);


% --- Outputs from this function are returned to the command line.
function varargout = fsmPostProc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushCurrentProject.
function pushCurrentProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushCurrentProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushProjectInfo.
function pushProjectInfo_Callback(hObject, eventdata, handles)
% hObject    handle to pushProjectInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


