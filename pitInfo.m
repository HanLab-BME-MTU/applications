function varargout = pitInfo(varargin)
% PITINFO M-file for pitInfo.fig
%      PITINFO, by itself, creates a new PITINFO or raises the existing
%      singleton*.
%
%      H = PITINFO returns the handle to a new PITINFO or the handle to
%      the existing singleton*.
%
%      PITINFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PITINFO.M with the given input arguments.
%
%      PITINFO('Property','Value',...) creates a new PITINFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pitInfo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pitInfo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pitInfo

% Last Modified by GUIDE v2.5 23-Feb-2010 20:34:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pitInfo_OpeningFcn, ...
                   'gui_OutputFcn',  @pitInfo_OutputFcn, ...
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

%OPENING FUNCTION
% --- Executes just before pitInfo is made visible.
function pitInfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pitInfo (see VARARGIN)
% Choose default command line output for pitInfo
% Update handles structure
handles = guidata(hObject);
handles.output = hObject;
% UIWAIT makes pitInfo wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% Update handles structure
guidata(hObject,handles)

%PLOT INTENSITY DATA
function makePlots(handles)
%get handles from endocytosisGUI
mainHandles = guidata(handles.endocytosisGUIHandle);
%update lifetime handles
set(handles.text1,'String',['pit lifetime: ' num2str(mainHandles.pitInfo.lifetime) ' seconds'])
% turn zoom on for axes of channel one
zoom(handles.axes1,'on')
%plot channel both chanels if there is data for channel2
if isfield(mainHandles.pitInfo,'intensityChannel2')
     [AX,H1,H2] = plotyy(handles.axes1,...
        1-mainHandles.pitInfo.initiationFrame:size(mainHandles.pitInfo.intensityReferenceChannel1,2) - mainHandles.pitInfo.initiationFrame,...
    mainHandles.pitInfo.intensityChannel1(mainHandles.pitInfo.pitID,:,1)-mainHandles.pitInfo.intensityReferenceChannel1(mainHandles.pitInfo.pitID,:,1),...
        1-mainHandles.pitInfo.initiationFrame:size(mainHandles.pitInfo.intensityReferenceChannel1,2) - mainHandles.pitInfo.initiationFrame,...
    mainHandles.pitInfo.intensityChannel2(mainHandles.pitInfo.pitID,:,1)-mainHandles.pitInfo.intensityReferenceChannel2(mainHandles.pitInfo.pitID,:,1));
set(H1,'Color','b','LineWidth',2)
set(H2,'Color','g','LineWidth',2)
ylabel(AX(1),'Channel1 intensity (background substracted)')
ylabel(AX(2),'Channel2 intensity (background substracted)')
else
    %plot channel 1 only
   plot(handles.axes1,...
       1-mainHandles.pitInfo.initiationFrame:size(mainHandles.pitInfo.intensityReferenceChannel1,2) - mainHandles.pitInfo.initiationFrame,...
    mainHandles.pitInfo.intensityChannel1(mainHandles.pitInfo.pitID,:,1)-mainHandles.pitInfo.intensityReferenceChannel1(mainHandles.pitInfo.pitID,:,1)...
    ,'g','LineWidth',2) 
end
%
handles.axes1_2 = AX(2);
handles.axes1_1 = AX(1);
% Update handles structure
guidata(handles.figure1,handles)


%OUTPUT
% --- Outputs from this function are returned to the command line.
function varargout = pitInfo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%TEXT FOR LIFETIME
% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%AXES FOR INTENSITY PLOTS
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
