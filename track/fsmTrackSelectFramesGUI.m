function varargout = fsmTrackSelectFramesGUI(varargin)
% FSMTRACKSELECTFRAMESGUI M-file for fsmTrackSelectFramesGUI.fig
%      FSMTRACKSELECTFRAMESGUI, by itself, creates a new FSMTRACKSELECTFRAMESGUI or raises the existing
%      singleton*.
%
%      H = FSMTRACKSELECTFRAMESGUI returns the handle to a new FSMTRACKSELECTFRAMESGUI or the handle to
%      the existing singleton*.
%
%      FSMTRACKSELECTFRAMESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMTRACKSELECTFRAMESGUI.M with the given input arguments.
%
%      FSMTRACKSELECTFRAMESGUI('Property','Value',...) creates a new FSMTRACKSELECTFRAMESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmTrackSelectFramesGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmTrackSelectFramesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmTrackSelectFramesGUI

% Last Modified by GUIDE v2.5 12-Jan-2004 09:31:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmTrackSelectFramesGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmTrackSelectFramesGUI_OutputFcn, ...
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


% --- Executes just before fsmTrackSelectFramesGUI is made visible.
function fsmTrackSelectFramesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmTrackSelectFramesGUI (see VARARGIN)

% Choose default command line output for fsmTrackSelectFramesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Read input parameters
first=varargin{1};
last=varargin{2};
minN=varargin{3};
titleGUI=char(varargin{4});

% Store minN
set(handles.pushOkay,'UserData',minN);

% Check input
if last-minN<first
    error('Incompatible values for minimum, maximum, and step.');
end
if last==first
    set(handles.SelectFramesGUI,'UserData',{first,last});
    return
end

% Set up GUI with the passed data
set(handles.editFirstFrame,'String',num2str(first));
set(handles.editLastFrame,'String',num2str(last));
set(findobj('Tag','SelectFramesGUI'),'Name',titleGUI);
sSteps=[1/(last-first) 1/(last-first)];
set(handles.sliderFirstFrame,'SliderStep',sSteps,'Max',last,'Min',first,'Value',first);
set(handles.sliderLastFrame,'SliderStep',sSteps,'Max',last,'Min',first,'Value',last);

% UIWAIT makes fsmTrackSelectFramesGUI wait for user response (see UIRESUME)
uiwait(handles.SelectFramesGUI);


% --- Outputs from this function are returned to the command line.
function varargout = fsmTrackSelectFramesGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
selection=get(handles.SelectFramesGUI,'UserData');
varargout{1}=selection{1};
varargout{2}=selection{2};
delete(hObject);

% --- Executes during object creation, after setting all properties.
function sliderFirstFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function sliderFirstFrame_Callback(hObject, eventdata, handles)
% hObject    handle to sliderFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minN=get(handles.pushOkay,'UserData');
valueSlider1=fix(get(handles.sliderFirstFrame,'Value'));
valueSlider2=fix(get(handles.sliderLastFrame,'Value'));
if valueSlider1>(valueSlider2-minN)
    valueSlider1=valueSlider2-minN;    
end
set(handles.sliderFirstFrame,'Value',valueSlider1);
set(handles.editFirstFrame,'String',num2str(valueSlider1));

% --- Executes during object creation, after setting all properties.
function editFirstFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editFirstFrame_Callback(hObject, eventdata, handles)
% hObject    handle to editFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFirstFrame as text
%        str2double(get(hObject,'String')) returns contents of editFirstFrame as a double
minN=get(handles.pushOkay,'UserData');
value=str2num(get(handles.editFirstFrame,'String'));
if ~isnumeric(value)
    return
end
valueSlider1=get(handles.sliderFirstFrame,'Value');
valueSlider2=get(handles.sliderLastFrame,'Value');
if value<get(handles.sliderFirstFrame,'Min') | value>(valueSlider2-minN) | value>get(handles.sliderFirstFrame,'Max')
    value=valueSlider1;    
end
set(handles.sliderFirstFrame,'Value',value);
set(handles.editFirstFrame,'String',num2str(value));


% --- Executes during object creation, after setting all properties.
function sliderLastFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLastFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function sliderLastFrame_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLastFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minN=get(handles.pushOkay,'UserData');
valueSlider1=fix(get(handles.sliderFirstFrame,'Value'));
valueSlider2=fix(get(handles.sliderLastFrame,'Value'));
if valueSlider2<(valueSlider1+minN)
    valueSlider2=valueSlider1+minN;    
end
set(handles.sliderLastFrame,'Value',valueSlider2);
set(handles.editLastFrame,'String',num2str(valueSlider2));


% --- Executes during object creation, after setting all properties.
function editLastFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLastFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editLastFrame_Callback(hObject, eventdata, handles)
% hObject    handle to editLastFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLastFrame as text
%        str2double(get(hObject,'String')) returns contents of editLastFrame as a double
minN=get(handles.pushOkay,'UserData');
value=str2num(get(handles.editLastFrame,'String'));
if ~isnumeric(value)
    return
end
valueSlider1=get(handles.sliderFirstFrame,'Value');
valueSlider2=get(handles.sliderLastFrame,'Value');
if value<get(handles.sliderLastFrame,'Min') | value<(valueSlider1+minN) | value>get(handles.sliderLastFrame,'Max')
    value=valueSlider2;    
end
set(handles.sliderLastFrame,'Value',value);
set(handles.editLastFrame,'String',num2str(value));


% --- Executes on button press in pushOkay.
function pushOkay_Callback(hObject, eventdata, handles)
% hObject    handle to pushOkay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% assignin('caller','uFirst',first);
% assignin('caller','uLast',last);
uFirst=get(handles.sliderFirstFrame,'Value');
uLast=get(handles.sliderLastFrame,'Value');
set(handles.SelectFramesGUI,'UserData',{uFirst,uLast});
uiresume(handles.SelectFramesGUI);


% --- Executes when user attempts to close SelectFramesGUI.
function SelectFramesGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to SelectFramesGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set output to {-1, -1} (-1 means failure - this is what will be passed back if the gui is closed)
set(handles.SelectFramesGUI,'UserData',{-1,-1});
uiresume(handles.SelectFramesGUI);

