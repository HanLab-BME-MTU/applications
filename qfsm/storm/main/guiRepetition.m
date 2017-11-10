function [varargout] = guiRepetition(varargin)
% GUIREPETITION MATLAB code for guiRepetition.fig
%      GUIREPETITION, by itself, creates a new GUIREPETITION or raises the existing
%      singleton*.
%
%      H = GUIREPETITION returns the handle to a new GUIREPETITION or the handle to
%      the existing singleton*.
%
%      GUIREPETITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIREPETITION.M with the given input arguments.
%
%      GUIREPETITION('Property','Value',...) creates a new GUIREPETITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiRepetition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiRepetition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiRepetition

% Last Modified by GUIDE v2.5 07-Feb-2012 13:45:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiRepetition_OpeningFcn, ...
                   'gui_OutputFcn',  @guiRepetition_OutputFcn, ...
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


% --- Executes just before guiRepetition is made visible.
function guiRepetition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiRepetition (see VARARGIN)

% Choose default command line output for guiRepetition
disp('Main: guiRepetition can only be called by guiMain!');

handles.output = hObject;

handles.roiSelector = varargin{1};
handles.currentConfig = varargin{2};
handles.state = 0;
handles.roiSelector.offset = handles.currentConfig.roiSize(1:2);
set(handles.edit1,'String',num2str(handles.roiSelector.repetitions(1)));
set(handles.edit5,'String',num2str(handles.roiSelector.repetitions(2)));
set(handles.edit6,'String',num2str(handles.roiSelector.offset(1)));
set(handles.edit7,'String',num2str(handles.roiSelector.offset(2)));

% Update handles structure
guidata(hObject,handles);

drawRepetitions(handles);

% UIWAIT makes guiRepetition wait for user response (see UIRESUME)
uiwait(handles.figure1);


function drawRepetitions(handles)
imagesc(handles.roiSelector.preview);
colormap(hot);
set(handles.axes1,'xticklabel',[]);
set(handles.axes1,'yticklabel',[]);
set(handles.axes1,'XTick',[]);
set(handles.axes1,'YTick',[]);
for x=1:handles.roiSelector.repetitions(1)
    for y=1:handles.roiSelector.repetitions(2)
        pos = handles.currentConfig.roiPosition(1:2) + ([x y]-1).*handles.roiSelector.offset;
        siz = handles.currentConfig.roiSize(1:2);
        position = [pos siz]*handles.roiSelector.previewScaleFactor;
        rectangle('Position',position,'EdgeColor',[1 0 0],'LineWidth',2);%,'FaceColor','r')
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = guiRepetition_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.state;
delete(handles.figure1);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
disp('e1')
handles.roiSelector.repetitions(1,1) = str2double(get(hObject,'String'));
drawRepetitions(handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
disp('e5')
handles.roiSelector.repetitions(1,2) = str2double(get(hObject,'String'));
drawRepetitions(handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
disp('e6')
handles.roiSelector.offset(1,1) = str2double(get(hObject,'String'));
drawRepetitions(handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
disp('e7')
handles.roiSelector.offset(1,2) = str2double(get(hObject,'String'));
drawRepetitions(handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb1')
handles.state = 0;
guidata(hObject,handles);
close(gcf);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb2')
handles.state = 1;
guidata(hObject,handles);
close(gcf);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume();
% delete(hObject);
