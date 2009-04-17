function varargout = plusTipTrackViz(varargin)
% PLUSTIPTRACKVIZ M-file for plusTipTrackViz.fig
%      PLUSTIPTRACKVIZ, by itself, creates a new PLUSTIPTRACKVIZ or raises the
%      existing
%      singleton*.
%
%      H = PLUSTIPTRACKVIZ returns the handle to a new PLUSTIPTRACKVIZ or the handle to
%      the existing singleton*.
%
%      PLUSTIPTRACKVIZ('CALLBACK',hObject,eventData,h,...) calls the local
%      function named CALLBACK in PLUSTIPTRACKVIZ.M with the given input arguments.
%
%      PLUSTIPTRACKVIZ('Property','Value',...) creates a new PLUSTIPTRACKVIZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipTrackViz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipTrackViz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipTrackViz

% Last Modified by GUIDE v2.5 16-Apr-2009 14:33:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipTrackViz_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipTrackViz_OutputFcn, ...
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


% --- Executes just before plusTipTrackViz is made visible.
function plusTipTrackViz_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% varargin   command line arguments to plusTipTrackViz (see VARARGIN)

% Choose default command line output for plusTipTrackViz
h.output = hObject;

h.homeDir=pwd;
h.projData=[];
h.tracksFinal=[];

h.roi=[];
h.timeRange=[1 inf];

h.doAvi=0;

h.indivTrack=[];
h.magCoef=[];
h.showTracks=1;
h.showDetect=1;

h.velLimit=inf;

h.img=[];
h.ask4select=0;
h.selectedTracks=[];
h.plotCurrentOnly=[];
h.movieInfo=[];

% Update h structure
guidata(hObject, h);

% UIWAIT makes plusTipTrackViz wait for user response (see UIRESUME)
% uiwait(h.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = plusTipTrackViz_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Get default command line output from h structure
varargout{1} = h.output;


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, h)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
plusTipGuiSwitch(hObject,eventdata,h,'resetButton');   

% --- Executes during object creation, after setting all properties.
function helpBoxAxes1_CreateFcn(hObject, eventdata, h)
% hObject    handle to helpBoxAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate helpBoxAxes1
img=imread('qIcon.jpg');
imagesc(img,'parent',hObject);
axis image
axis off
imHandle=get(hObject,'Children');
%info=get(hObject); info2=get(blah);
%set(blah,'HitTest','off')
set(imHandle,'ButtonDownFcn',@helpBoxAxes1_ButtonDownFcn)


% --- Executes on mouse press over axes background.
function helpBoxAxes1_ButtonDownFcn(hObject, eventdata, h)
% % hObject    handle to helpBoxAxes1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % h    structure with h and user data (see GUIDATA)
% % h.helpNote=1;
% % guidata(hObject, h);
% % helpNotes_Callback(hObject, eventdata, h)
% h.helpNotes
open plusTipTrackViz_README.txt


% --- Executes on button press in chooseProjData.
function chooseProjData_Callback(hObject, eventdata, h)
% hObject    handle to chooseProjData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
h=plusTipGuiSwitch(hObject,eventdata,h,'chooseProjData');
guidata(hObject, h);


% --- Executes on button press in selectSavedRoiPushbutton.
function selectSavedRoiPushbutton_Callback(hObject, eventdata, h)
% hObject    handle to selectSavedRoiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
h=plusTipGuiSwitch(hObject,eventdata,h,'selectSavedRoiPushbutton');
guidata(hObject, h);


function startFrame_Callback(hObject, eventdata, h)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrame as text
%        str2double(get(hObject,'String')) returns contents of startFrame as a double
h=plusTipGuiSwitch(hObject,eventdata,h,'startFrame');
guidata(hObject, h);

% --- Executes during object creation, after setting all properties.
function startFrame_CreateFcn(hObject, eventdata, h)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endFrame_Callback(hObject, eventdata, h)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrame as text
%        str2double(get(hObject,'String')) returns contents of endFrame as
%        a double
h=plusTipGuiSwitch(hObject,eventdata,h,'endFrame');  
guidata(hObject, h);

% --- Executes during object creation, after setting all properties.
function endFrame_CreateFcn(hObject, eventdata, h)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectTracksCheck.
function selectTracksCheck_Callback(hObject, eventdata, h)
% hObject    handle to selectTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectTracksCheck
h=plusTipGuiSwitch(hObject,eventdata,h,'selectTracksCheck');  
guidata(hObject, h);


% --- Executes on button press in showTracksCheck.
function showTracksCheck_Callback(hObject, eventdata, h)
% hObject    handle to showTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTracksCheck
h=plusTipGuiSwitch(hObject,eventdata,h,'showTracksCheck');  
guidata(hObject, h);


function speedLimitEdit_Callback(hObject, eventdata, h)
% hObject    handle to speedLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of
%        speedLimitEdit as a double
   h=plusTipGuiSwitch(hObject,eventdata,h,'speedLimitEdit');   
guidata(hObject, h);

% --- Executes during object creation, after setting all properties.
function speedLimitEdit_CreateFcn(hObject, eventdata, h)
% hObject    handle to speedLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function indivTrackNumbersEdit_Callback(hObject, eventdata, h)
% hObject    handle to indivTrackNumbersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of indivTrackNumbersEdit as text
%        str2double(get(hObject,'String')) returns contents of
%        indivTrackNumbersEdit as a double
h=plusTipGuiSwitch(hObject,eventdata,h,'indivTrackNumbersEdit');   
guidata(hObject, h);

% --- Executes during object creation, after setting all properties.
function indivTrackNumbersEdit_CreateFcn(hObject, eventdata, h)
% hObject    handle to indivTrackNumbersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotTracksPush.
function plotTracksPush_Callback(hObject, eventdata, h)
% hObject    handle to plotTracksPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
h=plusTipGuiSwitch(hObject,eventdata,h,'selectedTracksDisplay');   
h=plusTipGuiSwitch(hObject,eventdata,h,'plotTracksPush');
h=plusTipGuiSwitch(hObject,eventdata,h,'selectedTracksDisplay');   
guidata(hObject, h);


% --- Executes on button press in aviCheckTrackMov.
function aviCheckTrackMov_Callback(hObject, eventdata, h)
% hObject    handle to aviCheckTrackMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckTrackMov
h=plusTipGuiSwitch(hObject,eventdata,h,'aviCheckTrackMov');   
guidata(hObject, h);


function selectedTracksDisplay_Callback(hObject, eventdata, h)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectedTracksDisplay as text
%        str2double(get(hObject,'String')) returns contents of
%        selectedTracksDisplay as a double
h=plusTipGuiSwitch(hObject,eventdata,h,'selectedTracksDisplay');   
guidata(hObject, h);

% --- Executes during object creation, after setting all properties.
function selectedTracksDisplay_CreateFcn(hObject, eventdata, h)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called


% --- Executes on button press in speedMovieButton.
function speedMovieButton_Callback(hObject, eventdata, h)
% hObject    handle to speedMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
h=plusTipGuiSwitch(hObject,eventdata,h,'speedMovieButton');   
guidata(hObject, h);


% --- Executes on button press in trackMovieButton.
function trackMovieButton_Callback(hObject, eventdata, h)
% hObject    handle to trackMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
h=plusTipGuiSwitch(hObject,eventdata,h,'trackMovieButton');   
guidata(hObject, h);


% --- Executes on button press in aviCheckSpeedMov.
function aviCheckSpeedMov_Callback(hObject, eventdata, h)
% hObject    handle to aviCheckSpeedMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckSpeedMov
h=plusTipGuiSwitch(hObject,eventdata,h,'aviCheckSpeedMov');   
guidata(hObject, h);


% --- Executes on button press in detectionRadio1.
function detectionRadio1_Callback(hObject, eventdata, h)
% hObject    handle to detectionRadio1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio1
h.showDetect=1;
guidata(hObject, h);


% --- Executes on button press in detectionRadio2.
function detectionRadio2_Callback(hObject, eventdata, h)
% hObject    handle to detectionRadio2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio2
h.showDetect=2;
guidata(hObject, h);


% --- Executes on button press in detectionRadio3.
function detectionRadio3_Callback(hObject, eventdata, h)
% hObject    handle to detectionRadio3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio3
h.showDetect=3;
guidata(hObject, h);

% --- Executes on button press in detectionRadio4.
function detectionRadio4_Callback(hObject, eventdata, h)
% hObject    handle to detectionRadio4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio4
h.showDetect=0;
guidata(hObject, h);

% --- Executes when selected object is changed in radioButtonGroupDetection.
function radioButtonGroupDetection_SelectionChangeFcn(hObject, eventdata, h)
% hObject    handle to the selected object in radioButtonGroupDetection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% h    structure with h and user data (see GUIDATA)
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'detectionRadio1'
        detectionRadio1_Callback(hObject, eventdata, h)
    case 'detectionRadio2'
        detectionRadio2_Callback(hObject, eventdata, h)
    case 'detectionRadio3'
        detectionRadio3_Callback(hObject, eventdata, h)
    case 'detectionRadio3'
        detectionRadio4_Callback(hObject, eventdata, h)
end






