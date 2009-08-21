function varargout = plusTipTrackViz(varargin)
% PLUSTIPTRACKVIZ M-file for plusTipTrackViz.fig
%      PLUSTIPTRACKVIZ, by itself, creates a new PLUSTIPTRACKVIZ or raises the
%      existing
%      singleton*.
%
%      H = PLUSTIPTRACKVIZ returns the handle to a new PLUSTIPTRACKVIZ or the handle to
%      the existing singleton*.
%
%      PLUSTIPTRACKVIZ('CALLBACK',hObject,eventData,handles,...) calls the local
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

% Last Modified by GUIDE v2.5 13-Aug-2009 06:55:40


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
function plusTipTrackViz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipTrackViz (see VARARGIN)

% Choose default command line output for plusTipTrackViz
handles.output = hObject;

handles.getStr = 0;
handles.loadProjList = 0;
handles.loadProjList = 0;

handles.homeDir=pwd;
handles.projData=[];
handles.tracksFinal=[];

handles.roi=[];
handles.timeRange=[1 inf];

handles.doAvi=0;

handles.indivTrack=[];
handles.magCoef=[];
handles.showTracks=1;
handles.showDetect=1;

handles.velLimit=inf;

handles.img=[];
handles.ask4select=0;
handles.selectedTracks=[];
handles.plotCurrentOnly=[];
handles.movieInfo=[];


handles.xaxisScatter='growthSpeed';
handles.yaxisScatter='growthLifetime';
handles.xScatterPercent=50;
handles.yScatterPercent=50;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plusTipTrackViz wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = plusTipTrackViz_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plusTipGuiSwitch(hObject,eventdata,handles,'resetButton');


% --- Executes on button press in getProjPush.
function getProjPush_Callback(hObject, eventdata, handles)
% hObject    handle to getProjPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plusTipGuiSwitch(hObject,eventdata,handles,'getProjPush');

if ~isempty(handles.projList)
    % here we do NOT filter out any sub-directories
    a=struct2cell(handles.projList);
    a=a(2,:)';
    a=sort(a);
    
    % check for existence of projData in meta folder - if it's not there,
    % the project has not been tracked/processed, so we cannot visualize it
    b=zeros(length(a),1);
    for i=1:length(a)
        b(i)=exist([formatPath(a{i}) filesep 'meta' filesep 'projData.mat'],'file')==2;
    end
    a=a(logical(b));

    % allow only one project to be selected
    [selection,selectionList]=listSelectGUI(a,1,'move',1);

    % if a project was selected, save projData info and get data
    if ~isempty(selection)
        handles.dataDir=formatPath(selectionList{1,1});
        p=load([handles.dataDir filesep 'meta' filesep 'projData.mat']);
        handles.projData=p.projData;
    else
        msgbox('No projects selected or tracking has not been completed.')
        handles.dataDir=[];
        handles.projData=[];
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
    handles.projData=[];
end
guidata(hObject, handles);


% --- Executes on button press in selectSavedRoiPushbutton.
function selectSavedRoiPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectSavedRoiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectSavedRoiPushbutton');
guidata(hObject, handles);


function startFrame_Callback(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrame as text
%        str2double(get(hObject,'String')) returns contents of startFrame as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFrame');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function startFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endFrame_Callback(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrame as text
%        str2double(get(hObject,'String')) returns contents of endFrame as
%        a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrame');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function endFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectTracksCheck.
function selectTracksCheck_Callback(hObject, eventdata, handles)
% hObject    handle to selectTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectTracksCheck
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectTracksCheck');
guidata(hObject, handles);


% --- Executes on button press in showTracksCheck.
function showTracksCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTracksCheck
handles=plusTipGuiSwitch(hObject,eventdata,handles,'showTracksCheck');
guidata(hObject, handles);


function speedLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to speedLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of
%        speedLimitEdit as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'speedLimitEdit');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function speedLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function indivTrackNumbersEdit_Callback(hObject, eventdata, handles)
% hObject    handle to indivTrackNumbersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of indivTrackNumbersEdit as text
%        str2double(get(hObject,'String')) returns contents of
%        indivTrackNumbersEdit as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'indivTrackNumbersEdit');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function indivTrackNumbersEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to indivTrackNumbersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotTracksPush.
function plotTracksPush_Callback(hObject, eventdata, handles)
% hObject    handle to plotTracksPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectedTracksDisplay');
handles=plusTipGuiSwitch(hObject,eventdata,handles,'plotTracksPush');
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectedTracksDisplay');
guidata(hObject, handles);


% --- Executes on button press in aviCheckTrackMov.
function aviCheckTrackMov_Callback(hObject, eventdata, handles)
% hObject    handle to aviCheckTrackMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckTrackMov
handles=plusTipGuiSwitch(hObject,eventdata,handles,'aviCheckTrackMov');
guidata(hObject, handles);


function selectedTracksDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectedTracksDisplay as text
%        str2double(get(hObject,'String')) returns contents of
%        selectedTracksDisplay as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectedTracksDisplay');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function selectedTracksDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in speedMovieButton.
function speedMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to speedMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'speedMovieButton');
guidata(hObject, handles);


% --- Executes on button press in trackMovieButton.
function trackMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to trackMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'trackMovieButton');
guidata(hObject, handles);


% --- Executes on button press in aviCheckSpeedMov.
function aviCheckSpeedMov_Callback(hObject, eventdata, handles)
% hObject    handle to aviCheckSpeedMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckSpeedMov
handles=plusTipGuiSwitch(hObject,eventdata,handles,'aviCheckSpeedMov');
guidata(hObject, handles);


% --- Executes on button press in detectionRadio1.
function detectionRadio1_Callback(hObject, eventdata, handles)
% hObject    handle to detectionRadio1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio1
handles.showDetect=1;
guidata(hObject, handles);


% --- Executes on button press in detectionRadio2.
function detectionRadio2_Callback(hObject, eventdata, handles)
% hObject    handle to detectionRadio2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio2
handles.showDetect=2;
guidata(hObject, handles);


% --- Executes on button press in detectionRadio3.
function detectionRadio3_Callback(hObject, eventdata, handles)
% hObject    handle to detectionRadio3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio3
handles.showDetect=3;
guidata(hObject, handles);

% --- Executes on button press in detectionRadio4.
function detectionRadio4_Callback(hObject, eventdata, handles)
% hObject    handle to detectionRadio4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectionRadio4
handles.showDetect=0;
guidata(hObject, handles);

% --- Executes when selected object is changed in radioButtonGroupDetection.
function radioButtonGroupDetection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in radioButtonGroupDetection
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'detectionRadio1'
        detectionRadio1_Callback(hObject, eventdata, handles)
    case 'detectionRadio2'
        detectionRadio2_Callback(hObject, eventdata, handles)
    case 'detectionRadio3'
        detectionRadio3_Callback(hObject, eventdata, handles)
    case 'detectionRadio3'
        detectionRadio4_Callback(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function radioButtonGroupDetection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radioButtonGroupDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on mouse press over axes background.
function helpPic_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to helpPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open plusTipTrackViz_README.txt


% --- Executes during object creation, after setting all properties.
function helpPic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to helpPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate helpPic
img=imread('qIcon.jpg');
imagesc(img,'parent',hObject);
axis image
axis off
imHandle=get(hObject,'Children');
%info=get(hObject); info2=get(blah);
%set(blah,'HitTest','off')
set(imHandle,'ButtonDownFcn',@helpPic_ButtonDownFcn)


% --- Executes on button press in getQueryStr_Check.
function getQueryStr_Check_Callback(hObject, eventdata, handles)
% hObject    handle to getQueryStr_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getQueryStr_Check
handles.getStr=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in getHelpPush.
function getHelpPush_Callback(hObject, eventdata, handles)
% hObject    handle to getHelpPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open plusTipTrackViz_README.txt

% --- Executes during object creation, after setting all properties.
function getHelpPush_CreateFcn(hObject, eventdata, handles)
% hObject    handle to getHelpPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'CData',imread('help_icon.bmp'));


% --- Executes on button press in getProjListFile_check.
function getProjListFile_check_Callback(hObject, eventdata, handles)
% hObject    handle to getProjListFile_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getProjListFile_check
handles.loadProjList=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on selection change in xaxisScatterDrop.
function xaxisScatterDrop_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisScatterDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xaxisScatterDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xaxisScatterDrop
val = get(hObject,'Value');

% string labels in drop down menu
% 
% Growth speed
% Growth lifetime
% Growth displacement
% Shrinkage speed
% Shrinkage lifetime
% Shrinkage displacement
% Gap speed
% Gap lifetime
% Gap displacement

group = {...
    'growthSpeed',... 
    'growthLifetime',... 
    'growthDisp',... 
    'shrinkSpeed',...
    'shrinkLifetime',...
    'shrinkDisp'...
    'gapSpeed',...
    'gapLifetime',...
    'gapDisp',...
    };

handles.xaxisScatter=group{val};

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function xaxisScatterDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxisScatterDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xaxisScatterPerc_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisScatterPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxisScatterPerc as text
%        str2double(get(hObject,'String')) returns contents of xaxisScatterPerc as a double
handles.xScatterPercent = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function xaxisScatterPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxisScatterPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yaxisScatterPerc_Callback(hObject, eventdata, handles)
% hObject    handle to yaxisScatterPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaxisScatterPerc as text
%        str2double(get(hObject,'String')) returns contents of yaxisScatterPerc as a double
handles.yScatterPercent = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function yaxisScatterPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxisScatterPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yaxisScatterDrop.
function yaxisScatterDrop_Callback(hObject, eventdata, handles)
% hObject    handle to yaxisScatterDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns yaxisScatterDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yaxisScatterDrop
val = get(hObject,'Value');

% string labels in drop down menu
% 
% Growth speed
% Growth lifetime
% Growth displacement
% Shrinkage speed
% Shrinkage lifetime
% Shrinkage displacement
% Gap speed
% Gap lifetime
% Gap displacement

group = {...
    'growthSpeed',... 
    'growthLifetime',... 
    'growthDisp',... 
    'shrinkSpeed',...
    'shrinkLifetime',...
    'shrinkDisp'...
    'gapSpeed',...
    'gapLifetime',...
    'gapDisp',...
    };

handles.yaxisScatter=group{val};

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yaxisScatterDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxisScatterDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in quadScatterPlotPush.
function quadScatterPlotPush_Callback(hObject, eventdata, handles)
% hObject    handle to quadScatterPlotPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'projData')
    handles.projData=[];
end
param1=handles.xaxisScatter;
cutoff1=handles.xScatterPercent;
param2=handles.yaxisScatter;
cutoff2=handles.yScatterPercent;
[thresh1,thresh2]=plusTipParamPlot(param1,cutoff1,param2,cutoff2,handles.projData);

