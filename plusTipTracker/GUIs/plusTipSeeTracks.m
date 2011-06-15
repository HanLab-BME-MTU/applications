function varargout = plusTipSeeTracks(varargin)
% PLUSTIPSEETRACKS M-file for plusTipSeeTracks.fig
%      PLUSTIPSEETRACKS, by itself, creates a new PLUSTIPSEETRACKS or raises the
%      existing
%      singleton*.
%
%      H = PLUSTIPSEETRACKS returns the handle to a new PLUSTIPSEETRACKS or the handle to
%      the existing singleton*.
%
%      PLUSTIPSEETRACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPSEETRACKS.M with the given input arguments.
%
%      PLUSTIPSEETRACKS('Property','Value',...) creates a new PLUSTIPSEETRACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipSeeTracks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipSeeTracks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipSeeTracks

% Last Modified by GUIDE v2.5 11-Jun-2011 10:56:27


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plusTipSeeTracks_OpeningFcn, ...
    'gui_OutputFcn',  @plusTipSeeTracks_OutputFcn, ...
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


% --- Executes just before plusTipSeeTracks is made visible.
function plusTipSeeTracks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipSeeTracks (see VARARGIN)

% Choose default command line output for plusTipSeeTracks
handles.output = hObject;

% for "select projects" pushbutton
handles.projList=[]; % select projects pushbutton
handles.loadProjList = 0; % load projList checkbox
handles.getStr = 0; % narrow down list checkbox
handles.projData=[]; % if one project is selected, projData will be retrieved

% for "create groups" pushbutton
handles.autoGrp=1; % auto group from hierarchy
handles.groupList=[]; % also select groups pushbutton

% for "select saved ROI" pushbutton
handles.roi=[]; 

% for "choose frame range" edit boxes
handles.timeRangeDetect=[1 inf];

% check for movie type
handles.doAvi=0;

% for "track overlays" panel
handles.img=[];
handles.ask4select=0;
handles.selectedTracks=[];
handles.plotCurrentOnly=[];
handles.movieInfo=[];

% for "speed movies" panel
handles.velLimit=inf;

% for "track movies" panel
handles.indivTrack=[];
handles.magCoef=[];
handles.showTracks=1;
handles.showDetect=3;
handles.rawToo=1;

% for "sub-ROIs" panel
handles.subroiSelectType=0; % 0=manual, 1=center/singlePeriph, 2=center/quadPeriph
handles.subroiDistUnit='Microns';
handles.subroiDistVal=[];
handles.subroiTimeUnit='Fraction';
handles.subroiTimeVal=[];
handles.subroiExcludeRegion=0;

% for "quadrant scatter plots" panel
handles.xAxisScatParam='growthSpeed';
handles.xScatValType='Value';
handles.xScatInput=15;
handles.xAxisLim=[-inf inf];
handles.yAxisScatParam='growthLifetime';
handles.yScatValType='Value';
handles.yScatInput=10;
handles.yAxisLim=[-inf inf];
handles.remBegEnd=1;
handles.doBatchQuad=1;
handles.doPlotQuad=1;

%place image onto the axes, remove tick marks
pic=imread('pTT_logo_sm.png');
axes(handles.logoAxes);
image(pic);
axis off

set(handles.getHelpPush,'CData',imread('help_icon.png'),'Callback',...
    @(hObject,eventdata)open('plusTipSeeTracks.pdf'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plusTipSeeTracks wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = plusTipSeeTracks_OutputFcn(hObject, eventdata, handles)
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
    a=projList2Cell(handles.projList);
    a=a(:,1);

    % allow multiple projects to be selected
    if isempty(a)
        selection=[];
    else
        [selection,selectionList]=listSelectGUI(a,[],'move',1);
    end

    handles.dataDir=[];
    % if only one project was selected, save projData info and get data
    if isempty(selection)
        msgbox('No projects selected or tracking has not been completed.')
        handles.projList=[];
    elseif size(selection,1)==1
        handles.projList=handles.projList(selection,1);
        handles.dataDir=selectionList{1,1};
        p=load([handles.dataDir filesep 'meta' filesep 'projData.mat']);
        handles.projData=p.projData;
        set(handles.anDirEdit,'String',handles.projData.anDir);
        
        assignin('base','projData',handles.projData)
    else
        handles.projList=handles.projList(selection,1);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);
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
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFrameDetect');
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
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameDetect');
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


% --- Executes on button press in getQueryStr_Check.
function getQueryStr_Check_Callback(hObject, eventdata, handles)
% hObject    handle to getQueryStr_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getQueryStr_Check
handles.getStr=get(hObject,'Value');
guidata(hObject, handles);


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
% Growth speed (um/min)
% Growth lifetime (sec)
% Growth displacement (um)
% Bgap speed (um/min)
% Bgap lifetime (sec)
% Bgap displacement (um)
% Fgap speed (um/min)
% Fgap lifetime (sec)
% Fgap displacement (um)

group = {...
    'growthSpeed',...
    'growthLifetime',...
    'growthDisp',...
    'bgapSpeed',...
    'bgapLifetime',...
    'bgapDisp'...
    'fgapSpeed',...
    'fgapLifetime',...
    'fgapDisp',...
    };

handles.xAxisScatParam=group{val};

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



function xAxisScatterInput_Callback(hObject, eventdata, handles)
% hObject    handle to xAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xAxisScatterInput as text
%        str2double(get(hObject,'String')) returns contents of xAxisScatterInput as a double
handles.xScatInput = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xAxisScatterInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xAxisLim_Callback(hObject, eventdata, handles)
% hObject    handle to xAxisLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xAxisLim as text
%        str2double(get(hObject,'String')) returns contents of xAxisLim as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'xAxisLim');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xAxisLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xAxisLim (see GCBO)
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
% Growth speed (um/min)
% Growth lifetime (sec)
% Growth displacement (um)
% Bgap speed (um/min)
% Bgap lifetime (sec)
% Bgap displacement (um)
% Fgap speed (um/min)
% Fgap lifetime (sec)
% Fgap displacement (um)

group = {...
    'growthSpeed',...
    'growthLifetime',...
    'growthDisp',...
    'bgapSpeed',...
    'bgapLifetime',...
    'bgapDisp'...
    'fgapSpeed',...
    'fgapLifetime',...
    'fgapDisp',...
    };

handles.yAxisScatParam=group{val};

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


function yAxisScatterInput_Callback(hObject, eventdata, handles)
% hObject    handle to yAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yAxisScatterInput as text
%        str2double(get(hObject,'String')) returns contents of yAxisScatterInput as a double
handles.yScatInput = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function yAxisScatterInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yAxisLim_Callback(hObject, eventdata, handles)
% hObject    handle to yAxisLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yAxisLim as text
%        str2double(get(hObject,'String')) returns contents of yAxisLim as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'yAxisLim');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function yAxisLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yAxisLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in quadScatterPlotPush.
function quadScatterPlotPush_Callback(hObject, eventdata, handles)
% hObject    handle to quadScatterPlotPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'projData') && handles.doBatchQuad==0
    handles.projData.anDir=formatPath(handles.projData.anDir);
    saveDir=handles.projData.anDir;
    handles.groupList={'singleProject',handles.projData.anDir};
else
    saveDir=[];
end

xAxisInfo.name=handles.xAxisScatParam;
if strcmpi(handles.xScatValType,'percentile')
    xAxisInfo.splitPercentile=handles.xScatInput;
    xAxisInfo.splitValue=[];
else
    xAxisInfo.splitPercentile=[];
    xAxisInfo.splitValue=handles.xScatInput;
end
yAxisInfo.name=handles.yAxisScatParam;
if strcmpi(handles.yScatValType,'percentile')
    yAxisInfo.splitPercentile=handles.yScatInput;
    yAxisInfo.splitValue=[];
else
    yAxisInfo.splitPercentile=[];
    yAxisInfo.splitValue=handles.yScatInput;
end

xAxisInfo.minMax=handles.xAxisLim;
yAxisInfo.minMax=handles.yAxisLim;

plusTipQuadScatter(xAxisInfo,yAxisInfo,handles.groupList,handles.remBegEnd,...
    handles.timeRangeDetect,handles.doPlotQuad,saveDir);


% --- Executes on button press in selectOutputDirPush.
function selectOutputDirPush_Callback(hObject, eventdata, handles)
% hObject    handle to selectOutputDirPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectOutputDir');
guidata(hObject, handles);


% --- Executes on button press in dualPanelCheck.
function dualPanelCheck_Callback(hObject, eventdata, handles)
% hObject    handle to dualPanelCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dualPanelCheck
handles=plusTipGuiSwitch(hObject,eventdata,handles,'dualPanelCheck');
guidata(hObject, handles);


% --- Executes on button press in remTrackBegEnd.
function remTrackBegEnd_Callback(hObject, eventdata, handles)
% hObject    handle to remTrackBegEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of remTrackBegEnd
handles=plusTipGuiSwitch(hObject,eventdata,handles,'remBegEndCheck');
guidata(hObject, handles);


% --- Executes on selection change in xParamDrop.
function xParamDrop_Callback(hObject, eventdata, handles)
% hObject    handle to xParamDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xParamDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xParamDrop
val = get(hObject,'Value');
group = {'Value','Percentile'};
handles.xScatValType=group{val};

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xParamDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xParamDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yParamDrop.
function yParamDrop_Callback(hObject, eventdata, handles)
% hObject    handle to yParamDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns yParamDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yParamDrop
val = get(hObject,'Value');
group = {'Value','Percentile'};
handles.yScatValType=group{val};

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yParamDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yParamDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in subRoiPush.
function subRoiPush_Callback(hObject, eventdata, handles)
% hObject    handle to subRoiPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plusTipSubRoiTool(handles.projList,handles.subroiSelectType,...
    handles.subroiDistUnit,handles.subroiDistVal,...
    handles.subroiTimeUnit,handles.subroiTimeVal,...
    handles.roi,handles.subroiExcludeRegion);


function subroiDistValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to subroiDistValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subroiDistValEdit as text
%        str2double(get(hObject,'String')) returns contents of subroiDistValEdit as a double
handles.subroiDistVal=str2double(get(hObject,'String'));
if isnan(handles.subroiDistVal)
    handles.subroiDistVal=[];
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function subroiDistValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subroiDistValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function anDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to anDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anDirEdit as text
%        str2double(get(hObject,'String')) returns contents of anDirEdit as a double


% --- Executes during object creation, after setting all properties.
function anDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in summQuadPlotOnlyCheck.
function summQuadPlotOnlyCheck_Callback(hObject, eventdata, handles)
% hObject    handle to summQuadPlotOnlyCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of summQuadPlotOnlyCheck
handles.doPlotQuad=~get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in pickGroupsPush.
function pickGroupsPush_Callback(hObject, eventdata, handles)
% hObject    handle to pickGroupsPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'pickGroups');
assignin('base','groupList',handles.groupList);
guidata(hObject, handles);

% --- Executes on button press in autoGrpCheck.
function autoGrpCheck_Callback(hObject, eventdata, handles)
% hObject    handle to autoGrpCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoGrpCheck
handles.autoGrp=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in batchQuadCheck.
function batchQuadCheck_Callback(hObject, eventdata, handles)
% hObject    handle to batchQuadCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of batchQuadCheck
handles.doBatchQuad=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in selectGroupsPush.
function selectGroupsPush_Callback(hObject, eventdata, handles)
% hObject    handle to selectGroupsPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.groupList]=combineGroupListFiles(0);
assignin('base','groupList',handles.groupList);
guidata(hObject, handles);


% --- Executes on selection change in subroiDistUnitPop.
function subroiDistUnitPop_Callback(hObject, eventdata, handles)
% hObject    handle to subroiDistUnitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns subroiDistUnitPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subroiDistUnitPop
val = get(hObject,'Value');
group = {'Microns','Fraction'};
handles.subroiDistUnit=group{val};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subroiDistUnitPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subroiDistUnitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in subroiRadioPanel.
function subroiRadioPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in subroiRadioPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'subManualRadio'
        subManualRadio_Callback(hObject, eventdata, handles)
    case 'subAutoRadio'
        subAutoRadio_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in subroiManualRadio.
function subroiManualRadio_Callback(hObject, eventdata, handles)
% hObject    handle to subroiManualRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subroiManualRadio
handles.subroiSelectType=0;
if handles.subroiSelectType==0
    set(handles.subroiDistValEdit,'Enable','Off');
    set(handles.subroiDistValEdit,'String','');
    handles.subroiDistVal=[];
    
    set(handles.subroiDistUnitPop,'Enable','Off');
    set(handles.subroiAutoDivPeriphCheck,'Enable','Off');
else
    set(handles.subroiDistValEdit,'Enable','On');
    set(handles.subroiDistUnitPop,'Enable','On');
    set(handles.subroiAutoDivPeriphCheck,'Enable','On');
end
guidata(hObject, handles);


% --- Executes on button press in subroiAutoRadio.
function subroiAutoRadio_Callback(hObject, eventdata, handles)
% hObject    handle to subroiAutoRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subroiAutoRadio
handles.subroiSelectType=1;
if handles.subroiSelectType==0
    set(handles.subroiDistValEdit,'Enable','Off');
    set(handles.subroiDistValEdit,'String','');
    handles.subroiDistVal=[];
    
    set(handles.subroiDistUnitPop,'Enable','Off');
    set(handles.subroiAutoDivPeriphCheck,'Enable','Off');
else
    set(handles.subroiDistValEdit,'Enable','On');
    set(handles.subroiDistUnitPop,'Enable','On');
    set(handles.subroiAutoDivPeriphCheck,'Enable','On');
end

guidata(hObject, handles);


% --- Executes on button press in subroiAutoDivPeriphCheck.
function subroiAutoDivPeriphCheck_Callback(hObject, eventdata, handles)
% hObject    handle to subroiAutoDivPeriphCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subroiAutoDivPeriphCheck
handles.subroiSelectType=2;
guidata(hObject, handles);



function subroiTimeValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to subroiTimeValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subroiTimeValEdit as text
%        str2double(get(hObject,'String')) returns contents of subroiTimeValEdit as a double
handles.subroiTimeVal=str2double(get(hObject,'String'));
if isnan(handles.subroiTimeVal)
    handles.subroiTimeVal=[];
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subroiTimeValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subroiTimeValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in subroiTimeUnitPop.
function subroiTimeUnitPop_Callback(hObject, eventdata, handles)
% hObject    handle to subroiTimeUnitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns subroiTimeUnitPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subroiTimeUnitPop
val = get(hObject,'Value');
group = {'Fraction','Seconds'};
handles.subroiTimeUnit=group{val};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subroiTimeUnitPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subroiTimeUnitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in subroiExcludeCheck.
function subroiExcludeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to subroiExcludeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subroiExcludeCheck
handles.subroiExcludeRegion=get(hObject,'Value');
guidata(hObject, handles);
