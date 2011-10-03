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

% Last Modified by GUIDE v2.5 03-Oct-2011 13:08:50


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
handles.showDetect=3;

% for "sub-ROIs" panel
set(handles.subroiDistUnitPop,'String',{'Microns','Fraction'},'Value',1);
set(handles.subroiTimeUnitPop,'String',{'Fraction','Seconds'},'Value',1);

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


function endFrame_Callback(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrame as text
%        str2double(get(hObject,'String')) returns contents of endFrame as
%        a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameDetect');
guidata(hObject, handles);

% --- Executes on button press in selectTracksCheck.
function selectTracksCheck_Callback(hObject, eventdata, handles)
% hObject    handle to selectTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectTracksCheck
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectTracksCheck');
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


function indivTrackNumbersEdit_Callback(hObject, eventdata, handles)
% hObject    handle to indivTrackNumbersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of indivTrackNumbersEdit as text
%        str2double(get(hObject,'String')) returns contents of
%        indivTrackNumbersEdit as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'indivTrackNumbersEdit');
guidata(hObject, handles);


% --- Executes on button press in plotTracksPush.
function plotTracksPush_Callback(hObject, eventdata, handles)
% hObject    handle to plotTracksPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectedTracksDisplay');
handles=plusTipGuiSwitch(hObject,eventdata,handles,'plotTracksPush');
handles=plusTipGuiSwitch(hObject,eventdata,handles,'selectedTracksDisplay');
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

% --- Executes on button press in speedMovieButton.
function speedMovieButton_Callback(hObject, eventdata, handles)

doAvi=get(handles.aviCheckTrackMov,'Value');

plusTipSpeedMovie(handles.projData,handles.timeRangeDetect,handles.velLimit,...
    handles.roi,doAvi);


% --- Executes on button press in trackMovieButton.
function trackMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to trackMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read detection mode
detectedObject= get(get(handles.radioButtonGroupDetection,'SelectedObject'),'Tag');
showDetect = str2double(detectedObject(length('detectionRadio')+1:end));

% Read movie parameters
doAvi=get(handles.aviCheckTrackMov,'Value');
rawToo=get(handles.dualPanelCheck,'Value');
showTracks=get(handles.showTracksCheck,'Value');
magCoef =[];
plusTipTrackMovie(handles.projData,handles.indivTrack,handles.timeRangeDetect,...
    handles.roi,magCoef,showTracks,showDetect,doAvi,rawToo);

% --- Executes when selected object is changed in radioButtonGroupDetection.
function radioButtonGroupDetection_SelectionChangeFcn(hObject, eventdata, handles)

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


function xAxisScatterInput_Callback(hObject, eventdata, handles)
% hObject    handle to xAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xAxisScatterInput as text
%        str2double(get(hObject,'String')) returns contents of xAxisScatterInput as a double
handles.xScatInput = str2double(get(hObject,'String'));
guidata(hObject, handles);

function xAxisLim_Callback(hObject, eventdata, handles)
% hObject    handle to xAxisLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xAxisLim as text
%        str2double(get(hObject,'String')) returns contents of xAxisLim as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'xAxisLim');
guidata(hObject, handles);


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


function yAxisScatterInput_Callback(hObject, eventdata, handles)
% hObject    handle to yAxisScatterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yAxisScatterInput as text
%        str2double(get(hObject,'String')) returns contents of yAxisScatterInput as a double
handles.yScatInput = str2double(get(hObject,'String'));
guidata(hObject, handles);


function yAxisLim_Callback(hObject, eventdata, handles)
% hObject    handle to yAxisLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yAxisLim as text
%        str2double(get(hObject,'String')) returns contents of yAxisLim as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'yAxisLim');
guidata(hObject, handles);


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

% --- Executes on button press in subRoiPush.
function subRoiPush_Callback(hObject, eventdata, handles)

% Read selection type
if get(handles.subroiManualRadio,'Value'),
    subroiSelectType =0;
elseif get(handles.subroiAutoDivPeriphCheck,'Value'),
    subroiSelectType =2;
else
    subroiSelectType =1;
end
      
% Read distance units and value
props = get(handles.subroiDistUnitPop,{'String','Value'});
subroiDistUnit=props{1}{props{2}};
subroiDistVal=str2double(get(handles.subroiDistValEdit,'String'));
if isnan(subroiDistVal), subroiDistVal=[]; end

% Read time units and value
props = get(handles.subroiTimeUnitPop,{'String','Value'});
subroiTimeUnit=props{1}{props{2}};
subroiTimeVal=str2double(get(handles.subroiTimeValEdit,'String'));
if isnan(subroiTimeVal), subroiTimeVal=[]; end

subroiExcludeRegion = get(handles.subroiExcludeCheck,'Value');

plusTipSubRoiTool(handles.projList,subroiSelectType,...
    subroiDistUnit,subroiDistVal,subroiTimeUnit,subroiTimeVal,...
    handles.roi,subroiExcludeRegion);

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

% --- Executes when selected object is changed in subroiRadioPanel.
function subroiRadioPanel_SelectionChangeFcn(hObject, eventdata, handles)

if strcmp(get(hObject,'Tag'),'subroiManualRadio')
    set(handles.subroiDistValEdit,'Enable','Off','String','');
    set(handles.subroiDistUnitPop,'Enable','Off');
    set(handles.subroiAutoDivPeriphCheck,'Enable','Off');
else
        set(handles.subroiDistValEdit,'Enable','On');
    set(handles.subroiDistUnitPop,'Enable','On');
    set(handles.subroiAutoDivPeriphCheck,'Enable','On');
end

guidata(hObject, handles);
  
