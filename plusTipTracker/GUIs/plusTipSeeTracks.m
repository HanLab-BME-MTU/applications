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

% Last Modified by GUIDE v2.5 05-Oct-2011 09:02:45


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


%place image onto the axes, remove tick marks
pic=imread('pTT_logo_sm.jpg');
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
    selectionList=a(:,1);

    % allow multiple projects to be selected
    if ~isempty(selectionList)
        selection=listdlg('PromptString',...
            'Select the project you want to visualize','ListString',selectionList,...
            'SelectionMode','single','ListSize',[300 300]);
    end

    handles.dataDir=[];
    % if only one project was selected, save projData info and get data
    if isempty(selection)
        msgbox('No projects selected or tracking has not been completed.')
        handles.projList=[];
    else
        handles.projList=handles.projList(selection,1);
        handles.dataDir=selectionList{selection,1};
        p=load([handles.dataDir filesep 'meta' filesep 'projData.mat']);
        handles.projData=p.projData;
        set(handles.anDirEdit,'String',selectionList{selection,1});
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
velLimVal=get(handles.speedLimitEdit,'String');
if strcmpi(velLimVal,'max')
    velLimit=inf;
else
    velLimit=str2double(velLimVal);
end
        
%Read output directory
if isempty(handles.projData), errordlg('Please select a project'); return; end
saveDir = get(handles.edit_outputDir,'String');
if isempty(saveDir), errordlg('Please select an output directory'); return; end
handles.projData.saveDir = saveDir; % For compatibility with plusTipTrackMovie

plusTipSpeedMovie(handles.projData,handles.timeRangeDetect,velLimit,...
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
indivTrack=str2num(get(handles.indivTrackNumbersEdit,'String'))';
magCoef =[];

%Read output directory
if isempty(handles.projData), errordlg('Please select a project'); return;  end
saveDir = get(handles.edit_outputDir,'String');
if isempty(saveDir), errordlg('Please select an output directory');return;  end
handles.projData.saveDir = saveDir; % For compatibility with plusTipTrackMovie

plusTipTrackMovie(handles.projData,indivTrack,handles.timeRangeDetect,...
    handles.roi,magCoef,showTracks,showDetect,doAvi,rawToo);


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


% --- Executes on button press in selectOutputDirPush.
function selectOutputDirPush_Callback(hObject, eventdata, handles)

if isfield(handles,'dataDir')
    dirStart=handles.dataDir;
else
    dirStart=pwd;
end

outDir=uigetdir(dirStart,'Please select output directory for plots and movies');
if isequal(outDir,0), 
    set(handles.edit_outputDir,'String','');
end
set(handles.edit_outputDir,'String',outDir);
guidata(hObject, handles);


%---- Executes on button press in pushbutton_createMTdynamicsMaps.
function pushbutton_createMTdynamicsMaps_Callback(hObject, eventdata, handles)

%Read output directory
if isempty(handles.projData), errordlg('Please select a project'); return; end
saveDir = get(handles.edit_outputDir,'String');
if isempty(saveDir), errordlg('Please select an output directory'); return; end
handles.projData.saveDir = saveDir; % For compatibility with plusTipTrackMovie

remBegEnd=1;
% Read speed limit
value = get(handles.edit_speedLimMax,'String');
if strcmpi(value,'max'), speedLim=[]; else speedLim = str2double(value); end
if isnan(speedLim), errordlg('Please enter a valid maximum speed'); return; end

% Read lifetime limit
value = get(handles.edit_lifeLimMax,'String');
if strcmpi(value,'max'), lifeLim=[]; else lifeLim = str2double(value); end
if isnan(speedLim), errordlg('Please enter a valid maximum lifetime'); return; end


% Read displacement range
value = get(handles.edit_dispLimMax,'String');
if strcmpi(value,'max'), dispLim=[]; else dispLim = str2double(value); end
if isnan(dispLim), errordlg('Please enter a valid maximum displacement'); return; end

plusTipPlotResults(handles.projData,remBegEnd,handles.timeRangeDetect,...
    speedLim,lifeLim,dispLim,saveDir);
