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

% Last Modified by GUIDE v2.5 17-Oct-2011 11:18:36


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
handles.projData=[]; % if one project is selected, projData will be retrieved

% for "create groups" pushbutton
handles.groupList=[]; % also select groups pushbutton

% for "select saved ROI" pushbutton
handles.roi=[]; 

% for "track overlays" panel
handles.selectedTracks=[];

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
        set(handles.edit_currentProject,'String',handles.dataDir);
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

[FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat file');
if ~isequal(FileName,0)
    if ~isempty(strfind(FileName,'tif'))
        handles.roi=imread([PathName FileName]);
    elseif ~isempty(strfind(FileName,'mat'))
        p=load([PathName FileName]);
        handles.roi=p.roiYX;
    else
        errordlg('File chosen was not a roiYX.mat or roiMask.tif file. Please try again.','File Error');
    end
end
guidata(hObject, handles);


% --- Executes on button press in plotTracksPush.
function plotTracksPush_Callback(hObject, eventdata, handles)
% hObject    handle to plotTracksPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.projData), errordlg('Please select a project'); return; end

if isempty(handles.selectedTracks), state ='off'; else state='on'; end
set(handles.selectedTracksDisplay,'Enable',state,'String',num2str(handles.selectedTracks));

% Read time range
timeRange = getTimeRange(handles);
rawToo=0;
isStill=1;
img=[];
plotCurrentOnly=[];
movieInfo=[];
ask4select = get(handles.selectTracksCheck,'Value');
[handles.selectedTracks] = plusTipPlotTracks(handles.projData,[],...
    timeRange,img,ask4select,...
    plotCurrentOnly,handles.roi,movieInfo,rawToo,isStill);

if ~isempty(handles.selectedTracks)
    
    temp=vertcat(handles.selectedTracks{:});
    handles.selectedTracks=unique(temp(:,1));
    [l w]=size(handles.selectedTracks);
    
    if l>w
        handles.selectedTracks=handles.selectedTracks';
    end
    
end
        
if isempty(handles.selectedTracks), state ='off'; else state='on'; end
set(handles.selectedTracksDisplay,'Enable',state,'String',num2str(handles.selectedTracks));

guidata(hObject, handles);

% --- Executes on button press in speedMovieButton.
function speedMovieButton_Callback(hObject, eventdata, handles)

%Read project and output directory
if isempty(handles.projData), errordlg('Please select a project'); return; end
saveDir = get(handles.edit_outputDir,'String');
if isempty(saveDir), errordlg('Please select an output directory'); return; end
handles.projData.saveDir = saveDir; % For compatibility with plusTipTrackMovie

% Read time range
timeRange = getTimeRange(handles);

% Read movie specific parameters
doAvi=get(handles.aviCheckTrackMov,'Value');
velLimVal=get(handles.speedLimitEdit,'String');
if strcmpi(velLimVal,'max')
    velLimit=inf;
else
    velLimit=str2double(velLimVal);
end
        
plusTipSpeedMovie(handles.projData,timeRange,velLimit,handles.roi,doAvi);

% --- Executes on button press in trackMovieButton.
function trackMovieButton_Callback(hObject, eventdata, handles)

%Read project and output directory
if isempty(handles.projData), errordlg('Please select a project'); return; end
saveDir = get(handles.edit_outputDir,'String');
if isempty(saveDir), errordlg('Please select an output directory'); return; end
handles.projData.saveDir = saveDir; % For compatibility with plusTipTrackMovie

% Read time range
timeRange = getTimeRange(handles);

% Read detection mode
detectedObject= get(get(handles.radioButtonGroupDetection,'SelectedObject'),'Tag');
showDetect = str2double(detectedObject(length('detectionRadio')+1:end));

% Read movie parameters
doAvi=get(handles.aviCheckTrackMov,'Value');
rawToo=get(handles.dualPanelCheck,'Value');
showTracks=get(handles.showTracksCheck,'Value');
indivTrack=str2num(get(handles.indivTrackNumbersEdit,'String'))';
magCoef =[];

plusTipTrackMovie(handles.projData,indivTrack,timeRange,...
    handles.roi,magCoef,showTracks,showDetect,doAvi,rawToo);

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
 

timeRange = getTimeRange(handles);
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

plusTipPlotResults(handles.projData,remBegEnd,timeRange,...
    speedLim,lifeLim,dispLim,saveDir);


function timeRange = getTimeRange(handles)

% Read time range
sFVal=get(handles.startFrame,'String');
if strcmpi(sFVal,'min')
    timeRange(1)=1;
else
    timeRange(1)=str2double(sFVal);
end
eFVal=get(handles.endFrame,'String');
if strcmpi(eFVal,'max')
    timeRange(2)=Inf;
else
    timeRange(2)=str2double(eFVal);
end
if ~all(isposint(timeRange)),
    errordlg('Invalid frame  range. Please check again','Missing Input');
    return
end
