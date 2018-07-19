function varargout = plusTipMtVifDynamicAnalysisGUI(varargin) 
% PLUSTIPMTVIFDYNAMICANALYSISGUI M-file for plusTipMtVifDynamicAnalysisGUI.fig
%      PLUSTIPMTVIFDYNAMICANALYSISGUI, by itself, creates a new PLUSTIPMTVIFDYNAMICANALYSISGUI or raises the
%      existing
%      singleton*.
%
%      H = PLUSTIPMTVIFDYNAMICANALYSISGUI returns the handle to a new PLUSTIPMTVIFDYNAMICANALYSISGUI or the handle to
%      the existing singleton*.
%
%      PLUSTIPMTVIFDYNAMICANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPMTVIFDYNAMICANALYSISGUI.M with the given input arguments.
%
%      PLUSTIPMTVIFDYNAMICANALYSISGUI('Property','Value',...) creates a new PLUSTIPMTVIFDYNAMICANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipMtVifDynamicAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipMtVifDynamicAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipMtVifDynamicAnalysisGUI

% Last Modified by GUIDE v2.5 27-Jan-2012 15:48:55


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plusTipMtVifDynamicAnalysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @plusTipMtVifDynamicAnalysisGUI_OutputFcn, ...
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


% --- Executes just before plusTipMtVifDynamicAnalysisGUI is made visible.
function plusTipMtVifDynamicAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipMtVifDynamicAnalysisGUI (see VARARGIN)

% Choose default command line output for plusTipMtVifDynamicAnalysisGUI
handles.output = hObject;

% for "select projects" pushbutton
handles.Channels = cell(1,2);
handles.Channels{1} = struct('projList',[],'projData',[],'dataDir',[],'edit_currentProject',[],'MD',[]);
handles.Channels{2} = struct('projList',[],'projData',[],'dataDir',[],'edit_currentProject',[],'MD',[]);

handles.projList= []; % select projects pushbutton
projData = []; % if one project is selected, projData will be retrieved
handles.Channels{1}.projList = [];
handles.Channels{1}.projData = [];
handles.Channels{1}.MD = [];
handles.Channels{2}.projList = [];
handles.Channels{2}.projData = [];
handles.Channels{2}.MD = [];

% for "create groups" pushbutton
handles.groupList=[]; % also select groups pushbutton
handles.Channels{1}.groupList=[];
handles.Channels{2}.groupList=[];

% for "select saved ROI" pushbutton
handles.roi=[]; 
handles.Channels{1}.roi=[]; 
handles.Channels{2}.roi=[]; 

% for "track overlays" panel
handles.selectedTracks=[];
handles.Channels{1}.selectedTracks=[];
handles.Channels{2}.selectedTracks=[];

%place image onto the axes, remove tick marks
pic=imread('pTT_logo_sm.png');
axes(handles.logoAxes);
image(pic);
axis off

set(handles.getHelpPush,'CData',imread('help_icon.png'),'Callback',...
    @(hObject,eventdata)open('plusTipSeeTracks.pdf'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plusTipMtVifDynamicAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = plusTipMtVifDynamicAnalysisGUI_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in select_mt_channel.
function select_mt_channel_Callback(hObject, eventdata, handles)
% hObject    handle to select_mt_channel (see GCBO)
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
        set(handles.edit_Channel_1_proj,'String',handles.dataDir);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);

handles.Channels{1}.projList = handles.projList;
handles.Channels{1}.dataDir = handles.dataDir;
handles.Channels{1}.projData = handles.projData;
handles.Channels{1}.edit_currentProject = handles.edit_Channel_1_proj;

guidata(hObject, handles);


% --- Executes on button press in selectSavedRoi_MT_Pushbutton.
function selectSavedRoi_MT_Pushbutton_Callback(hObject, eventdata, handles)

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

handles.Channels{1}.roi = handles.roi;

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
% doAvi=get(handles.aviCheckTrackMov,'Value');
% velLimVal=get(handles.speedLimitEdit,'String');
% if strcmpi(velLimVal,'max')
%     velLimit=inf;
% else
%     velLimit=str2double(velLimVal);
% end
doAvi=0;
velLimit=Inf;

% Read red channel index
value = get(handles.edit_displayred,'String');
if strcmpi(value,'none'), handles.projData.redChannelInd=0; else handles.projData.redChannelInd = str2double(value); end
if isnan(handles.projData.redChannelInd), errordlg('Please enter a valid red channel index'); return; end

% Read red channel index
value = get(handles.edit_displaygreen,'String');
if strcmpi(value,'none'), handles.projData.greenChannelInd=0; else handles.projData.greenChannelInd = str2double(value); end
if isnan(handles.projData.greenChannelInd), errordlg('Please enter a valid green channel index'); return; end

% Read red channel index
value = get(handles.edit_displayblue,'String');
if strcmpi(value,'none'), handles.projData.blueChannelInd=0; else handles.projData.blueChannelInd = str2double(value); end
if isnan(handles.projData.blueChannelInd), errordlg('Please enter a valid blue channel index'); return; end

if (handles.projData.greenChannelInd==0&&handles.projData.blueChannelInd==0)
    handles.projData.blueChannelInd = handles.projData.redChannelInd;
    handles.projData.greenChannelInd = handles.projData.redChannelInd;
end

if (handles.projData.redChannelInd==0&&handles.projData.blueChannelInd==0)
    handles.projData.blueChannelInd =handles.projData.greenChannelInd;
    handles.projData.redChannelInd = handles.projData.greenChannelInd;
end

% Read detection mode
displayFeatureChannelinput= get(get(handles.radioDisplayPlusTipGroup,'SelectedObject'),'Tag');
displayFeatureChannelNo = str2double(displayFeatureChannelinput(length('radio_displayFeauture')+1:end));
displayFeatureChannelNo = int8(displayFeatureChannelNo);

handles.projData.DisplayFeatureChannel(1)= displayFeatureChannelNo/10;
handles.projData.DisplayFeatureChannel(2)= mod(displayFeatureChannelNo,10);

guidata(hObject, handles);

% handles.projData = handles.Channels{1}.projData;
% handles.projData.imDir = handles.Channels{2}.projData.imDir;
% handles.projData.anDir = handles.Channels{2}.projData.anDir;
% handles.projList = handles.Channels{2}.projList;

plusTipDualchannelTrackingMovie(handles.projData, handles.Channels, timeRange,velLimit,handles.roi,doAvi);

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
    handles.projDataroi,magCoef,showTracks,showDetect,doAvi,rawToo);

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


%---- Executes on button press in pushbutton_setdisplay.
function pushbutton_setdisplay_Callback(hObject, eventdata, handles)

% Read red channel index
value = get(handles.edit_displayred,'String');
if strcmpi(value,'none'), handles.projData.redChannelInd=0; else handles.projData.redChannelInd = str2double(value); end
if isnan(handles.projData.redChannelInd), errordlg('Please enter a valid red channel index'); return; end

% Read red channel index
value = get(handles.edit_displaygreen,'String');
if strcmpi(value,'none'), handles.projData.greenChannelInd=0; else handles.projData.greenChannelInd = str2double(value); end
if isnan(handles.projData.greenChannelInd), errordlg('Please enter a valid green channel index'); return; end

% Read red channel index
value = get(handles.edit_displayblue,'String');
if strcmpi(value,'none'), handles.projData.blueChannelInd=0; else handles.projData.blueChannelInd = str2double(value); end
if isnan(handles.projData.blueChannelInd), errordlg('Please enter a valid blue channel index'); return; end

% Read detection mode
displayFeatureChannelinput= get(get(handles.radioDisplayPlusTipGroup,'SelectedObject'),'Tag');
displayFeatureChannelNo = str2double(displayFeatureChannelinput(length('radio_displayFeauture')+1:end));
displayFeatureChannelNo = int8(displayFeatureChannelNo);

handles.projData.DisplayFeatureChannel(1)= displayFeatureChannelNo/10;
handles.projData.DisplayFeatureChannel(2)= mod(displayFeatureChannelNo,10);

guidata(hObject, handles);


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



function edit_Channel_2_proj_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Channel_2_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Channel_2_proj as text
%        str2double(get(hObject,'String')) returns contents of edit_Channel_2_proj as a double


% --- Executes during object creation, after setting all properties.
function edit_Channel_2_proj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Channel_2_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_channel_VIF.
function pushbutton_select_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to select_channel_VIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.projData2), errordlg('Please select a project'); return; end

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

% --- Executes on button press in selectSavedRoi_VIF_Pushbutton.
function selectSavedRoi_VIF_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectSavedRoi_VIF_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

handles.Channels{2}.roi = handles.roi;

guidata(hObject, handles);


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in select_channel_VIF.
function select_channel_VIF_Callback(hObject, eventdata, handles)
% hObject    handle to select_channel_VIF (see GCBO)
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
        set(handles.edit_Channel_2_proj,'String',handles.dataDir);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);

handles.Channels{2}.projList = handles.projList;
handles.Channels{2}.dataDir = handles.dataDir;
handles.Channels{2}.projData = handles.projData;
handles.Channels{2}.edit_currentProject = handles.edit_Channel_2_proj;
handles.Channels{2}.MD = load([handles.dataDir(1:end-6),'/movieData.mat']);

guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over select_channel_VIF.
function select_channel_VIF_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to select_channel_VIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_displayred_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displayred (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displayred as text
%        str2double(get(hObject,'String')) returns contents of edit_displayred as a double



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_displayredfeature.
function radiobutton_displayredfeature_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_displayredfeature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_displayredfeature


% --- Executes on button press in radiobutton_displaygreenfeature.
function radiobutton_displaygreenfeature_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_displaygreenfeature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_displaygreenfeature


% --- Executes on button press in pushbutton_PdfSlaveGivenMaster.
function pushbutton_PdfSlaveGivenMaster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PdfSlaveGivenMaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plusTipPdfSlaveGivenMaster(handles.projData, handles.Channels);


% --- Executes on button press in pushbutton_setmaster.
function pushbutton_setmaster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setmaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Read detection mode
displayMasterChannelTag = get(get(handles.radioWhichIsMasterGroup,'SelectedObject'),'Tag');
displayMasterChannelNo = str2num(displayMasterChannelTag(length('radioMaster')+1:end));

handles.projData.displayMasterChannel= displayMasterChannelNo;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton_PdfSlaveGivenMaster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_PdfSlaveGivenMaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_VIF_Direction_map.
function pushbutton_VIF_Direction_map_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_VIF_Direction_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

steerable_filter_segmentation(handles.Channels{2}.MD);

handles.track_direction_profiles = get_direction_profile(handles.Channels{1}.projData,handles.Channels{2}.MD);
