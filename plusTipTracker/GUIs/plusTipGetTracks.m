function varargout = plusTipGetTracks(varargin)
% PLUSTIPGETTRACKS M-file for plusTipGetTracks.fig
%      PLUSTIPGETTRACKS, by itself, creates a new PLUSTIPGETTRACKS or raises the existing
%      singleton*.
%
%      H = PLUSTIPGETTRACKS returns the handle to a new PLUSTIPGETTRACKS or the handle to
%      the existing singleton*.
%
%      PLUSTIPGETTRACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPGETTRACKS.M with the given input arguments.
%
%      PLUSTIPGETTRACKS('Property','Value',...) creates a new PLUSTIPGETTRACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipGetTracks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipGetTracks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipGetTracks
% 
%
% adding space to test SVN
%
% Last Modified by GUIDE v2.5 11-Jun-2011 09:59:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipGetTracks_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipGetTracks_OutputFcn, ...
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


% --- Executes just before plusTipGetTracks is made visible.
function plusTipGetTracks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipGetTracks (see VARARGIN)

% Choose default command line output for plusTipGetTracks
handles.output = hObject;
handles.getStr = 0;
handles.loadProjList = 0;
handles.projList = [];

handles.doDetect=0;
handles.doTrack=0; 
handles.doMeta=0; 

% DETECTION parameters
handles.timeRangeDetect = [1 inf];
handles.bitDepth  = 16; 
handles.savePlots = 1;


% TRACKING parameters
handles.timeWindow=[];
handles.minTrackLen=3;
handles.minRadius=[];
handles.maxRadius=[];
handles.maxFAngle=30;
handles.maxBAngle=10;
handles.maxShrinkFactor=1.5;
handles.fluctRad=1;
handles.timeRangeTrack = [1 inf];

% META parameters
handles.secPerFrame=[];
handles.pixSizeNm=[];
handles.timeRangePost=[1 inf];
handles.doHist=1;

%place image onto the axes, remove tick marks
set(handles.figure1,'CurrentAxes',handles.logoAxes);
image(imread('pTT_logo_sm.png'));
axis off

for i=1:3
    hslider= handles.(['slider_' num2str(i)]);
    props = get(hslider,{'Min','Max'});
    set(hslider,'SliderStep',[.1 .5]/(props{2}-props{1}));
end

updateDetection(hObject, eventdata, handles);

set(handles.getHelpPush,'CData',imread('help_icon.png'),'Callback',...
    @(hObject,eventdata)open('plusTipGetTracks.pdf'));

userData= get(handles.figure1,'UserData');
userData.previewGUI=-1;
% Update handles structure
set(handles.figure1,'UserData',userData)
guidata(hObject,handles);

% UIWAIT makes plusTipGetTracks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipGetTracks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in getProjPush.
function getProjPush_Callback(hObject, eventdata, handles)
% hObject    handle to getProjPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plusTipGuiSwitch(hObject,eventdata,handles,'getProjPush');

if ~isempty(handles.projList)
    
    % here we filter out any sub-directories
    a=projList2Cell(handles.projList);
    a=a(:,1);
    b=cellfun(@isempty, strfind(a,'sub'));
    a=a(b);
        
    % allow multiple projects to be selected
    [selection,selectionList]=listSelectGUI(a,[],'move',1);

    if ~isempty(selection)
        handles.projList=handles.projList(selection,1);
    else
        msgbox('No projects selected.')
        handles.projList=[];
    end
else
    msgbox('No projects selected.')
    handles.projList=[];
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);
guidata(hObject, handles);

if ~isempty(handles.projList)
    [file,path] = uiputfile('projList.mat','Find a location to save your project list',...
        [pwd filesep 'projList.mat']);
    if isequal(file,0), return; end
    save([path file],'-struct','handles','projList');
end

% --- Executes on button press in getQueryStr_Check.
function getQueryStr_Check_Callback(hObject, eventdata, handles)
% hObject    handle to getQueryStr_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.getStr=get(hObject,'Value');
guidata(hObject, handles);


function startFrameDetect_Callback(hObject, eventdata, handles)
% hObject    handle to startFrameDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrameDetect as text
%        str2double(get(hObject,'String')) returns contents of startFrameDetect as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFrameDetect');
guidata(hObject, handles);


function endFrameDetect_Callback(hObject, eventdata, handles)
% hObject    handle to endFrameDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrameDetect as text
%        str2double(get(hObject,'String')) returns contents of endFrameDetect as
%        a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameDetect');  
guidata(hObject, handles);


function startFramePost_Callback(hObject, eventdata, handles)
% hObject    handle to startFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFramePost as text
%        str2double(get(hObject,'String')) returns contents of startFramePost as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFramePost');
guidata(hObject, handles);

function endFramePost_Callback(hObject, eventdata, handles)
% hObject    handle to endFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFramePost as text
%        str2double(get(hObject,'String')) returns contents of endFramePost as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFramePost');
guidata(hObject, handles);

function startFrameTrack_Callback(hObject, eventdata, handles)
% hObject    handle to startFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrameTrack as text
%        str2double(get(hObject,'String')) returns contents of startFrameTrack as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFrameTrack');
guidata(hObject, handles);


function endFrameTrack_Callback(hObject, eventdata, handles)
% hObject    handle to endFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrameTrack as text
%        str2double(get(hObject,'String')) returns contents of endFrameTrack as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameTrack');
guidata(hObject, handles);


function bitDepthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bitDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bitDepthEdit as text
%        str2double(get(hObject,'String')) returns contents of bitDepthEdit as a double
handles.bitDepth=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes on button press in savePlotCheck.
function savePlotCheck_Callback(hObject, eventdata, handles)
% hObject    handle to savePlotCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savePlotCheck
handles.savePlots=get(hObject,'Value');
guidata(hObject, handles);

function frameRateEdit_Callback(hObject, eventdata, handles)
% hObject    handle to frameRateEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameRateEdit as text
%        str2double(get(hObject,'String')) returns contents of frameRateEdit as a double
handles.secPerFrame=str2double(get(hObject,'String'));
guidata(hObject, handles);

function pixSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pixSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of pixSizeEdit as a double
handles.pixSizeNm=str2double(get(hObject,'String'));
guidata(hObject, handles);


function timeWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to timeWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of timeWindowEdit as a double
handles.timeWindow=str2double(get(hObject,'String'));
guidata(hObject, handles);

function minTrackLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minTrackLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minTrackLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of minTrackLengthEdit as a double
handles.minTrackLen=str2double(get(hObject,'String'));
guidata(hObject, handles);

function minRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of minRadiusEdit as a double
handles.minRadius=str2double(get(hObject,'String'));
guidata(hObject, handles);

function maxRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of maxRadiusEdit as a double
handles.maxRadius=str2double(get(hObject,'String'));
guidata(hObject, handles);


function maxFAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxFAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxFAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of maxFAngleEdit as a double
handles.maxFAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);


function maxBAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxBAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxBAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of maxBAngleEdit as a double
handles.maxBAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);


function maxShrinkFactorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxShrinkFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxShrinkFactorEdit as text
%        str2double(get(hObject,'String')) returns contents of maxShrinkFactorEdit as a double
handles.maxShrinkFactor=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes on button press in doDetectCheck.
function updateDetection(hObject, eventdata, handles)
% hObject    handle to doDetectCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doDetectCheck
handles.doDetect=get(handles.doDetectCheck,'Value');
detectionChildren=get(handles.detectionPanel,'Children');
isPanel = strcmp(get(detectionChildren,'Type'),'uipanel');
panelTags = get(detectionChildren(isPanel),'Tag');
panelIndx = cellfun(@(x) str2double(x(length('uipanel_method')+1:end)),...
    panelTags);
methodIndx = get(handles.popupmenu_detectionMethod,'Value');
methodPanelTag = panelTags{panelIndx==methodIndx};
set(handles.(methodPanelTag),'Visible','on');
otherPanelTag = panelTags(panelIndx~=methodIndx);
cellfun(@(x) set(handles.(x),'Visible','off'),otherPanelTag);
    
if handles.doDetect==1
    set(detectionChildren(~isPanel),'Enable','on');
    set(handles.startFrameDetect,'Enable','on');
    set(handles.endFrameDetect,'Enable','on');    
    if get(handles.checkbox_custom,'Value'),enableState = 'on'; else enableState = 'off'; end
    set(get(handles.(methodPanelTag),'Children'),'Enable',enableState);
    otherPanelTag = panelTags(panelIndx~=methodIndx);
    cellfun(@(x) set(handles.(x),'Visible','off'),otherPanelTag);
    cellfun(@(x) set(get(handles.(x),'Children'),'Enable','off'),otherPanelTag);
else
    set(detectionChildren(~isPanel),'Enable','off');
    arrayfun(@(x) set(get(x,'Children'),'Enable','off'),...
        detectionChildren(isPanel));
    set(handles.startFrameDetect,'Enable','off');
    set(handles.endFrameDetect,'Enable','off');
end
guidata(hObject, handles);


% --- Executes on button press in trackingCheck.
function trackingCheck_Callback(hObject, eventdata, handles)
% hObject    handle to trackingCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trackingCheck
handles.doTrack=get(hObject,'Value');
trackingHandles=get(handles.trackingPanel,'Children');
if handles.doTrack==1
    set(trackingHandles,'Enable','on')
    set(handles.startFrameTrack,'Enable','on')
    set(handles.endFrameTrack,'Enable','on')
    
else
    set(trackingHandles,'Enable','off')
    set(handles.startFrameTrack,'Enable','off')
    set(handles.endFrameTrack,'Enable','off')
end
guidata(hObject, handles);


% --- Executes on button press in metaCheck.
function metaCheck_Callback(hObject, eventdata, handles)
% hObject    handle to metaCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of metaCheck
handles.doMeta=get(hObject,'Value');
if handles.doMeta==1
    set(handles.startFramePost,'Enable','on')
    set(handles.endFramePost,'Enable','on')
    set(handles.frameRateEdit,'Enable','on')
    set(handles.pixSizeEdit,'Enable','on')
    set(handles.histCheck,'Enable','on')
else
    set(handles.startFramePost,'Enable','off')
    set(handles.endFramePost,'Enable','off')
    set(handles.frameRateEdit,'Enable','off')
    set(handles.pixSizeEdit,'Enable','off')
    set(handles.histCheck,'Enable','off')
end
guidata(hObject, handles);


% --- Executes on button press in startPush.
function startPush_Callback(hObject, eventdata, handles)
% hObject    handle to startPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.projList)
    errordlg('Please select project(s) first','No project error');
    return
end

numProj=size(handles.projList,1);

if handles.doTrack==1
    if isempty(handles.timeWindow) || isempty(handles.minRadius) || isempty(handles.maxRadius)
        errordlg('Please check tracking parameters.','Missing Input');
        return
    end       
end
if handles.doMeta==1
    if isempty(handles.secPerFrame) || isempty(handles.pixSizeNm)
        errordlg('Please check post-processing parameters.','Missing Input');
        return
    end       
end


for i=1:numProj
    %try
    % detection
    if handles.doDetect==1
        tic
        disp(['Detecting project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir]);
        detectionMethod = get(handles.popupmenu_detectionMethod,'Value');
        switch detectionMethod
            case 1
                if get(handles.checkbox_custom, 'Value')
                    handles.scales(1) = get(handles.slider_1, 'Value');
                    handles.scales(2) = get(handles.slider_2, 'Value');
                    handles.multFactor4Thresh = get(handles.slider_3, 'Value');
                    optional_arguments ={handles.scales,handles.multFactor4Thresh};
                else                    
                    optional_arguments={};
                end
                plusTipCometDetector(handles.projList(i),...
                        handles.timeRangeDetect,handles.bitDepth,...
                        handles.savePlots,optional_arguments{:});
            case 2
                psfSigma = str2double(get(handles.edit_psfSigma,'String'));
                alpha = str2double(get(handles.edit_alpha,'String'));
                displayFirstImage = get(handles.checkbox_displayFirstImage,'Value');
                plusTipAnisoGaussianCometDetector(handles.projList(i),psfSigma,...
                    handles.timeRangeDetect,handles.bitDepth,handles.savePlots,...
                    'alpha',alpha,'displayFirstImage',displayFirstImage);
            otherwise
                error('Unrecognized detection method')
        end
        toc
    end

    % tracking
    if handles.doTrack==1
        tic
        disp(['Tracking project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
        plusTipCometTracker(handles.projList(i),handles.timeWindow,...
            handles.minTrackLen,handles.minRadius,handles.maxRadius,...
            handles.maxFAngle,handles.maxBAngle,handles.maxShrinkFactor,...
            handles.fluctRad,handles.timeRangeTrack);
        toc
    end
    
    % post-processing
    if handles.doMeta==1
        tic
        disp(['Post-processing project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
        [projData]=plusTipPostTracking(handles.projList(i),...
            handles.secPerFrame,handles.pixSizeNm,handles.timeRangePost,handles.doHist);
        toc
    end
%     catch
%         disp(['Problem with ' handles.projList(i).anDir])
%     end
end
disp('Finished!')

% --- Executes on button press in histCheck.
function histCheck_Callback(hObject, eventdata, handles)
% hObject    handle to histCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of histCheck
handles.doHist=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in setupRoiPush.
function setupRoiPush_Callback(hObject, eventdata, handles)
overwriteROIs=0;
doCrop=0;
selectRoiTag = get(get(handles.uipanel_roi,'SelectedObject'),'Tag');
selectRoi = str2double(selectRoiTag(end));
setupRoiDirectories(selectRoi,overwriteROIs,doCrop);

% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plusTipGuiSwitch(hObject,eventdata,handles,'resetButton');   


% --- Executes on button press in getProjListFile_check.
function getProjListFile_check_Callback(hObject, eventdata, handles)
% hObject    handle to getProjListFile_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getProjListFile_check
handles.loadProjList=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in pushbutton_preview.
function pushbutton_preview_Callback(hObject, eventdata, handles)

userData= get(handles.figure1,'UserData');
if ~isempty(userData.previewGUI) && ishandle(userData.previewGUI)
   delete(userData.previewGUI) 
end
userData.previewGUI = detectionPreviewGUI ('mainFig', handles.figure1);
set(handles.figure1,'UserData',userData)


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
if ~isempty(userData.previewGUI) && ishandle(userData.previewGUI)
   delete(userData.previewGUI) 
end

function editWatershed(hObject, eventdata, handles)

% Retrieve tag and style and read the new and old value
props = get(hObject,{'Tag','Style'});
id = str2double(props{1}(end));
if strcmp(props{2},'slider')
    value = get(hObject,'Value');
    oldvalue = str2double(get(handles.(['edit_detect_' num2str(id)]),'String'));
else
    value = str2double(get(hObject,'String'));
    oldvalue = get(handles.(['slider_' num2str(id)]),'Value');
end

% DO a series of test and reset value to old value in case of failure
nanTest = isnan(value);
minmaxTest = value<get(handles.(['slider_' num2str(id)]),'Min') ||...
        value>get(handles.(['slider_' num2str(id)]),'Max');
diffTest1 = (id==2 && value<=get(handles.slider_1, 'Value')); 
diffTest2 = (id==1 && value>=get(handles.slider_2, 'Value')); 
if nanTest || minmaxTest || diffTest1 || diffTest2,  value=oldvalue; end

% Update the slider and the edit_box
set(handles.(['slider_' num2str(id)]),'Value',value);
set(handles.(['edit_detect_' num2str(id)]),'String',value);

function fluctRadEdit_Callback(hObject, eventdata, handles)
handles.fluctRad=str2double(get(hObject,'String'));
guidata(hObject, handles);
