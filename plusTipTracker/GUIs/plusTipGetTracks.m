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
% Last Modified by GUIDE v2.5 17-Oct-2011 10:19:39

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
handles.projList = [];
handles.nProjLists = 0;

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

% Allow the user to save the project list if not from a single projList
isSingleProjList = handles.nProjLists == 1 && ~get(handles.getQueryStr_Check,'Value');
if ~isempty(handles.projList) && ~isSingleProjList
    if ~isempty(handles.strList)
        defaultListName = ['projList', sprintf('_%s',handles.strList{:}) '.mat'];
    else
        defaultListName='projList.mat';
    end
    [file,path] = uiputfile('projList.mat','Find a location to save your project list',...
        [pwd filesep defaultListName]);
    if isequal(file,0), return; end
    save([path file],'-struct','handles','projList');
end


% --- Executes on button press in detectionCheck.
function updateDetection(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of detectionCheck
doDetect=get(handles.detectionCheck,'Value');
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
    
if doDetect
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

trackingHandles=get(handles.trackingPanel,'Children');
if get(hObject,'Value')
    set(trackingHandles,'Enable','on')
    set(handles.startFrameTrack,'Enable','on')
    set(handles.endFrameTrack,'Enable','on')
else
    set(trackingHandles,'Enable','off')
    set(handles.startFrameTrack,'Enable','off')
    set(handles.endFrameTrack,'Enable','off')
end

% --- Executes on button press in metaCheck.
function metaCheck_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
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
doDetect=get(handles.detectionCheck,'Value');
doTrack=get(handles.trackingCheck,'Value');
doMeta=get(handles.metaCheck,'Value');

% Read detection parameters
if doDetect
    detectionParams=[];
    
    % Read detection time range
    sFVal=get(handles.startFrameDetect,'String');
    if strcmpi(sFVal,'min')
        detectionParams.timeRangeDetect(1)=1;
    else
        detectionParams.timeRangeDetect(1)=str2double(sFVal);
    end  
    eFVal=get(handles.endFrameDetect,'String');
    if strcmpi(eFVal,'max')
        detectionParams.timeRangeDetect(2)=Inf;
    else
        detectionParams.timeRangeDetect(2)=str2double(eFVal);
    end
    if ~all(isposint(detectionParams.timeRangeDetect)),
        errordlg('Invalid detection parameter range. Please check again','Missing Input');
        return
    end
    
    % Read detection common parameters
    detectionParams.bitDepth = str2double(get(handles.bitDepth,'String'));
    detectionParams.savePlots = get(handles.savePlotCheck,'value');
    detectionParams.addArgs={};
    
    % Read detection method specific parameters      
    detectionMethod = get(handles.popupmenu_detectionMethod,'Value');
    if get(handles.checkbox_custom, 'Value') && detectionMethod ==1;
        scales(1) = get(handles.slider_1, 'Value');
        scales(2) = get(handles.slider_2, 'Value');
        multFactor4Thresh = get(handles.slider_3, 'Value');
        detectionParams.addArgs ={scales,multFactor4Thresh};
    end
    if detectionMethod==2
        detectionParams.psfSigma = str2double(get(handles.edit_psfSigma,'String'));
    end
    if get(handles.checkbox_custom, 'Value') && detectionMethod ==2;
        alpha = str2double(get(handles.edit_alpha,'String'));
        displayFirstImage = get(handles.checkbox_displayFirstImage,'Value');
        detectionParams.addArgs ={'alpha',alpha,'displayFirstImage',displayFirstImage};
    end
end

% Read tracking parameters
if doTrack
    trackParams=[];
    
    % Read tracking time range
    sFVal=get(handles.startFrameTrack,'String');
    if strcmpi(sFVal,'min')
        trackParams.timeRangeTrack(1)=1;
    else
        trackParams.timeRangeTrack(1)=str2double(sFVal);
    end  
    eFVal=get(handles.endFrameTrack,'String');
    if strcmpi(eFVal,'max')
        trackParams.timeRangeTrack(2)=Inf;
    else
        trackParams.timeRangeTrack(2)=str2double(eFVal);
    end
    if ~all(isposint(trackParams.timeRangeDetect)),
        errordlg('Invalid tracking parameter range. Please check again','Missing Input');
        return
    end
    
    % Read tracking numeric parameters
    trackParamsNames = {'timeWindow','minTrackLength','minRadius','maxRadius',...
        'maxFAngle','maxBAngle','maxShrinkFactor','fluctRad'};
    for j=1:numel(trackParamsNames);
        value = str2double(get(handles.(trackParamsNames{j}),'String'));
        if isnan(value) || value < 0
            errordlg('Invalid tracking parameter value. Please check again','Missing Input');
            return
        end
        trackParams.(trackParamsNames{j})=value;
    end        
end

% Read meta parameters
if doMeta
    metaParams=[];
    % Read post-processing time range
    sFVal=get(handles.startFramePost,'String');
    if strcmpi(sFVal,'min')
        metaParams.timeRangePost(1)=1;
    else
        metaParams.timeRangePost(1)=str2double(sFVal);
    end  
    eFVal=get(handles.endFramePost,'String');
    if strcmpi(eFVal,'max')
        metaParams.timeRangePost(2)=Inf;
    else
        metaParams.timeRangePost(2)=str2double(eFVal);
    end
    if ~all(isposint(metaParams.timeRangeDetect)),
        errordlg('Invalid post-processing parameter range. Please check again','Missing Input');
        return
    end
    
    % Read post-processing parameters
    metaParams.secPerFrame = str2double(get(handles.frameRateEdit,'String'));
    metaParams.pixSizeNm = str2double(get(handles.pixSizeEdit,'String'));
    metaParams.doHist=get(handles.histCheck,'Value');
    if isnan(metaParams.secPerFrame) || isnan(metaParams.pixSizeNm)
        errordlg('Invalid post-processing parameter value. Please check again','Missing Input');
        return
    end       
end


for i=1:numProj
    %try
    % detection
    if doDetect
        tic
        disp(['Detecting project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir]);
        detectionMethod = get(handles.popupmenu_detectionMethod,'Value');
        switch detectionMethod
            case 1
                plusTipCometDetector(handles.projList(i),...
                    detectionParams.timeRangeDetect,detectionParams.bitDepth,...
                    detectionParams.savePlots,detectionParams.addArgs{:});
            case 2
                plusTipAnisoGaussianCometDetector(handles.projList(i),...
                    detectionParams.psfSigma,detectionParams.timeRangeDetect,...
                    detectionParams.bitDepth,detectionParams.savePlots,...
                    detectionParams.addArgs{:});
            otherwise
                error('Unrecognized detection method')
        end
        toc
    end

    % tracking
    if doTrack==1
        tic
        disp(['Tracking project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
        plusTipCometTracker(handles.projList(i),trackParams.timeWindow,...
            trackParams.minTrackLength,trackParams.minRadius,trackParams.maxRadius,...
            trackParams.maxFAngle,trackParams.maxBAngle,trackParams.maxShrinkFactor,...
            trackParams.fluctRad,trackParams.timeRangeTrack);
        toc
    end
    
    % post-processing
    if doMeta
        tic
        disp(['Post-processing project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
        plusTipPostTracking(handles.projList(i),metaParams.secPerFrame,...
            metaParams.pixSizeNm,metaParams.timeRangePost,metaParams.doHist);
        toc
    end
end
disp('Finished!')

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
