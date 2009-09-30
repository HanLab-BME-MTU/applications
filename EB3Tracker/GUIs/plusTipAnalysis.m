function varargout = plusTipAnalysis(varargin)
% PLUSTIPANALYSIS M-file for plusTipAnalysis.fig
%      PLUSTIPANALYSIS, by itself, creates a new PLUSTIPANALYSIS or raises the existing
%      singleton*.
%
%      H = PLUSTIPANALYSIS returns the handle to a new PLUSTIPANALYSIS or the handle to
%      the existing singleton*.
%
%      PLUSTIPANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPANALYSIS.M with the given input arguments.
%
%      PLUSTIPANALYSIS('Property','Value',...) creates a new PLUSTIPANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipAnalysis
% 
%
% adding space to test SVN
%
% Last Modified by GUIDE v2.5 29-Aug-2009 17:58:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipAnalysis_OutputFcn, ...
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


% --- Executes just before plusTipAnalysis is made visible.
function plusTipAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipAnalysis (see VARARGIN)

% Choose default command line output for plusTipAnalysis
handles.output = hObject;
handles.selectRoi = 1;
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







% Update handles structure
guidata(hObject,handles);

% UIWAIT makes plusTipAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipAnalysis_OutputFcn(hObject, eventdata, handles) 
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
    a=struct2cell(handles.projList);
    if isempty(strfind(a{1,1},'roi_'))
        a=a(2,:)';
    else
        a=a(1,:)';
    end
    a=sort(a);
    b=cellfun(@isempty, strfind(a,'sub'));
    a=a(b);
    
    % don't look for tracking results as in plusTipTrackViz - need to
    % collect untracked projects for tracking!
    
    % allow multiple projects to be selected
    [selection,selectionList]=listSelectGUI(a,[],'move',1);

    % if a project was selected, save projData info and get data
    if ~isempty(selection)
        handles.projList=handles.projList(selection,1);
    else
        handles.projList=[];
    end
else
    msgbox('No projects selected.')
    handles.projList=[];
end
guidata(hObject, handles);


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

% --- Executes during object creation, after setting all properties.
function startFrameDetect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrameDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endFrameDetect_Callback(hObject, eventdata, handles)
% hObject    handle to endFrameDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrameDetect as text
%        str2double(get(hObject,'String')) returns contents of endFrameDetect as
%        a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameDetect');  
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function endFrameDetect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrameDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startFramePost_Callback(hObject, eventdata, handles)
% hObject    handle to startFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFramePost as text
%        str2double(get(hObject,'String')) returns contents of startFramePost as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFramePost');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function startFramePost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endFramePost_Callback(hObject, eventdata, handles)
% hObject    handle to endFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFramePost as text
%        str2double(get(hObject,'String')) returns contents of endFramePost as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFramePost');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function endFramePost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFramePost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startFrameTrack_Callback(hObject, eventdata, handles)
% hObject    handle to startFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrameTrack as text
%        str2double(get(hObject,'String')) returns contents of startFrameTrack as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'startFrameTrack');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function startFrameTrack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endFrameTrack_Callback(hObject, eventdata, handles)
% hObject    handle to endFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrameTrack as text
%        str2double(get(hObject,'String')) returns contents of endFrameTrack as a double
handles=plusTipGuiSwitch(hObject,eventdata,handles,'endFrameTrack');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function endFrameTrack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrameTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bitDepthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bitDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bitDepthEdit as text
%        str2double(get(hObject,'String')) returns contents of bitDepthEdit as a double
handles.bitDepth=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function bitDepthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bitDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

% --- Executes during object creation, after setting all properties.
function frameRateEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameRateEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pixSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of pixSizeEdit as a double
handles.pixSizeNm=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pixSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timeWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to timeWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of timeWindowEdit as a double
handles.timeWindow=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function timeWindowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minTrackLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minTrackLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minTrackLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of minTrackLengthEdit as a double
handles.minTrackLen=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minTrackLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minTrackLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of minRadiusEdit as a double
handles.minRadius=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minRadiusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of maxRadiusEdit as a double
handles.maxRadius=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxRadiusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxFAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxFAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxFAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of maxFAngleEdit as a double
handles.maxFAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxFAngleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxFAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxBAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxBAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxBAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of maxBAngleEdit as a double
handles.maxBAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxBAngleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxBAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxShrinkFactorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxShrinkFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxShrinkFactorEdit as text
%        str2double(get(hObject,'String')) returns contents of maxShrinkFactorEdit as a double
handles.maxShrinkFactor=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxShrinkFactorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxShrinkFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doDetectCheck.
function doDetectCheck_Callback(hObject, eventdata, handles)
% hObject    handle to doDetectCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doDetectCheck
handles.doDetect=get(hObject,'Value');
if handles.doDetect==1
    set(handles.bitDepthEdit,'Enable','on')
    set(handles.startFrameDetect,'Enable','on')
    set(handles.endFrameDetect,'Enable','on')
    set(handles.savePlotCheck,'Enable','on')
else
    set(handles.bitDepthEdit,'Enable','off')
    set(handles.startFrameDetect,'Enable','off')
    set(handles.endFrameDetect,'Enable','off')
    set(handles.savePlotCheck,'Enable','off')
end
guidata(hObject, handles);


% --- Executes on button press in trackingCheck.
function trackingCheck_Callback(hObject, eventdata, handles)
% hObject    handle to trackingCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trackingCheck
handles.doTrack=get(hObject,'Value');
if handles.doTrack==1
    set(handles.timeWindowEdit,'Enable','on')
    set(handles.minTrackLengthEdit,'Enable','on')
    set(handles.minRadiusEdit,'Enable','on')
    set(handles.maxRadiusEdit,'Enable','on')
    set(handles.maxFAngleEdit,'Enable','on')
    set(handles.maxBAngleEdit,'Enable','on')
    set(handles.maxShrinkFactorEdit,'Enable','on')
    set(handles.fluctRadEdit,'Enable','on')
    set(handles.startFrameTrack,'Enable','on')
    set(handles.endFrameTrack,'Enable','on')
    
else
    set(handles.timeWindowEdit,'Enable','off')
    set(handles.minTrackLengthEdit,'Enable','off')
    set(handles.minRadiusEdit,'Enable','off')
    set(handles.maxRadiusEdit,'Enable','off')
    set(handles.maxFAngleEdit,'Enable','off')
    set(handles.maxBAngleEdit,'Enable','off')
    set(handles.maxShrinkFactorEdit,'Enable','off')
    set(handles.fluctRadEdit,'Enable','off')
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


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


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
            disp(['Detecting project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            plusTipCometDetector(handles.projList(i),...
                handles.timeRangeDetect,handles.bitDepth,handles.savePlots);
        end

        % tracking
        if handles.doTrack==1
            disp(['Tracking project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            plusTipTracker(handles.projList(i),handles.timeWindow,...
                handles.minTrackLen,handles.minRadius,handles.maxRadius,...
                handles.maxFAngle,handles.maxBAngle,handles.maxShrinkFactor,...
                handles.fluctRad,handles.timeRangeTrack);
        end
        
        % post-processing
        if handles.doMeta==1
            disp(['Post-processing project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            [projData]=plusTipPostTracking(handles.projList(i),...
                handles.secPerFrame,handles.pixSizeNm,handles.timeRangePost,handles.doHist);
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
% hObject    handle to setupRoiPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setupRoiDirectories(handles.selectRoi,1);


% --- Executes on button press in selectRoiCheck.
function selectRoiCheck_Callback(hObject, eventdata, handles)
% hObject    handle to selectRoiCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectRoiCheck
handles.selectRoi=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plusTipGuiSwitch(hObject,eventdata,handles,'resetButton');   



function fluctRadEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fluctRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fluctRadEdit as text
%        str2double(get(hObject,'String')) returns contents of fluctRadEdit as a double
handles.fluctRad=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fluctRadEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fluctRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function fluctRad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fluctRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in getProjListFile_check.
function getProjListFile_check_Callback(hObject, eventdata, handles)
% hObject    handle to getProjListFile_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getProjListFile_check
handles.loadProjList=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in getHelpPush.
function getHelpPush_Callback(hObject, eventdata, handles)
% hObject    handle to getHelpPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open plusTipAnalysis_README.txt

% --- Executes during object creation, after setting all properties.
function getHelpPush_CreateFcn(hObject, eventdata, handles)
% hObject    handle to getHelpPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'CData',imread('help_icon.png'));



