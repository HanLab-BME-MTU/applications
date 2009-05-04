function varargout = plusTipTrackerGUI(varargin)
% PLUSTIPTRACKERGUI M-file for plusTipTrackerGUI.fig
%      PLUSTIPTRACKERGUI, by itself, creates a new PLUSTIPTRACKERGUI or raises the existing
%      singleton*.
%
%      H = PLUSTIPTRACKERGUI returns the handle to a new PLUSTIPTRACKERGUI or the handle to
%      the existing singleton*.
%
%      PLUSTIPTRACKERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPTRACKERGUI.M with the given input arguments.
%
%      PLUSTIPTRACKERGUI('Property','Value',...) creates a new PLUSTIPTRACKERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipTrackerGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipTrackerGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipTrackerGUI

% Last Modified by GUIDE v2.5 04-May-2009 10:38:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipTrackerGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipTrackerGUI_OutputFcn, ...
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


% --- Executes just before plusTipTrackerGUI is made visible.
function plusTipTrackerGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipTrackerGUI (see VARARGIN)

% Choose default command line output for plusTipTrackerGUI
handles.output = hObject;
handles.selectRoi = 1;
handles.getStr = 1;

handles.doDetect=0;
handles.doTrack=0; 
handles.doMeta=0; 

% DETECTION parameters
handles.timeRange = [1 inf];
handles.bitDepth  = 14; % change according to your camera (should be 12, 14, or 16)
handles.savePlots = 1;  % 1 to save overlay plots of detection results in subfolder; 0 if not (may run faster)

% TRACKING parameters
handles.timeWindow=[];
handles.minTrackLen=3;
handles.minRadius=[];
handles.maxRadius=[];
handles.maxFAngle=45;
handles.maxBDist=1.5;

% META parameters
handles.secPerFrame=[];
handles.pixSizeNm=[];
handles.doHist=1;

% Update handles structure
guidata(hObject,handles);

% UIWAIT makes plusTipTrackerGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipTrackerGUI_OutputFcn(hObject, eventdata, handles) 
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
if handles.getStr==1
    handles.strList=inputGUI;
    [projList,projPathParsed]=getProj(handles.strList);
else
    handles.strList=[];
    [projList,projPathParsed]=getProj;
end
if ~isempty(projList)
    a=struct2cell(projList);
    a=a(2,:)';
    a=sort(a);
    b=cellfun(@isempty, strfind(a,'sub'));
    [selection,selectionList]=listSelectGUI(a(b),[],'move');
    handles.projList=projList(selection,1);
else
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

% --- Executes on button press in doDetectCheck.
function doDetectCheck_Callback(hObject, eventdata, handles)
% hObject    handle to doDetectCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doDetectCheck
handles.doDetect=get(hObject,'Value');
if handles.doDetect==1
    set(handles.bitDepthEdit,'Enable','on')
    set(handles.startFrame,'Enable','on')
    set(handles.endFrame,'Enable','on')
    set(handles.savePlotCheck,'Enable','on')
else
    set(handles.bitDepthEdit,'Enable','off')
    set(handles.startFrame,'Enable','off')
    set(handles.endFrame,'Enable','off')
    set(handles.savePlotCheck,'Enable','off')
end
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



function maxAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of maxAngleEdit as a double
handles.maxFAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxAngleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxBDistEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxBDistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxBDistEdit as text
%        str2double(get(hObject,'String')) returns contents of maxBDistEdit as a double
handles.maxBDist=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxBDistEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxBDistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
    set(handles.maxAngleEdit,'Enable','on')
    set(handles.maxBDistEdit,'Enable','on')
else
    set(handles.timeWindowEdit,'Enable','off')
    set(handles.minTrackLengthEdit,'Enable','off')
    set(handles.minRadiusEdit,'Enable','off')
    set(handles.maxRadiusEdit,'Enable','off')
    set(handles.maxAngleEdit,'Enable','off')
    set(handles.maxBDistEdit,'Enable','off')
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
    set(handles.frameRateEdit,'Enable','on')
    set(handles.pixSizeEdit,'Enable','on')
    set(handles.histCheck,'Enable','on')
else
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
numProj=size(handles.projList,1);
for i=1:numProj
%    try
        % detection
        if handles.doDetect==1
            disp(['Detecting project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            eb1SpotDetector(handles.projList(i),handles.timeRange,handles.bitDepth,handles.savePlots);
        end

        % tracking
        if handles.doTrack==1
            disp(['Tracking project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            plusTipTracker(handles.projList(i),handles.timeWindow,handles.minTrackLen,...
                handles.minRadius,handles.maxRadius,handles.maxFAngle,handles.maxBDist);
        end
        if handles.doMeta==1
            disp(['Post-processing project ' num2str(i) filesep num2str(numProj) ': ' handles.projList(i).anDir])
            [projData]=metaEB3analysis(handles.projList(i),handles.secPerFrame,handles.pixSizeNm);
            if handles.doHist==1
                popHist(projData);
            end
        end
%     catch
%         disp(['Problem with ' handles.projList(i).anDir])
%     end
end

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

% --- Executes on mouse press over axes background.
function helpPic_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to helpPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open plusTipAnalysis_README.txt


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

