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

% Last Modified by GUIDE v2.5 11-Apr-2009 13:45:22

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




% --- Executes on button press in chooseProjData.
function chooseProjData_Callback(hObject, eventdata, handles)
% hObject    handle to chooseProjData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load projData from meta
[FileName,PathName] = uigetfile('*.mat','Select projData from project meta folder');
if ~isequal(FileName,0)
    cd(PathName)
    cd ..
    handles.dataDir=pwd;


    p=load([PathName filesep FileName]);
    handles.projData=p.projData;

    % load tracksFinal from track
    anDir=formatPath(handles.projData.anDir);
    trackDir=[anDir filesep 'track'];
    [listOfFiles]=searchFiles('.mat',[],trackDir,0);
    if ~isempty(listOfFiles)
        load([listOfFiles{1,2} filesep listOfFiles{1,1}])
        if ~exist('tracksFinal','var')
            error('--trackMovie: tracksFinal missing...');
        end
    else
        error('--ebTrackViz: tracksFinal missing...');
    end
    handles.tracksFinal=tracksFinal;
end
guidata(hObject, handles);




% --- Executes on button press in selectSavedRoiPushbutton.
function selectSavedRoiPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectSavedRoiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat or roiMask.tif');
if ~isequal(FileName,0)
    if ~isempty(strfind(FileName,'tif'))
        handles.roi=imread([PathName FileName]);
    else
        p=load([PathName FileName]);
        handles.roi=p.roiYX;
    end
end
guidata(hObject, handles);





function startFrame_Callback(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrame as text
%        str2double(get(hObject,'String')) returns contents of startFrame as a double
sFVal=get(hObject,'String');
if isequal(lower(sFVal),'min')
    handles.timeRange(1)=1;
else
    handles.timeRange(1)=str2double(sFVal);
end
    
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
%        str2double(get(hObject,'String')) returns contents of endFrame as a double
eFVal=get(hObject,'String');
if isequal(lower(eFVal),'max')
    handles.timeRange(2)=handles.projData.numFrames;
else
    handles.timeRange(2)=str2double(eFVal);
end
    
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
val=get(hObject,'Value');
if val==0
    handles.ask4select=0;
else
    handles.ask4select=1;
end
guidata(hObject, handles);




% --- Executes on button press in showTracksCheck.
function showTracksCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showTracksCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTracksCheck
val=get(hObject,'Value');
if val==0
    handles.showTracks=0;
else
    handles.showTracks=1;
end
guidata(hObject, handles);




% --- Executes on button press in showDetectCheck.
function showDetectCheck_Callback(hObject, eventdata, handles)
% hObject    handle to showDetectCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showDetectCheck
val=get(hObject,'Value');
if val==0
    handles.showDetect=0;
else
    handles.showDetect=1;
end
guidata(hObject, handles);




function speedLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to speedLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of speedLimitEdit as a double
velLimVal=get(hObject,'String');
if isequal(lower(velLimVal),'max')
    handles.velLimit=inf;
else
    handles.velLimit=str2double(velLimVal);
end
    
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
%        str2double(get(hObject,'String')) returns contents of indivTrackNumbersEdit as a double
userInput=get(hObject,'String');
handles.indivTrack=str2num(userInput)';

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

selectedTracksDisplay_Callback(hObject, eventdata, handles)
if ~isempty(handles.projData)
    cd(formatPath(handles.projData.imDir))
end
[handles.selectedTracks] = plotTracks2D_EB3(handles.tracksFinal,...
    handles.timeRange,handles.img,handles.ask4select,...
    handles.plotCurrentOnly,handles.roi,handles.movieInfo);

cd(handles.homeDir)

if ~isempty(handles.selectedTracks)
    
    temp=vertcat(handles.selectedTracks{:});
    handles.selectedTracks=unique(temp(:,1));
    [l w]=size(handles.selectedTracks);

    if l>w
        handles.selectedTracks=handles.selectedTracks';
    end
    
end
selectedTracksDisplay_Callback(hObject, eventdata, handles)
guidata(hObject, handles);



% --- Executes on button press in aviCheckTrackMov.
function aviCheckTrackMov_Callback(hObject, eventdata, handles)
% hObject    handle to aviCheckTrackMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckTrackMov
val=get(hObject,'Value');
if val==0
    handles.doAvi=0;
else
    handles.doAvi=1;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function selectedTracksDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function selectedTracksDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to selectedTracksDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectedTracksDisplay as text
%        str2double(get(hObject,'String')) returns contents of selectedTracksDisplay as a double
hObject=handles.selectedTracksDisplay;
if isempty(handles.selectedTracks)
    set(hObject,'Visible','Off');
else
    set(hObject,'Visible','On');
end
set(hObject,'String',num2str(handles.selectedTracks));
guidata(hObject, handles);



% --- Executes on button press in speedMovieButton.
function speedMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to speedMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
featVelMovie(handles.projData,handles.timeRange,handles.velLimit,handles.roi,handles.doAvi);
guidata(hObject, handles);


% --- Executes on button press in trackMovieButton.
function trackMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to trackMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trackMovie(handles.projData,handles.indivTrack,handles.timeRange,handles.roi,handles.magCoef,handles.showTracks,handles.showDetect,handles.doAvi);
guidata(hObject, handles);


% --- Executes on button press in aviCheckSpeedMov.
function aviCheckSpeedMov_Callback(hObject, eventdata, handles)
% hObject    handle to aviCheckSpeedMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aviCheckSpeedMov
val=get(hObject,'Value');
if val==0
    handles.doAvi=0;
else
    handles.doAvi=1;
end
guidata(hObject, handles);


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closeGUI = handles.figure1; %handles.figure1 is the GUI figure
 
guiPosition = get(handles.figure1,'Position'); %get the position of the GUI
guiName = get(handles.figure1,'Name'); %get the name of the GUI
eval(guiName) %call the GUI again
 
close(closeGUI); %close the old GUI
set(gcf,'Position',guiPosition); %set the position for the new GUI


% --- Executes during object creation, after setting all properties.
function helpBoxAxes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to helpBoxAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate helpBoxAxes1
img=imread('qIcon.jpg');
imagesc(img,'parent',hObject);
axis image
axis off
imHandle=get(hObject,'Children');
%info=get(hObject); info2=get(blah);
%set(blah,'HitTest','off')
set(imHandle,'ButtonDownFcn',@helpBoxAxes1_ButtonDownFcn)


%function mouseClickCallback(hObject, eventdata, handles)
%disp('blah')


% --- Executes on mouse press over axes background.
function helpBoxAxes1_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to helpBoxAxes1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % handles.helpNote=1;
% % guidata(hObject, handles);
% % helpNotes_Callback(hObject, eventdata, handles)
% handles.helpNotes
open plusTipTrackViz_README.txt