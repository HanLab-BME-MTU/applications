function varargout = labelgui2(varargin)
% LABELGUI2 M-file for labelgui2.fig
%      LABELGUI2, by itself, creates a new LABELGUI2 or raises the existing
%      singleton*.
%
%      H = LABELGUI2 returns the handle to a new LABELGUI2 or the handle to
%      the existing singleton*.
%
%      LABELGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABELGUI2.M with the given input arguments.
%
%      LABELGUI2('Property','Value',...) creates a new LABELGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before labelgui2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to labelgui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help labelgui2

% Last Modified by GUIDE v2.5 02-Mar-2007 14:17:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labelgui2_OpeningFcn, ...
                   'gui_OutputFcn',  @labelgui2_OutputFcn, ...
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


% --- Executes just before labelgui2 is made visible.
function labelgui2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labelgui2 (see VARARGIN)

% Choose default command line output for labelgui2
handles.output = hObject;

% initialize fields
handles.movieWindowH = [];
handles.currentDir = [];
% if you're adding fields here, please don't forget to update
% LG_launchMovieWindow
handles.positions.LG_movieWindow = [];
handles.positions.LG_movieDataFigure = [];
handles.positions.LG_reAssignGUI = [];
handles.positions.LG_testRatiosFigure = [];
handles.positions.LG_intensityFigure = [];
handles.positions.LG_distanceFigure = [];
handles.positions.LG_displacementFigure = [];

% set flag if called from runCtBatch (or similar)
if ~isempty(varargin) && ~isempty(varargin{1})
    handles.launchedFromOutside = varargin{1};
else
    handles.launchedFromOutside = 0;
end

% set position
screenSize = get(0,'MonitorPositions');
% get labelgui-position in pixels
units = get(hObject,'Units');
set(hObject,'Units','Pixels');
navigatorPosition = get(hObject,'Position');
% x = 5, y = screenHeight - 5. Be careful with chars/pix
navigatorPosition(1) = 5;
navigatorPosition(2) = screenSize(4) - navigatorPosition(4);
set(hObject,'Position',navigatorPosition);
set(hObject,'Units',units);
% only now adjust screenHeight
navigatorPosition = get(hObject,'Position');
navigatorPosition(2) = navigatorPosition(2) - 3.5;
set(hObject,'Position',navigatorPosition);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = labelgui2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in LG_navi_zSliceTitle_txt.
function LG_navi_zSliceTitle_txt_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_zSliceTitle_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get value. If checked, turn on slider and edt. Otherwise, turn off
isChecked = get(hObject,'Value');
if isChecked
    set(handles.LG_navi_zSliceSlider_sli,'Enable','on')
    set(handles.LG_navi_zSliceNumber_edt,'Visible','on')
    % plot slice
    LG_gotoFrame([], round(get(handles.LG_navi_zSliceSlider_sli,'Value')))
else
set(handles.LG_navi_zSliceSlider_sli,'Enable','off')
    set(handles.LG_navi_zSliceNumber_edt,'Visible','off')
    % plot entire frame
    LG_gotoFrame(round(get(handles.LG_navi_timepointSlider_sli,'Value')))
end




% --- Executes on slider movement.
function LG_navi_zSliceSlider_sli_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_zSliceSlider_sli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get time
currentSlice = round(get(hObject,'Value'));

% goto frame
LG_gotoFrame([],currentSlice);

% --- Executes during object creation, after setting all properties.
function LG_navi_zSliceSlider_sli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_zSliceSlider_sli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function LG_navi_zSliceNumber_edt_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_zSliceNumber_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get string, check for nan and out-of-range
timeString = get(hObject,'String');
currentSlice = str2double(timeString);
if isnan(currentSlice)
    h = errordlg(sprintf('%s is not a valid slice number!',timeString));
    uiwait(h)
    return
end

if currentSlice < 1
    currentSlice = 1;
elseif currentSlice > get(handles.LG_navi_zSliceSlider_sli,'Max')
    currentSlice = get(handles.LG_navi_zSliceSlider_sli,'Max');
end

% gotoFrame
LG_gotoFrame([],currentSlice);
    
% --- Executes during object creation, after setting all properties.
function LG_navi_zSliceNumber_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_zSliceNumber_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes during object creation, after setting all properties.
function LG_navi_flagName_pd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_flagName_pd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function LG_navi_timepointSlider_sli_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_timepointSlider_sli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get time
currentTime = round(get(hObject,'Value'));

% goto frame
LG_gotoFrame(currentTime);

% --- Executes during object creation, after setting all properties.
function LG_navi_timepointSlider_sli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_timepointSlider_sli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function LG_navi_timepointNumber_edt_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_timepointNumber_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get string, check for nan and out-of-range
timeString = get(hObject,'String');
currentTime = str2double(timeString);
if isnan(currentTime)
    h = errordlg(sprintf('%s is not a valid timepoint!',timeString));
    uiwait(h)
    return
end

if currentTime < 1
    currentTime = 1;
elseif currentTime > get(handles.LG_navi_timepointSlider_sli,'Max')
    currentTime = get(handles.LG_navi_timepointSlider_sli,'Max');
end

% gotoFrame
LG_gotoFrame(currentTime);

% --- Executes during object creation, after setting all properties.
function LG_navi_timepointNumber_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_timepointNumber_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in LG_navi_menuLoadMovie_navi_flagLast_pb.
function LG_navi_flagLast_pb_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_flagLast_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LG_gotoFlag(-1);


% --- Executes on button press in LG_navi_menuLoadMovie_navi_flagNext_pb.
function LG_navi_flagNext_pb_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuLoadMovie_navi_flagNext_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LG_gotoFlag(+1);

% --------------------------------------------------------------------










% --------------------------------------------------------------------
function LG_navi_menuShowIntensities_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuShowIntensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function LG_navi_menuShowDistances_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuShowDistances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LG_menuView_Callback(hObject, eventdata, handles)
% hObject    handle to LG_menuView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LG_navi_menuShowDisplacements_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuShowDisplacements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function LG_navi_menuShowMovieData_Callback(hObject, eventdata, handles)
% hObject    handle to LG_navi_menuShowMovieData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


