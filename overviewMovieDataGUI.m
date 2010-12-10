function varargout = overviewMovieDataGUI(varargin)
% OVERVIEWMOVIEDATAGUI M-file for overviewMovieDataGUI.fig
%      OVERVIEWMOVIEDATAGUI, by itself, creates a new OVERVIEWMOVIEDATAGUI or raises the existing
%      singleton*.
%
%      H = OVERVIEWMOVIEDATAGUI returns the handle to a new OVERVIEWMOVIEDATAGUI or the handle to
%      the existing singleton*.
%
%      OVERVIEWMOVIEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OVERVIEWMOVIEDATAGUI.M with the given input arguments.
%
%      OVERVIEWMOVIEDATAGUI('Property','Value',...) creates a new OVERVIEWMOVIEDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before overviewMovieDataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to overviewMovieDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help overviewMovieDataGUI

% Last Modified by GUIDE v2.5 10-Jun-2010 12:49:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @overviewMovieDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @overviewMovieDataGUI_OutputFcn, ...
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


% --- Executes just before overviewMovieDataGUI is made visible.
function overviewMovieDataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Useful tools
%       userData.mainFig - handle of main figure
%       

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

% Get main figure handle and MovieData obj
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};

userData_main = get(userData.mainFig, 'UserData');
userData.MD = userData_main.MD;
MD = userData_main.MD;

% Notify main figure that the current process setting panel is open
setappdata(userData.mainFig,'overviewFlag',1);

% -------------------- Write MovieData Information ---------------------

set(handles.text_1, 'String', [MD.movieDataPath_ MD.movieDataFileName_]);

hr = num2str(mod(MD.createTime_(4),12));
min = num2str(MD.createTime_(5));
if MD.createTime_(4) >=12
    t = 'PM';
else
    t = 'AM';
end
if MD.createTime_(5) == 0
    min = '00';
end
if MD.createTime_(4) == 0
    hr = '00';
end
set(handles.text_2, 'String', ...
    [num2str(MD.createTime_(2)),'/',num2str(MD.createTime_(3)),'/',...
    num2str(MD.createTime_(1)),'  ',hr,':',min,' ',t]);

set(handles.text_3, 'String', num2str(MD.pixelSize_));

set(handles.text_4, 'String', num2str(MD.timeInterval_));

set(handles.text_5, 'String', ['Width:  ',num2str(MD.imSize_(1)), ...
    '      Height:  ',num2str(MD.imSize_(2))]);

set(handles.text_6, 'String', num2str(MD.nFrames_));

set(handles.text_7, 'String', num2str(length(MD.packages_)));

set(handles.text_8, 'String', num2str(length(MD.channelPath_)));

set(handles.listbox_channel, 'String', MD.channelPath_);

set(handles.edit_notes, 'String', MD.notes_);

% ----------------------------------------------------------------------

% Update handles structure
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% UIWAIT makes overviewMovieDataGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = overviewMovieDataGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'Userdata');
userData.MD.setNotes(get(handles.edit_notes, 'String'));
delete(handles.figure1);


% --- Executes on selection change in listbox_channel.
function listbox_channel_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_channel


% --- Executes during object creation, after setting all properties.
function listbox_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_notes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_notes as text
%        str2double(get(hObject,'String')) returns contents of edit_notes as a double


% --- Executes during object creation, after setting all properties.
function edit_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
