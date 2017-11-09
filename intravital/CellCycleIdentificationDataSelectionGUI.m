function varargout = CellCycleIdentificationDataSelectionGUI(varargin)
% CELLCYCLEIDENTIFICATIONDATASELECTIONGUI MATLAB code for CellCycleIdentificationDataSelectionGUI.fig
%      CELLCYCLEIDENTIFICATIONDATASELECTIONGUI, by itself, creates a new CELLCYCLEIDENTIFICATIONDATASELECTIONGUI or raises the existing
%      singleton*.
%
%      H = CELLCYCLEIDENTIFICATIONDATASELECTIONGUI returns the handle to a new CELLCYCLEIDENTIFICATIONDATASELECTIONGUI or the handle to
%      the existing singleton*.
%
%      CELLCYCLEIDENTIFICATIONDATASELECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLCYCLEIDENTIFICATIONDATASELECTIONGUI.M with the given input arguments.
%
%      CELLCYCLEIDENTIFICATIONDATASELECTIONGUI('Property','Value',...) creates a new CELLCYCLEIDENTIFICATIONDATASELECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellCycleIdentificationDataSelectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellCycleIdentificationDataSelectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellCycleIdentificationDataSelectionGUI

% Last Modified by GUIDE v2.5 11-Jul-2014 13:26:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellCycleIdentificationDataSelectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellCycleIdentificationDataSelectionGUI_OutputFcn, ...
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


% --- Executes just before CellCycleIdentificationDataSelectionGUI is made visible.
function CellCycleIdentificationDataSelectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellCycleIdentificationDataSelectionGUI (see VARARGIN)

% Choose default command line output for CellCycleIdentificationDataSelectionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellCycleIdentificationDataSelectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figCellCycleIdentificationDataSelectionGUI);


% --- Outputs from this function are returned to the command line.
function varargout = CellCycleIdentificationDataSelectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnOk.
function btnOk_Callback(hObject, eventdata, handles)
% hObject    handle to btnOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnSelectFucciCdt1File.
function btnSelectFucciCdt1File_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectFucciCdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editFUCCICdt1File_Callback(hObject, eventdata, handles)
% hObject    handle to editFUCCICdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFUCCICdt1File as text
%        str2double(get(hObject,'String')) returns contents of editFUCCICdt1File as a double


% --- Executes during object creation, after setting all properties.
function editFUCCICdt1File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFUCCICdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSelectFucciGemininFile.
function btnSelectFucciGemininFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editFucciGemininFile_Callback(hObject, eventdata, handles)
% hObject    handle to editFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFucciGemininFile as text
%        str2double(get(hObject,'String')) returns contents of editFucciGemininFile as a double


% --- Executes during object creation, after setting all properties.
function editFucciGemininFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSelectNuclearMarkerFile.
function btnSelectNuclearMarkerFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editNuclearMarkerFile_Callback(hObject, eventdata, handles)
% hObject    handle to editNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNuclearMarkerFile as text
%        str2double(get(hObject,'String')) returns contents of editNuclearMarkerFile as a double


% --- Executes during object creation, after setting all properties.
function editNuclearMarkerFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
