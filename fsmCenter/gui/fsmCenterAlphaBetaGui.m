function varargout = fsmCenterAlphaBetaGui(varargin)
% FSMCENTERALPHABETAGUI M-file for fsmCenterAlphaBetaGui.fig
%      FSMCENTERALPHABETAGUI, by itself, creates a new FSMCENTERALPHABETAGUI or raises the existing
%      singleton*.
%
%      H = FSMCENTERALPHABETAGUI returns the handle to a new FSMCENTERALPHABETAGUI or the handle to
%      the existing singleton*.
%
%      FSMCENTERALPHABETAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMCENTERALPHABETAGUI.M with the given input arguments.
%
%      FSMCENTERALPHABETAGUI('Property','Value',...) creates a new FSMCENTERALPHABETAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmCenterAlphaBetaGui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmCenterAlphaBetaGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmCenterAlphaBetaGui

% Last Modified by GUIDE v2.5 06-May-2004 17:56:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmCenterAlphaBetaGui_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmCenterAlphaBetaGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before fsmCenterAlphaBetaGui is made visible.
function fsmCenterAlphaBetaGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmCenterAlphaBetaGui (see VARARGIN)

% Choose default command line output for fsmCenterAlphaBetaGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fsmCenterAlphaBetaGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Read parameter experiments from fsmExpParams.txt
userDir=fsmCenter_getUserSettings;
if isempty(userDir)
    fsmExpParamPath=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
else
    fsmExpParamPath=[userDir,filesep,'fsmExpParams.txt'];
    if exist(fsmExpParamPath)~=2
        % No database found in user-defined directory
        % Reverting to default database
        fsmExpParamPath=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
    end
end
if exist(fsmExpParamPath)~=2
    uiwait(msgbox('Could not find experiment database! Get fsmExpParams.txt from the repository and restart SpeckTackle.','Error','modal'));
    return
end

% Fill parameters structure
fsmExpParam=fsmGuiScanDataBase(fsmExpParamPath);

% Fill the scroll-down menu in the user interface
labels=cell(length(fsmExpParam)+1,1);
labels(1,1)={'Select experiment'};
for i=1:length(fsmExpParam)
    labels(i+1,1)={fsmExpParam(i).label};
end
set(handles.expPopup,'String',labels);

% Attach the experiment database to expPopup userdata
set(handles.expPopup,'UserData',fsmExpParam);

% Attach a list of z values to handles.editQuantile
zValues=[1.15 1.29 1.45 1.645 1.96 2.58];
set(handles.editQuantile,'UserData',zValues);

% --- Outputs from this function are returned to the command line.
function varargout = fsmCenterAlphaBetaGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function expPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in expPopup.
function expPopup_Callback(hObject, eventdata, handles)
% hObject    handle to expPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns expPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from expPopup
fsmExpParam=get(handles.expPopup,'UserData');
exp=get(handles.expPopup,'Value');
% Update description
if exp==1
    set(handles.textDescr,'String','Experiment description');
else
    set(handles.textDescr,'String',fsmExpParam(exp-1).description);
end


% --- Executes on button press in pushStart.
function pushStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read noise parameters for chosen experiment
fsmExpParam=get(handles.expPopup,'UserData');
exp=get(handles.expPopup,'Value');
if exp==1
    uiwait(msgbox('Please select an experiment you want to optimize.','Error','modal'));
    return
end
noiseParams=zeros(1,5);
noiseParams(2:4)=fsmExpParam(exp).noiseParams; % [alpha beta I0]
noiseParams(5)=str2num(get(handles.editQuantile,'String'));
noiseParams(1)=noiseParams(5)/fsmExpParam(exp).gaussRatio;
bitDepth=fsmExpParam(exp).bitDepth;

% Startin optimization
% [noiseParameter,actualP]=fsmAlphaBetaOptimization1D(optimImageStack,noiseParams,noiseParams(5),pixelDepth);
noiseParameter=noiseParams;

fprintf(1,'\n\nTo add this experiment to your database:\n');
fprintf(1,'(1) Click on ''Edit experiment parameters'' in fsmCenter.\n');
fprintf(1,'(2) Copy/paste this record at the end of your experiment settings file.\n');
fprintf(1,'-------------------------------------------------------------------\n');
newLabel=[fsmExpParam(exp).label,' - OPTIMIZED FOR QUANTILE = ',num2str(noiseParams(5))];
fprintf(1,'LABEL\t\t\t%s\n',newLabel);
fprintf(1,'DESCRIPTION\t\t%s\n',fsmExpParam(exp).description);
fprintf(1,'BIT DEPTH\t\t"%s"\n',num2str(bitDepth));
fprintf(1,'NOISE PARAMS\t"%1.8f %1.8f %1.8f"\n',noiseParameter(1),noiseParameter(2),noiseParameter(3));
fprintf(1,'GAUSS RATIO\t\t"%1.2f"\n',fsmExpParam(exp).gaussRatio);
fprintf(1,'#\n');
fprintf(1,'-------------------------------------------------------------------\n');
fprintf(1,'(Don''t forget the ''#'')\n');


% --- Executes on button press in radioConf75.
function radioConf75_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf75
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',1);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editQuantile,'String',zValues(1));

% --- Executes on button press in radioConf80.
function radioConf80_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf80
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',1);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editQuantile,'String',zValues(2));


% --- Executes on button press in radioConf85.
function radioConf85_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf85
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',1);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editQuantile,'String',zValues(3));


% --- Executes on button press in radioConf90.
function radioConf90_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf90
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',1);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editQuantile,'String',zValues(4));


% --- Executes on button press in radioConf95.
function radioConf95_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf95
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',1);
set(handles.radioConf99,'Value',0);
set(handles.editQuantile,'String',zValues(5));


% --- Executes on button press in radioConf99.
function radioConf99_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf99
zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',1);
set(handles.editQuantile,'String',zValues(6));


% --- Executes during object creation, after setting all properties.
function editQuantile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editQuantile_Callback(hObject, eventdata, handles)
% hObject    handle to editQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editQuantile as text
%        str2double(get(hObject,'String')) returns contents of editQuantile as a double

zValues=get(handles.editQuantile,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',1);
set(handles.editQuantile,'String',zValues(6));

