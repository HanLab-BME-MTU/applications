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

% Last Modified by GUIDE v2.5 11-May-2004 12:01:50

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
% uiwait(handles.NoiseParameterOptimizer);

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

% Attach a list of z values to handles.editProb
prob=[75 80 85 90 95 99];
set(handles.editProb,'UserData',prob);

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
% Discard the first line ("Select experiment")
exp=exp-1;

% Prepare input
noiseParams=zeros(1,5);
noiseParams(2:4)=fsmExpParam(exp).noiseParams; % [alpha beta I0]
noiseParams(5)=str2num(get(handles.editQuantile,'String'));
noiseParams(1)=noiseParams(5)/fsmExpParam(exp).gaussRatio;
bitDepth=fsmExpParam(exp).bitDepth;
prob=str2num(get(handles.editProb,'String'))/100;

% Loading images
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select first image');
if~((isa(fName,'char') & isa(dirName,'char')))
    return
end

% Get file names
outFileList=getFileStackNames([dirName,fName]);
n=length(outFileList);

% The user can decide the number of images to be cropped
prompt={'Specify the number of images to be used for the optimization.'};
dlg_title='User input requested';
num_lines=1;
default=3; % Default value is 3 frames
def={num2str(default)}; 
answer=fix(str2num(char(inputdlg(prompt,dlg_title,num_lines,def))));

% Check the selected number
if isempty(answer)
    disp('Aborting.\n');
    return
end

if answer<1 | answer>n
    fprintf(1,'Invalid number of images specified. Using the default value (%d).\n',default);
else
    % Crop outFileList
    n=answer;
    outFileList=outFileList(1:n);
end

% Starting optimization
[noiseParameter,actualP]=fsmAlphaBetaOptimization1D(outFileList,noiseParams,prob,bitDepth);

% Inform user
string=['Attained probability: ',num2str(100*actualP),'%'];
uiwait(msgbox(string,'Info','modal'));

% Print the result to console and ask the user to copy/paste the record into the database
fprintf(1,'\n\nTo add this experiment to your database:\n');
fprintf(1,'(1) Click on ''Edit experiment parameters'' in fsmCenter.\n');
fprintf(1,'(2) Copy/paste this record at the end of your experiment settings file.\n');
fprintf(1,'-------------------------------------------------------------------\n');
newLabel=[fsmExpParam(exp).label,' - OPTIMIZED FOR PROBABILITY = ',num2str(100*prob),'% (QUANTILE = ',num2str(noiseParams(5)),')'];
fprintf(1,'LABEL\t\t\t%s\n',newLabel);
fprintf(1,'DESCRIPTION\t\t%s\n',fsmExpParam(exp).description);
fprintf(1,'BIT DEPTH\t\t"%s"\n',num2str(bitDepth));
fprintf(1,'NOISE PARAMS\t"%1.8f %1.8f %1.8f %d"\n',noiseParameter(1),noiseParameter(2),noiseParameter(3),1); % '1' means that this is an optimized record
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
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',1);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editProb,'String',prob(1));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(1)/100)));

% --- Executes on button press in radioConf80.
function radioConf80_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf80
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',1);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editProb,'String',prob(2));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(2)/100)));


% --- Executes on button press in radioConf85.
function radioConf85_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf85
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',1);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editProb,'String',prob(3));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(3)/100)));

% --- Executes on button press in radioConf90.
function radioConf90_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf90
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',1);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);
set(handles.editProb,'String',prob(4));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(4)/100)));


% --- Executes on button press in radioConf95.
function radioConf95_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf95
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',1);
set(handles.radioConf99,'Value',0);
set(handles.editProb,'String',prob(5));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(5)/100)));


% --- Executes on button press in radioConf99.
function radioConf99_Callback(hObject, eventdata, handles)
% hObject    handle to radioConf99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioConf99
prob=get(handles.editProb,'UserData');
set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',1);
set(handles.editProb,'String',prob(6));
set(handles.editQuantile,'String',num2str(quantileFromConfProb(prob(6)/100)));


% --- Executes during object creation, after setting all properties.
function editProb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editProb_Callback(hObject, eventdata, handles)
% hObject    handle to editProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editProb as text
%        str2double(get(hObject,'String')) returns contents of editProb as a double

updateConfProbs(handles);



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
q=str2num(get(handles.editQuantile,'String'));
if q<=0
    q=1.96;
    set(handles.editQuantile,'String',num2str(1.96));
end
p=confProbFromQuantile(q);
set(handles.editProb,'String',num2str(100*p));
updateConfProbs(handles);

function updateConfProbs(handles)

set(handles.radioConf75,'Value',0);
set(handles.radioConf80,'Value',0);
set(handles.radioConf85,'Value',0);
set(handles.radioConf90,'Value',0);
set(handles.radioConf95,'Value',0);
set(handles.radioConf99,'Value',0);

prob=get(handles.editProb,'UserData');
userProb=str2num(get(handles.editProb,'String'));
if userProb>=100 | userProb<=0
    userProb=95;
end
indx=find(prob==userProb);
if ~isempty(indx)
    switch indx
        case 1, set(handles.radioConf75,'Value',1);
        case 2, set(handles.radioConf80,'Value',1);
        case 3, set(handles.radioConf85,'Value',1);
        case 4, set(handles.radioConf90,'Value',1);
        case 5, set(handles.radioConf95,'Value',1);
        case 6, set(handles.radioConf99,'Value',1);
        otherwise
            error('This should not happen');
    end
end
set(handles.editQuantile,'String',num2str(quantileFromConfProb(userProb/100)));


% --------------------------------------------------------------------
function menuExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close NoiseParameterOptimizer.
function NoiseParameterOptimizer_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to NoiseParameterOptimizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
hF=findall(0,'Tag','NoiseParameterOptimizer');
choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
switch choice,
    case 'Yes', delete(hF);
    case 'No', return;
end
