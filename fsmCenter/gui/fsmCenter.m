function varargout = fsmCenter(varargin)
% FSMCENTER M-file for fsmCenter.fig
%      FSMCENTER, by itself, creates a new FSMCENTER or raises the existing
%      singleton*.
%
%      H = FSMCENTER returns the handle to a new FSMCENTER or the handle to
%      the existing singleton*.
%
%      FSMCENTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMCENTER.M with the given input arguments.
%
%      FSMCENTER('Property','Value',...) creates a new FSMCENTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmCenter_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmCenter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmCenter

% Last Modified by GUIDE v2.5 27-Aug-2004 11:19:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmCenter_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmCenter_OutputFcn, ...
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

% --- Executes just before fsmCenter is made visible.
function fsmCenter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmCenter (see VARARGIN)

% Choose default command line output for fsmCenter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fsmCenter wait for user response (see UIRESUME)
% uiwait(handles.fsmCenter);

% Checks for user settings
userDir=fsmCenter_getUserSettings;
if isempty(userDir)
    uiwait(msgbox('No user directory has been defined yet. Please click on "User settings" and follow the instructions.','Info','modal'));
end

% --- Outputs from this function are returned to the command line.
function varargout = fsmCenter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% createFcn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function editBitDepth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBitDepth (see GCBO)
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
function editSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSigma (see GCBO)
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
function editFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFileName (see GCBO)
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
function editY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
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
function editX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
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
function editTSampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTSampling (see GCBO)
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
function editLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLabel (see GCBO)
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
function editSpectrumSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSpectrumSigma (see GCBO)
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
function editFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrame (see GCBO)
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
function editScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
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
function editD0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD0 (see GCBO)
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
function editImageNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImageNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Callbacks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushFSM_Callback(hObject, eventdata, handles)
fsmGuiMain;

function pushCrop_Callback(hObject, eventdata, handles)
y0S=get(handles.editCalY0,'String');
yS=get(handles.editCalY,'String');
x0S=get(handles.editCalX0,'String');
xS=get(handles.editCalX,'String');
if any([isempty(y0S) isempty(yS) isempty(x0S) isempty(xS)])
    area=[];
else
    y0=fix(str2num(y0S));
    y=fix(str2num(yS));
    x0=fix(str2num(x0S));
    x=fix(str2num(xS));
    if (y<=y0 | x<=x0 | y0<0 | x0<0) % We don't know the size of the image yet
        uiwait(msgbox('Invalid region specified. You will need to draw.','error','modal'));
        area=[];
    else
        area=[y0 x0 y x];
    end
end
retArea=cropStack(area);
if isempty(retArea)
    % Interrupted by user
    return
end
if isempty(area)
    msg=['Selected area: y = ',num2str(retArea(1)),' : ',num2str(retArea(3)),'; x = ',num2str(retArea(2)),' : ',num2str(retArea(4))];
    choice=myQuestdlg(msg,'Info','Fill','Discard','Fill');
    if strcmp(choice,'Fill')
        set(handles.editCalY0,'String',num2str(retArea(1)));
        set(handles.editCalY,'String',num2str(retArea(3)));
        set(handles.editCalX0,'String',num2str(retArea(2)));
        set(handles.editCalX,'String',num2str(retArea(4)));
    end
end

function pushCal_Callback(hObject, eventdata, handles)
userDir=fsmCenter_getUserSettings;
if ~isempty(userDir)
    dataFile=[userDir,filesep,'fsmExpParams.txt'];
else
    dataFile=which('fsmExpParams.txt'); % Default fsmExpParams.txt file
    % Check
    fsmMainPath=which('fsmMain.m');
    [pathExpParams,fname,no,ext]=getFilenameBody(dataFile);
    [pathFsmMain,fname,no,ext]=getFilenameBody(fsmMainPath);
    if ~strcmp(pathExpParams,pathFsmMain)
        uiwait(msgbox('Something is wrong with your configuration. Plase click on "User settings" in fsmCenter and follow the instructions. Then try again.','Error','modal'));
        return
    end
end
[I0,sDN,GaussRatio,status]=fsmCalcNoiseParam([],str2num(get(handles.editBitDepth,'string')),str2num(get(handles.editSigma,'string')),dataFile);
if status==-1
    uiwait(msgbox('The noise parameters could not be saved to your user database. They were returned to the MATLAB console. Make sure that a file called ''fsmExpParams.txt'' exists in your user directory or in the root directory of the qFSM package. Click on ''User settings'' in fsmCenter and follow the instructions.','Info','modal'));
elseif status==0
    disp('Calibration interrupted by the user.');
elseif status==1
    uiwait(msgbox('The noise parameters were successfully added to your user database. They will be added to the ''Camera calibration'' popup menu the next time you''ll start SpeckTackle. If you need to modify the entered record, click on ''Edit experiment database'' in fsmCenter and scroll down to the end to the file.','Info','modal'));
else
    error('The function fsmCalcNoiseParam returned an invalid status.');
end

function editBitDepth_Callback(hObject, eventdata, handles)

function editSigma_Callback(hObject, eventdata, handles)

function editFileName_Callback(hObject, eventdata, handles)

function pushSpectral_Callback(hObject, eventdata, handles)
analysis=[get(handles.checkPowerSpectrum,'Value') get(handles.checkMeanScores,'Value')];
if analysis==[0 0]
    uiwait(msgbox('Please select some analysis to perform.','Info','modal'));
    return
end
moveROI=get(handles.checkMoveROI,'Value');
choice=find([get(handles.radioGrid,'Value') get(handles.radioPoly,'Value') get(handles.radioRect,'Value')]);
switch choice
    case 1, roi='g';
    case 2, roi='p';
    case 3, roi='r';
    otherwise error('bad selection for ROI');
end
gridSize=[str2num(get(handles.editY,'String')) str2num(get(handles.editX,'String'))];
sigma=str2num(get(handles.editTimeSigma,'String'));
sigmaFreq=str2num(get(handles.editSpectrumSigma,'String'));
tSampling=str2num(get(handles.editTSampling,'String'));
label=str2num(get(handles.editLabel,'String'));
if roi=='g' % The b/w mask is used only to get a subset of the grid
    if get(handles.checkLoadBwMask,'Value')==1
        % Load wavesBwMask.mat
        [fName,dirName] = uigetfile(...
            {'wavesBwMask.mat;','wavesBwMask.mat';
            '*.*','All Files (*.*)'},...
            'Select wavesBwMask.mat');
        if ~(isa(fName,'char') & isa(dirName,'char'))
            bwMask=[];
        else
            load([dirName,filesep,fName]);
            if exist('bwMask')~=1
                uiwait(errordlg('The chosen file does not appear to contain a valid b/w mask.','Error','modal'));
                return;
            end
        end
    else
        bwMask=[];
    end
else
    bwMask=[];
end
[scores,allScores,bwMask,py,px]=fsmWaveAnalysis(analysis,roi,gridSize,sigma,sigmaFreq,label,tSampling,moveROI,bwMask);
if ~isempty(bwMask)
    assignin('base','scores',scores);
    assignin('base','allScores',allScores);
    assignin('base','bwMask',bwMask);
    assignin('base','py',py);
    assignin('base','px',px);
    if roi=='g' % The b/w mask is used only to get a subset of the grid
        if get(handles.checkSaveBwMask,'Value')==1
            path=uigetdir('','Please select a directory where to save the b/w mask'); 
            if path~=0
                save([path,filesep,'wavesBwMask.mat'],'bwMask');
            end
        end
    end
end

function editY_Callback(hObject, eventdata, handles)

function editX_Callback(hObject, eventdata, handles)

function editTSampling_Callback(hObject, eventdata, handles)
set(handles.editTSampling,'String',num2str(abs(str2num(get(handles.editTSampling,'String')))));

function editLabel_Callback(hObject, eventdata, handles)
if str2num(get(handles.editLabel,'String'))<0 | str2num(get(handles.editLabel,'String'))>1
    uiwait(msgbox('Please enter a value between 0 and 1.','Help','modal'));
    set(handles.editLabel,'String','0.3')
end

function radioGrid_Callback(hObject, eventdata, handles)
set(handles.radioGrid,'Value',1);
set(handles.radioPoly,'Value',0);
set(handles.radioRect,'Value',0);
%set(handles.checkMoveROI,'Enable','off');
set(handles.textGrid,'Enable','on');
set(handles.textY,'Enable','on');
set(handles.textX,'Enable','on');
set(handles.editY,'Enable','on');
set(handles.editX,'Enable','on');
set(handles.textBwMask,'Enable','on');
set(handles.checkLoadBwMask,'Enable','on');
set(handles.checkSaveBwMask,'Enable','on');

function radioPoly_Callback(hObject, eventdata, handles)
set(handles.radioGrid,'Value',0);
set(handles.radioPoly,'Value',1);
set(handles.radioRect,'Value',0);
%set(handles.checkMoveROI,'Enable','on');
set(handles.textGrid,'Enable','off');
set(handles.textY,'Enable','off');
set(handles.textX,'Enable','off');
set(handles.editY,'Enable','off');
set(handles.editX,'Enable','off');
set(handles.textBwMask,'Enable','off');
set(handles.checkLoadBwMask,'Enable','off');
set(handles.checkSaveBwMask,'Enable','off');

function radioRect_Callback(hObject, eventdata, handles)
set(handles.radioGrid,'Value',0);
set(handles.radioPoly,'Value',0);
set(handles.radioRect,'Value',1);
%set(handles.checkMoveROI,'Enable','on');
set(handles.textGrid,'Enable','off');
set(handles.textY,'Enable','off');
set(handles.textX,'Enable','off');
set(handles.editY,'Enable','off');
set(handles.editX,'Enable','off');
set(handles.textBwMask,'Enable','off');
set(handles.checkLoadBwMask,'Enable','off');
set(handles.checkSaveBwMask,'Enable','off');

function pushMeta_Callback(hObject, eventdata, handles)
stk2tif(get(handles.editFileName,'string'));

function editSpectrumSigma_Callback(hObject, eventdata, handles)
% Make sure the entered value is > 0
set(handles.editSpectrumSigma,'String',num2str(abs(fix(str2num(get(handles.editSpectrumSigma,'String'))))));
if ~mod(str2num(get(handles.editSpectrumSigma,'String')),2)
    uiwait(msgbox('Please enter an integer, ODD value > 0.','Help','modal'));
    set(handles.editSpectrumSigma,'String','1');
end

function pushUser_Callback(hObject, eventdata, handles)
fsmCenter_setUserSettings;

function pushVector_Callback(hObject, eventdata, handles)
nAvg=str2num(get(handles.editVectorAnalysis,'String'));
displ=[get(handles.checkRaw,'Value') get(handles.checkInterp,'Value') get(handles.checkNoise,'Value') get(handles.checkError,'Value') get(handles.checkDisplayImg,'Value')];
roi=[get(handles.checkROI,'Value') get(handles.checkLoadROI,'Value') get(handles.checkSaveROI,'Value')];
scale=str2num(get(handles.editScale,'String'));
d0=str2num(get(handles.editD0,'String'));
useDiv=get(handles.checkDiv,'Value');
output=find([get(handles.radioCMap,'Value') get(handles.radioCircle,'Value')]);
displROI=get(handles.checkDisplayROI,'Value');
invertROI=get(handles.checkInvertROI,'Value');
fsmVectorAnalysis(nAvg,roi,displ,scale,d0,useDiv,output,displROI,invertROI);

function editFrame_Callback(hObject, eventdata, handles)

function checkAutoScale_Callback(hObject, eventdata, handles)

function editScale_Callback(hObject, eventdata, handles)

function editD0_Callback(hObject, eventdata, handles)

function checkDiv_Callback(hObject, eventdata, handles)

function checkDisplayImg_Callback(hObject, eventdata, handles)

function checkRaw_Callback(hObject, eventdata, handles)

function checkInterp_Callback(hObject, eventdata, handles)

function checkNoise_Callback(hObject, eventdata, handles)

function checkError_Callback(hObject, eventdata, handles)
if get(handles.checkError,'Value')==1
    set(handles.radioCircle,'Enable','on');
    set(handles.radioCMap,'Enable','on');
else
    set(handles.radioCircle,'Enable','off');
    set(handles.radioCMap,'Enable','off');
end
    
function radioCMap_Callback(hObject, eventdata, handles)
set(handles.radioCMap,'Value',1);
set(handles.radioCircle,'Value',0);

function radioCircle_Callback(hObject, eventdata, handles)
set(handles.radioCMap,'Value',0);
set(handles.radioCircle,'Value',1);

function checkROI_Callback(hObject, eventdata, handles)
set(handles.checkLoadROI,'Value',0);

function checkSaveROI_Callback(hObject, eventdata, handles)
if get(handles.checkLoadROI,'Value')==1 | get(handles.checkROI,'Value')==0
    set(handles.checkSaveROI,'Value',0);
end

function checkLoadROI_Callback(hObject, eventdata, handles)
set(handles.checkROI,'Value',0);
set(handles.checkSaveROI,'Value',0);

function checkdisplayROI_Callback(hObject, eventdata, handles)

function pushBatchJobs_Callback(hObject, eventdata, handles)
fsmBatchJobs;

function pushImShow_Callback(hObject, eventdata, handles)
figure;fsmCenterCB_loadNewImage;

function checkDisplayROI_Callback(hObject, eventdata, handles)

function checkInvertROI_Callback(hObject, eventdata, handles)

function pushEditExpParams_Callback(hObject, eventdata, handles)
userDir=fsmCenter_getUserSettings;
if ~isempty(userDir)
    edit([userDir,filesep,'fsmExpParams.txt']);
else
    edit('fsmExpParams.txt'); % Default fsmExpParams.txt file
end

function editImageNumber_Callback(hObject, eventdata, handles)

function pushScalarLow_Callback(hObject, eventdata, handles)

function pushClearSA_Callback(hObject, eventdata, handles)

function checkMoveROI_Callback(hObject, eventdata, handles)

function pushHelpKin_Callback(hObject, eventdata, handles)
%uiwait(msgbox('Sorry, help not available yet.','Help','modal'));
web(['file:///' which('qFSM_default.html')]);

function pushHelpFlow_Callback(hObject, eventdata, handles)
%uiwait(msgbox('Sorry, help not available yet.','Help','modal'));
web(['file:///' which('qFSM_default.html')]);

function checkPowerSpectrum_Callback(hObject, eventdata, handles)
if get(handles.checkPowerSpectrum,'Value')==1
    set(handles.textSpectrumSigma,'Enable','on');
    set(handles.editSpectrumSigma,'Enable','on');
    set(handles.textTSampling,'Enable','on');
    set(handles.editTSampling,'Enable','on');
    set(handles.textLabel,'Enable','on');
    set(handles.editLabel,'Enable','on');
else
    set(handles.textSpectrumSigma,'Enable','off');
    set(handles.editSpectrumSigma,'Enable','off');
    set(handles.textTSampling,'Enable','off');
    set(handles.editTSampling,'Enable','off');
    set(handles.textLabel,'Enable','off');
    set(handles.editLabel,'Enable','off');
end

function checkMeanScores_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function EditTimeSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditTimeSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function EditTimeSigma_Callback(hObject, eventdata, handles)
% hObject    handle to EditTimeSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditTimeSigma as text
%        str2double(get(hObject,'String')) returns contents of EditTimeSigma as a double


% --------------------------------------------------------------------
function mainMenu_Callback(hObject, eventdata, handles)
% hObject    handle to mainMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function editTimeSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTimeSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editTimeSigma_Callback(hObject, eventdata, handles)
% hObject    handle to editTimeSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTimeSigma as text
%        str2double(get(hObject,'String')) returns contents of editTimeSigma as a double

% --- Executes on button press in pushHelpSpeed.
function pushSpeedMaps_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkSMMask,'Value')==1
    % Load userROI.mat
    [fName,dirName] = uigetfile(...
        {'userROI.mat;','userROI.mat';
        '*.*','All Files (*.*)'},...
        'Select userROI.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        choice=questdlg('You didn''t pick any file.','Error','Continue without ROI','Exit','Exit');
        switch choice,
            case 'Continue without ROI'; userROIbw=[]; filePicked=0;
            case 'Exit', return;
        end % switch
    else
        filePicked=1;
    end
    if filePicked==1
        % Load the file
        load([dirName,fName]);
    end
    if filePicked==1 & exist('userROIbw')~=1
        choice=questdlg('The loaded file does not seem to contain a valid ROI.','Error','Continue without ROI','Exit','Exit');
        switch choice,
            case 'Continue without ROI', userROIbw=[];
            case 'Exit', return;
        end % switch
    end
else
    userROIbw=[];
end
gridSize=[str2num(get(handles.editYSP,'string')) str2num(get(handles.editXSP,'string'))];
n=str2num(get(handles.editFrameSP,'string'));
d0=str2num(get(handles.editSMCL,'string'));
loadMPM=get(handles.radioSMMPM,'value');
overlayVectors=get(handles.checkSMOverlay,'value');
sampling=str2num(get(handles.editSMSampling,'string'));
pixelSize=str2num(get(handles.editSMPixel,'string'));
maxSpeed=str2num(get(handles.editSMMax,'string'));
segment=get(handles.checkSMSegment,'value');
bitDepth=str2num(get(handles.editSMBitDepth,'String'));
outputdir=fsmSpeedMaps(gridSize,n,d0,loadMPM,sampling,pixelSize,overlayVectors,userROIbw,maxSpeed,segment,bitDepth);
if ~isempty(outputdir)
    % Maps have been saved to disk
    msg=['Speed maps have been saved to ',outputdir,'.'];
    uiwait(msgbox(msg,'Help','modal'));
end

% --- Executes on button press in pushHelpSpeed.
function pushHelpSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiwait(msgbox('Sorry, help not available yet.','Help','modal'));
web(['file:///' which('qFSM_default.html')]);

% --- Executes during object creation, after setting all properties.
function editYSP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editYSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editYSP_Callback(hObject, eventdata, handles)
% hObject    handle to editYSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editYSP as text
%        str2double(get(hObject,'String')) returns contents of editYSP as a double


% --- Executes during object creation, after setting all properties.
function editXSP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editXSP_Callback(hObject, eventdata, handles)
% hObject    handle to editXSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXSP as text
%        str2double(get(hObject,'String')) returns contents of editXSP as a double


% --- Executes during object creation, after setting all properties.
function editFrameSP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editFrameSP_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameSP as text
%        str2double(get(hObject,'String')) returns contents of editFrameSP as a double
nAvg=fix(str2num(get(handles.editFrameSP,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
end
set(handles.editFrameSP,'String',num2str(nAvg));


% --- Executes during object creation, after setting all properties.
function editSMCL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSMCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSMCL_Callback(hObject, eventdata, handles)
% hObject    handle to editSMCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSMCL as text
%        str2double(get(hObject,'String')) returns contents of editSMCL as a double


% --- Executes on button press in radioSMMPM.
function radioSMMPM_Callback(hObject, eventdata, handles)
% hObject    handle to radioSMMPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioSMMPM
set(handles.radioSMMPM,'Value',1);
set(handles.radioSMMD,'Value',0);

% --- Executes on button press in radioSMMD.
function radioSMMD_Callback(hObject, eventdata, handles)
% hObject    handle to radioSMMD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioSMMD
set(handles.radioSMMPM,'Value',0);
set(handles.radioSMMD,'Value',1);


% --- Executes on button press in checkSMPoly.
function checkSMPoly_Callback(hObject, eventdata, handles)
% hObject    handle to checkSMPoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSMPoly


% --- Executes during object creation, after setting all properties.
function editSMSampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSMSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSMSampling_Callback(hObject, eventdata, handles)
% hObject    handle to editSMSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSMSampling as text
%        str2double(get(hObject,'String')) returns contents of editSMSampling as a double


% --- Executes during object creation, after setting all properties.
function editSMPixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSMPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSMPixel_Callback(hObject, eventdata, handles)
% hObject    handle to editSMPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSMPixel as text
%        str2double(get(hObject,'String')) returns contents of editSMPixel as a double


% --- Executes during object creation, after setting all properties.
function editSMMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSMMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSMMax_Callback(hObject, eventdata, handles)
% hObject    handle to editSMMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSMMax as text
%        str2double(get(hObject,'String')) returns contents of editSMMax as a double


% --- Executes on button press in checkSMOverlay.
function checkSMOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to checkSMOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSMOverlay


% --- Executes on button press in pushPRStart.
function pushPRStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushPRStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushPRHelp.
function pushPRHelp_Callback(hObject, eventdata, handles)
% hObject    handle to pushPRHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushPlotTraj.
function pushPlotTraj_Callback(hObject, eventdata, handles)
% hObject    handle to pushPlotTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
method=find([get(handles.radioPTMethod1,'Value') get(handles.radioPTMethod2,'Value')]);
colorCode=get(handles.checkTrajPlots,'Value');
fsmPlotTrajectories(method,colorCode);

% --- Executes on button press in pushHelpPlotTraj.
function pushHelpPlotTraj_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpPlotTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web(['file:///' which('qFSM_plotTrajectories.html')]);


% --- Executes on button press in radioPTMethod1.
function radioPTMethod1_Callback(hObject, eventdata, handles)
% hObject    handle to radioPTMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioPTMethod1
set(handles.radioPTMethod1,'Value',1);
set(handles.radioPTMethod2,'Value',0);

% --- Executes on button press in radioPTMethod2.
function radioPTMethod2_Callback(hObject, eventdata, handles)
% hObject    handle to radioPTMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioPTMethod2
set(handles.radioPTMethod1,'Value',0);
set(handles.radioPTMethod2,'Value',1);


% --- Executes during object creation, after setting all properties.
function editCalY0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCalY0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editCalY0_Callback(hObject, eventdata, handles)
% hObject    handle to editCalY0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCalY0 as text
%        str2double(get(hObject,'String')) returns contents of editCalY0 as a double


% --- Executes during object creation, after setting all properties.
function editCalY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCalY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editCalY_Callback(hObject, eventdata, handles)
% hObject    handle to editCalY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCalY as text
%        str2double(get(hObject,'String')) returns contents of editCalY as a double


% --- Executes during object creation, after setting all properties.
function editCalX0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCalX0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editCalX0_Callback(hObject, eventdata, handles)
% hObject    handle to editCalX0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCalX0 as text
%        str2double(get(hObject,'String')) returns contents of editCalX0 as a double


% --- Executes during object creation, after setting all properties.
function editCalX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCalX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editCalX_Callback(hObject, eventdata, handles)
% hObject    handle to editCalX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCalX as text
%        str2double(get(hObject,'String')) returns contents of editCalX as a double


% --- Executes on button press in pushKinMaps.
function pushKinMaps_Callback(hObject, eventdata, handles)
% hObject    handle to pushKinMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load the first image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Load image');
if(isa(fName,'char') & isa(dirName,'char'))
    
    % Image information - needed for image size
    info=imfinfo([dirName,fName]);
    
else
    return % Returns an error (status=0)
end

% Read parameters
n=str2num(get(handles.editFrameTN,'String'));
sigma=str2num(get(handles.editSigmaTN,'String'));

% Call function
[polyMap,depolyMap,netMapRGB,outputdir]=fsmKineticMaps([],[info.Height info.Width],[-1 n],sigma);

if isempty(outputdir)
    
    if ~isempty(polyMap) & ~isempty(depolyMap) & ~isempty(netMapRGB)
        % Show maps and return them to Matlab base workspace
        figure; imshow(polyMap,[]); title('Polymerization map'); assignin('base','polyMap',polyMap);
        figure; imshow(-depolyMap,[]); title('Depolymerization map'); assignin('base','depolyMap',depolyMap);
        figure; imshow(netMapRGB,[]); title('Net assembly map'); assignin('base','netMapRGB',netMapRGB);
        assignin('base','netMap',polyMap+depolyMap);
    end
    
else
    
    if ~isempty(polyMap) & ~isempty(depolyMap) & ~isempty(netMapRGB)
        % Maps have been saved to disk
        msg=['Turnover maps have been saved to ',outputdir,'.'];
        uiwait(msgbox(msg,'Help','modal'));
    else
        % Interrupted by user in the last step
    end
    
end


% --- Executes on button press in pushHelpTurnover.
function pushHelpTurnover_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpTurnover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web(['file:///' which('qFSM_default.html')]);


% --- Executes during object creation, after setting all properties.
function editFrameTN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editFrameTN_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameTN as text
%        str2double(get(hObject,'String')) returns contents of editFrameTN as a double
nAvg=fix(str2num(get(handles.editFrameTN,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
end
set(handles.editFrameTN,'String',num2str(nAvg));

% --- Executes during object creation, after setting all properties.
function editSigmaTN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSigmaTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSigmaTN_Callback(hObject, eventdata, handles)
% hObject    handle to editSigmaTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSigmaTN as text
%        str2double(get(hObject,'String')) returns contents of editSigmaTN as a double


% --- Executes on button press in pushCalReset.
function pushCalReset_Callback(hObject, eventdata, handles)
% hObject    handle to pushCalReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.editCalY0,'String','');
set(handles.editCalY,'String','');
set(handles.editCalX0,'String','');
set(handles.editCalX,'String','');

% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close fsmCenter.
function fsmCenter_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fsmCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
fsmH=findall(0,'Tag','fsmCenter'); % Get the handle of fsmCenter
choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
switch choice,
    case 'Yes', delete(fsmH);
    case 'No', return;
end % switch


% --- Executes during object deletion, before destroying properties.
function fsmCenter_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fsmCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function helpFsmCenter_Callback(hObject, eventdata, handles)
% hObject    handle to helpFsmCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% stat=web(['file:///' which('fsmCenterHelp.html') ],'-browser');
% stat=web(['file:///' which('fsmCenterHelp.html') ]);


% --- Executes on button press in checkTrajPlots.
function checkTrajPlots_Callback(hObject, eventdata, handles)
% hObject    handle to checkTrajPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkTrajPlots


% --- Executes on button press in checkSMMask.
function checkSMMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkSMMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSMMask


% --- Executes on button press in pushConvertSA.
function pushConvertSA_Callback(hObject, eventdata, handles)
% hObject    handle to pushConvertSA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load speckleArray.mat
[fName,dirName] = uigetfile(...
    {'speckleArray.mat;','speckleArray.mat';
    '*.*','All Files (*.*)'},...
    'Select speckleArray.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return
end
load([dirName,fName]);
if exist('speckleArray')~=1
    % The loaded speckleArray is invalid
    uiwait(msgbox('The loaded speckleArray is invalid.','error','modal'));
    return
end
% Convert
[speckleArray,status]=fsmConvertSpeckleArray(speckleArray);
if isempty(speckleArray)
    if status==1
        % The loaded speckleArray is invalid
        uiwait(msgbox('The loaded speckleArray is invalid.','error','modal'));
        return
    elseif status==2
        % The loaded speckleArray was already in the new format
        uiwait(msgbox('The loaded speckleArray is already in the new format.','Info','modal'));
        return
    else
        error('if speckleArray is empty, status has to be either 1 or 2');
    end
end
% Save
quit=0;
rounds=0;
while quit==0
    rounds=rounds+1;
    if rounds==1
        path=uigetdir(dirName,'Please select a directory where to save the new speckleArray');
    else
        path=uigetdir(dirName,'Please specifiy a different directory or leave whithout saving');
    end
    if path~=0
        if exist([path,filesep,'speckleArray.mat'])==2
            choice=questdlg('A speckleArray.mat file already exists in this directory. Do you want to overwrite it?','User input requested','Yes','No','Yes');
            switch choice,
                case 'Yes', save([path,filesep,'speckleArray.mat'],'speckleArray'); quit=1; uiwait(msgbox('Converted speckleArray.mat successfully saved.','Info','modal'));
                case 'No', % Do nothing
            end % switch
        else
            save([path,filesep,'speckleArray.mat'],'speckleArray'); quit=1; uiwait(msgbox('Converted speckleArray.mat successfully saved.','Info','modal'));
        end
    else
        disp('Aborted.');
        quit=1;
    end
end




% --------------------------------------------------------------------
function menuToolsCalib_Callback(hObject, eventdata, handles)
% hObject    handle to menuToolsCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuToolsConvert_Callback(hObject, eventdata, handles)
% hObject    handle to menuToolsConvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushOptimCal.
function pushOptimCal_Callback(hObject, eventdata, handles)
% hObject    handle to pushOptimCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fsmCenterAlphaBetaGui;


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkLoadBwMask.
function checkLoadBwMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkLoadBwMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkLoadBwMask


% --- Executes on button press in checkSaveBwMask.
function checkSaveBwMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkSaveBwMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSaveBwMask


% --- Executes on button press in checkSMSegment.
function checkSMSegment_Callback(hObject, eventdata, handles)
% hObject    handle to checkSMSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSMSegment
if get(handles.checkSMSegment,'Value')==1
    set(handles.textSMBitDepth,'Enable','on');
    set(handles.editSMBitDepth,'Enable','on');
else
    set(handles.textSMBitDepth,'Enable','off');
    set(handles.editSMBitDepth,'Enable','off');
end    

% --------------------------------------------------------------------
function MenuToolsCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to MenuToolsCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuToolsOptCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to MenuToolsOptCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function editSMBitDepth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSMBitDepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSMBitDepth_Callback(hObject, eventdata, handles)
% hObject    handle to editSMBitDepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSMBitDepth as text
%        str2double(get(hObject,'String')) returns contents of editSMBitDepth as a double




% --- Executes on button press in pushAssignSpecklesToClasses.
function pushAssignSpecklesToClasses_Callback(hObject, eventdata, handles)
% hObject    handle to pushAssignSpecklesToClasses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load speckleArray.mat
[fName,dirName] = uigetfile(...
    {'speckleArray.mat;','speckleArray.mat';
    '*.*','All Files (*.*)'},...
    'Select speckleArray.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return
end
load([dirName,fName]);
if exist('speckleArray')~=1
    % The loaded speckleArray is invalid
    uiwait(msgbox('The loaded speckleArray is invalid.','error','modal'));
    return
end

% Check that the speckleArray is not in the old format
if ~(length(speckleArray)==1 & length(speckleArray(1).timepoint)>1)
    disp('This speckleArray structure is in the old format. Use ''Convert speckleArray'' in fsmCenter to update it and then retry.');
    return
end

% Assign speckles to classes
speckleClasses=fsmAssignSpecklesToClasses(speckleArray);

% Save
quit=0;
rounds=0;
while quit==0
    rounds=rounds+1;
    if rounds==1
        path=uigetdir(dirName,'Please select a directory where to save speckleClasses.mat');
    else
        path=uigetdir(dirName,'Please specifiy a different directory or click on cancel to leave whithout saving');
    end
    if path~=0
        if exist([path,filesep,'speckleClasses.mat'])==2
            choice=questdlg('A speckleClasses.mat file already exists in this directory. Do you want to overwrite it?','User input requested','Yes','No','Yes');
            switch choice,
                case 'Yes', save([path,filesep,'speckleClasses.mat'],'speckleClasses'); quit=1; uiwait(msgbox('speckleClasses.mat successfully saved.','Info','modal'));
                case 'No', % Do nothing
            end % switch
        else
            save([path,filesep,'speckleClasses.mat'],'speckleClasses'); quit=1; uiwait(msgbox('speckleClasses.mat successfully saved.','Info','modal'));
        end
    else
        disp('Aborted.');
        quit=1;
    end
end


% --- Executes on button press in pushHelpAssignClasses.
function pushHelpAssignClasses_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpAssignClasses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web(['file:///' which('qFSM_speckleClasses.html')]);


% --- Executes on button press in pushHelpTurnoverMaps.
function pushHelpTurnoverMaps_Callback(hObject, eventdata, handles)
% hObject    handle to pushHelpTurnoverMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web(['file:///' which('qFSM_turnoverMaps.html')]);




% --- Executes during object creation, after setting all properties.
function editVectorAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVectorAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editVectorAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to editVectorAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVectorAnalysis as text
%        str2double(get(hObject,'String')) returns contents of editVectorAnalysis as a double
nAvg=fix(str2num(get(handles.editVectorAnalysis,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
end
set(handles.editVectorAnalysis,'String',num2str(nAvg));



