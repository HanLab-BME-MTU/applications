function varargout = fsmPostProc(varargin)
% FSMPOSTPROC M-file for fsmPostProc.fig
%      FSMPOSTPROC, by itself, creates a new FSMPOSTPROC or raises the existing
%      singleton*.
%
%      H = FSMPOSTPROC returns the handle to a new FSMPOSTPROC or the handle to
%      the existing singleton*.
%
%      FSMPOSTPROC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMPOSTPROC.M with the given input arguments.
%
%      FSMPOSTPROC('Property','Value',...) creates a new FSMPOSTPROC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmPostProc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmPostProc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmPostProc

% Last Modified by GUIDE v2.5 01-Sep-2004 14:28:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmPostProc_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmPostProc_OutputFcn, ...
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

% --- Executes just before fsmPostProc is made visible.
function fsmPostProc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmPostProc (see VARARGIN)

% Choose default command line output for fsmPostProc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Check whether fsmCenter is running - if not, open it
hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if isempty(hfsmC)
    hfsmC=fsmCenter;
end
% Get current project
handlesFsmCenter=guidata(hfsmC);
projDir=get(handlesFsmCenter.fsmCenter,'UserData');
set(handles.textCurrentProject,'String',projDir);


% UIWAIT makes fsmPostProc wait for user response (see UIRESUME)
% uiwait(handles.fsmPostProc);

% --- Outputs from this function are returned to the command line.
function varargout = fsmPostProc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CREATEFCN
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editScale_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editBitDepth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editFileName_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editTSampling_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editLabel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSpectrumSigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editFrame_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function editD0_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editImageNumber_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editVectorAnalysis_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function EditTimeSigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editTimeSigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editYSP_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editXSP_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editFrameSP_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSMCL_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSMSampling_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSMPixel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSMMax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editCalX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editCalX0_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editCalY0_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editCalY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editFrameTN_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSigmaTN_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editSMBitDepth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  TURNOVER AND CYCLE ANALYSIS
% %
% %    (and related functions/callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function checkLoadBwMask_Callback(hObject, eventdata, handles)

function checkSaveBwMask_Callback(hObject, eventdata, handles)

function checkMeanScores_Callback(hObject, eventdata, handles)

function checkMoveROI_Callback(hObject, eventdata, handles)

function editTimeSigma_Callback(hObject, eventdata, handles)

% Grid size
function editY_Callback(hObject, eventdata, handles)
function editX_Callback(hObject, eventdata, handles)

% Sampling interval
function editTSampling_Callback(hObject, eventdata, handles)
set(handles.editTSampling,'String',num2str(abs(str2num(get(handles.editTSampling,'String')))));

% Label peaks with power > user setting
function editLabel_Callback(hObject, eventdata, handles)
if str2num(get(handles.editLabel,'String'))<0 | str2num(get(handles.editLabel,'String'))>1
    uiwait(msgbox('Please enter a value between 0 and 1.','Help','modal'));
    set(handles.editLabel,'String','0.3')
end

% Polygon type
function radioGrid_Callback(hObject, eventdata, handles)
set(handles.radioGrid,'Value',1);
set(handles.radioPoly,'Value',0);
set(handles.radioRect,'Value',0);
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
set(handles.textGrid,'Enable','off');
set(handles.textY,'Enable','off');
set(handles.textX,'Enable','off');
set(handles.editY,'Enable','off');
set(handles.editX,'Enable','off');
set(handles.textBwMask,'Enable','off');
set(handles.checkLoadBwMask,'Enable','off');
set(handles.checkSaveBwMask,'Enable','off');

function editSpectrumSigma_Callback(hObject, eventdata, handles)
% Make sure the entered value is > 0
set(handles.editSpectrumSigma,'String',num2str(abs(fix(str2num(get(handles.editSpectrumSigma,'String'))))));
if ~mod(str2num(get(handles.editSpectrumSigma,'String')),2)
    uiwait(msgbox('Please enter an integer, ODD value > 0.','Help','modal'));
    set(handles.editSpectrumSigma,'String','1');
end

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

function pushHelpKin_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_default.html')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  FLOW DISPLAY AND ANALYSIS
% %
% %    (and related functions/callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushVector_Callback(hObject, eventdata, handles)
projDir=getProjDir(handles);
if isempty(projDir)
    return
end
nAvg=str2num(get(handles.editVectorAnalysis,'String'));
displ=[get(handles.checkRaw,'Value') get(handles.checkInterp,'Value') get(handles.checkNoise,'Value') get(handles.checkError,'Value') get(handles.checkDisplayImg,'Value')];
roi=[get(handles.checkROI,'Value') get(handles.checkLoadROI,'Value') get(handles.checkSaveROI,'Value')];
scale=str2num(get(handles.editScale,'String'));
d0=str2num(get(handles.editD0,'String'));
useDiv=get(handles.checkDiv,'Value');
output=find([get(handles.radioCMap,'Value') get(handles.radioCircle,'Value')]);
displROI=get(handles.checkDisplayROI,'Value');
invertROI=get(handles.checkInvertROI,'Value');
fsmVectorAnalysis(projDir,nAvg,roi,displ,scale,d0,useDiv,output,displROI,invertROI);

function checkInvertROI_Callback(hObject, eventdata, handles)

function checkDisplayROI_Callback(hObject, eventdata, handles)

function checkRaw_Callback(hObject, eventdata, handles)

function checkInterp_Callback(hObject, eventdata, handles)

function checkNoise_Callback(hObject, eventdata, handles)

function checkDisplayImg_Callback(hObject, eventdata, handles)

function checkDiv_Callback(hObject, eventdata, handles)

function checkSMOverlay_Callback(hObject, eventdata, handles)

function checkSMMask_Callback(hObject, eventdata, handles)

function editScale_Callback(hObject, eventdata, handles)

function editD0_Callback(hObject, eventdata, handles)

function editVectorAnalysis_Callback(hObject, eventdata, handles)
nAvg=fix(str2num(get(handles.editVectorAnalysis,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
elseif nAvg>0
    if mod(nAvg,2)==0
        warndlg('The number of frames must be ODD.','Warning','modal')
        nAvg=1;
    end
end
set(handles.editVectorAnalysis,'String',num2str(nAvg));

function checkROI_Callback(hObject, eventdata, handles)
set(handles.checkLoadROI,'Value',0);

function checkSaveROI_Callback(hObject, eventdata, handles)
if get(handles.checkLoadROI,'Value')==1 | get(handles.checkROI,'Value')==0
    set(handles.checkSaveROI,'Value',0);
end

function checkLoadROI_Callback(hObject, eventdata, handles)
set(handles.checkROI,'Value',0);
set(handles.checkSaveROI,'Value',0);

function radioCMap_Callback(hObject, eventdata, handles)
set(handles.radioCMap,'Value',1);
set(handles.radioCircle,'Value',0);

function radioCircle_Callback(hObject, eventdata, handles)
set(handles.radioCMap,'Value',0);
set(handles.radioCircle,'Value',1);

function checkError_Callback(hObject, eventdata, handles)
if get(handles.checkError,'Value')==1
    set(handles.radioCircle,'Enable','on');
    set(handles.radioCMap,'Enable','on');
else
    set(handles.radioCircle,'Enable','off');
    set(handles.radioCMap,'Enable','off');
end

function pushHelpFlow_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_default.html')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  SPEED MAPS
% %
% %    (and related functions/callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushSpeedMaps_Callback(hObject, eventdata, handles)
if get(handles.checkSMMask,'Value')==1
    % Load userROI.mat
    [fName,dirName] = uigetfile(...
        {'userROI.mat;','userROI.mat';
        '*.*','All Files (*.*)'},...
        'Select userROI.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        choice=questdlg('You didn''t pick any file.','Error','Continue without ROI','Exit','Exit');
        switch choice,
            case 'Continue without ROI'; userROIpoly=[]; filePicked=0;
            case 'Exit', return;
        end % switch
    else
        filePicked=1;
    end
    if filePicked==1
        % Load the file
        load([dirName,fName]);
    end
    if filePicked==1 & exist('userROIpoly')~=1
        choice=questdlg('The loaded file does not seem to contain a valid ROI.','Error','Continue without ROI','Exit','Exit');
        switch choice,
            case 'Continue without ROI', userROIpoly=[];
            case 'Exit', return;
        end % switch
    end
else
    userROIpoly=[];
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
outputdir=fsmSpeedMaps(gridSize,n,d0,loadMPM,sampling,pixelSize,overlayVectors,userROIpoly,maxSpeed,segment,bitDepth);
if ~isempty(outputdir)
    % Maps have been saved to disk
    msg=['Speed maps have been saved to ',outputdir,'.'];
    uiwait(msgbox(msg,'Help','modal'));
end

function pushHelpSpeed_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_default.html')]);

function editFrameSP_Callback(hObject, eventdata, handles)
nAvg=fix(str2num(get(handles.editFrameSP,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
end
set(handles.editFrameSP,'String',num2str(nAvg));

function radioSMMPM_Callback(hObject, eventdata, handles)
set(handles.radioSMMPM,'Value',1);
set(handles.radioSMMD,'Value',0);

function radioSMMD_Callback(hObject, eventdata, handles)
set(handles.radioSMMPM,'Value',0);
set(handles.radioSMMD,'Value',1);

function checkSMSegment_Callback(hObject, eventdata, handles)
if get(handles.checkSMSegment,'Value')==1
    set(handles.textSMBitDepth,'Enable','on');
    set(handles.editSMBitDepth,'Enable','on');
else
    set(handles.textSMBitDepth,'Enable','off');
    set(handles.editSMBitDepth,'Enable','off');
end    

function editYSP_Callback(hObject, eventdata, handles)

function editXSP_Callback(hObject, eventdata, handles)

function editSMCL_Callback(hObject, eventdata, handles)

function editSMBitDepth_Callback(hObject, eventdata, handles)

function editSMSampling_Callback(hObject, eventdata, handles)

function editSMPixel_Callback(hObject, eventdata, handles)

function editSMMax_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  TURNOVER MAPS
% %
% %     (and related callbacks)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushHelpTurnoverMaps_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_turnoverMaps.html')]);

function pushKinMaps_Callback(hObject, eventdata, handles)
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

function pushHelpTurnover_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_default.html')]);

function editFrameTN_Callback(hObject, eventdata, handles)
nAvg=fix(str2num(get(handles.editFrameTN,'String')));
if isempty(nAvg)
    nAvg=1;
end
if nAvg<0
    nAvg=1;
end
set(handles.editFrameTN,'String',num2str(nAvg));

function editSigmaTN_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  PLOT TRAJECTORIES
% %
% %    (and related functions/callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushPlotTraj_Callback(hObject, eventdata, handles)
method=find([get(handles.radioPTMethod1,'Value') get(handles.radioPTMethod2,'Value')]);
colorCode=get(handles.checkTrajPlots,'Value');
fsmPlotTrajectories(method,colorCode);

function pushHelpPlotTraj_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_plotTrajectories.html')]);

function radioPTMethod1_Callback(hObject, eventdata, handles)
set(handles.radioPTMethod1,'Value',1);
set(handles.radioPTMethod2,'Value',0);

function radioPTMethod2_Callback(hObject, eventdata, handles)
set(handles.radioPTMethod1,'Value',0);
set(handles.radioPTMethod2,'Value',1);

function checkTrajPlots_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ASSIGN SPECKLES TO CLASSES
% %
% %     (and related callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushAssignSpecklesToClasses_Callback(hObject, eventdata, handles)
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

function pushHelpAssignClasses_Callback(hObject, eventdata, handles)
web(['file:///' which('qFSM_speckleClasses.html')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  EXIT CALLBACK
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CloseRequestFcn
function fsmPostProc_CloseRequestFcn(hObject, eventdata, handles)
fsmH=findall(0,'Tag','fsmPostProc'); % Get the handle of fsmPostProc
choice=questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');
switch choice,
    case 'Yes', delete(fsmH);
    case 'No', return;
end


% Menu
function menuHelp_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ACCESSORY FUNCTIONS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function projDir=getProjDir(handles)
% Chech whether a project exists
projDir=get(handles.textCurrentProject,'String');
if isempty(projDir)
    warndlg('No project defined. Please create/load a project in fsmCenter.','Warning','modal');
    return
else
    % Check that the project directory exists
    if ~isdir(projDir)
        warndlg('The directory spedified does not exist. Please check your project.','Warning','modal');
        return
    end
end
% TODO: check for sub-directories 'tack'
projDir=[projDir,filesep,'tack',filesep];
disp('Check for sub-directories ''tack''');
