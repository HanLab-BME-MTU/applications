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
%      FSMCENTER('Property','Value',...) creates a new FSMCENTER or raises
%      the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmCenter_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to fsmCenter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmCenter

% Last Modified by GUIDE v2.5 10-Mar-2005 18:22:34

% Begin initialization code - DO NOT EDIT
% Check Matlab version: we restrict FSMCenter to Matlab 2006a or newer

matlabVersion=ver('MATLAB');
if str2double(matlabVersion.Version) < 7.2
    errordlg('FSM Center requires Matlab 2006a or newer.', 'FSM Center');
    return;
end

gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fsmCenter_OpeningFcn, ...
    'gui_OutputFcn',  @fsmCenter_OutputFcn, ...
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

% --- Executes just before fsmCenter is made visible.
function fsmCenter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmCenter (see VARARGIN)

% Choose default command line output for fsmCenter
handles.output = hObject;

settings = get(hObject,'UserData');

if isempty(settings) || ~isfield(settings,'projDir') || ...
   (isfield(settings,'projDir') && isempty(settings.projDir))
   handles.physiParamFile = 'fsmPhysiParam.mat';

   handles.physiParam = getDefFsmPhysiParam;
   handles = updateFsmGuiPhysiParamEdit(handles);

   %Disable all fsm software package and physical parameters entrance when there
   %is no project setup.
   handles = enableFsmPackages(handles, 'off');
   
   % Disable menus
   set(handles.menuCloseProject,'Enable', 'off');
   set(handles.menuProjectDescription, 'Enable', 'off');
   set(handles.menuProjectStatus, 'Enable', 'off');
end

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



function handles = updateFsmGuiPhysiParamEdit(handles)

physiParam = handles.physiParam;

set(handles.editNA,'String',num2str(physiParam.NA));
set(handles.editWaveLen,'String',num2str(physiParam.waveLen));
set(handles.editBitDepth,'String',num2str(physiParam.bitDepth));
set(handles.editPixelSize,'String',num2str(physiParam.pixelSize));
set(handles.editFrameInterval,'String',num2str(physiParam.frameInterval));
set(handles.editCaliSigma,'String',sprintf('%1.3f',physiParam.psfSigma));


function handles = enableFsmPackages(handles,fsmState)

if strcmp(fsmState,'off')
   set(handles.pushFSM,'Enable','off');
   set(handles.pushBatchJobs,'Enable','off');
   set(handles.pushImKymoAnalysis,'Enable','off');
   set(handles.pushEdgeTracker,'Enable','off');
   set(handles.pushFsmTransition,'Enable','off');
   set(handles.pushPostProc,'Enable','off');
else
   set(handles.pushFSM,'Enable','on');
   set(handles.pushBatchJobs,'Enable','on');
   set(handles.pushImKymoAnalysis,'Enable','on');
   set(handles.pushEdgeTracker,'Enable','on');
   set(handles.pushFsmTransition,'Enable','on');
   set(handles.pushPostProc,'Enable','on');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  CREATEFCN
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editField_CreateFcn(hObject, eventdata, handles)
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

function editCaliSigma_CreateFcn(hObject, eventdata, handles)
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

function editCalX0_CreateFcn(hObject, eventdata, handles)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  SPECKTACKLE
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushFSM_Callback(hObject, eventdata, handles)

hFsmGuiMain=findall(0,'Tag','fsmGuiMain','Name','SpeckTackle');
if isempty(hFsmGuiMain)
    fsmGuiMain;
elseif ishandle(hFsmGuiMain)
    figure(hFsmGuiMain);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  TOOLS
% %
% %     (and related callbacks)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushCrop_Callback(hObject, eventdata, handles)
y0S=get(handles.editCalY0,'String');
yS=get(handles.editCalY,'String');
x0S=get(handles.editCalX0,'String');
xS=get(handles.editCalX,'String');
if any([isempty(y0S) isempty(yS) isempty(x0S) isempty(xS)])
    area=[];
else
    y0=fix(str2double(y0S));
    y=fix(str2double(yS));
    x0=fix(str2double(x0S));
    x=fix(str2double(xS));
    if (y<=y0 || x<=x0 || y0<0 || x0<0) % We don't know the size of the image yet
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

function editCalY0_Callback(hObject, eventdata, handles)

function editCalY_Callback(hObject, eventdata, handles)

function editCalX0_Callback(hObject, eventdata, handles)

function editCalX_Callback(hObject, eventdata, handles)

function editCaliSigma_Callback(hObject, eventdata, handles)

function editFileName_Callback(hObject, eventdata, handles)

function pushUser_Callback(hObject, eventdata, handles)
fsmCenter_setUserSettings;

function pushBatchJobs_Callback(hObject, eventdata, handles)
fsmBatchJobs;

function pushImShow_Callback(hObject, eventdata, handles)
if get (handles.checkViewSequence,'Value')
    figure;fsmCenterCB_loadNewSequence;
else
    figure;fsmCenterCB_loadNewImage;
end

function checkViewSequence_Callback(hObject, eventdata, handles)

function pushEditExpParams_Callback(hObject, eventdata, handles)
userDir=fsmCenter_getUserSettings;
if ~isempty(userDir)
    edit([userDir,filesep,'fsmExpParams.txt']);
else
    edit('fsmExpParams.txt'); % Default fsmExpParams.txt file
end

function pushMeta_Callback(hObject, eventdata, handles)
stk2tif(get(handles.editFileName,'string'));

function pushConvertSA_Callback(hObject, eventdata, handles)
% Load speckleArray.mat
[fName,dirName] = uigetfile(...
    {'speckleArray.mat;','speckleArray.mat';
    '*.*','All Files (*.*)'},...
    'Select speckleArray.mat');
if ~(isa(fName,'char') && isa(dirName,'char'))
    return
end
load([dirName,fName]);
if exist('speckleArray', 'var')~=1
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
        if exist([path,filesep,'speckleArray.mat'], 'file') == 2
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

function pushCal_Callback(hObject, eventdata, handles)
userDir=fsmCenter_getUserSettings;
if ~isempty(userDir)
    dataFile=[userDir,filesep,'fsmExpParams.txt'];
else
    dataFile=which('fsmExpParams.txt'); % Default fsmExpParams.txt file
    % Check
    fsmMainPath=which('fsmMain.m');
    pathExpParams = getFilenameBody(dataFile);
    pathFsmMain = getFilenameBody(fsmMainPath);
    if ~strcmp(pathExpParams,pathFsmMain)
        uiwait(msgbox('Something is wrong with your configuration. Plase click on "User settings" in fsmCenter and follow the instructions. Then try again.','Error','modal'));
        return
    end
end
[I0,sDN,GaussRatio,status]=fsmCalcNoiseParam([],str2num(get(handles.editBitDepth,'string')),str2num(get(handles.editCaliSigma,'string')),dataFile);
if status==-1
    uiwait(msgbox('The noise parameters could not be saved to your user database. They were returned to the MATLAB console. Make sure that a file called ''fsmExpParams.txt'' exists in your user directory or in the root directory of the qFSM package. Click on ''User settings'' in fsmCenter and follow the instructions.','Info','modal'));
elseif status==0
    disp('Calibration interrupted by the user.');
elseif status==1
    uiwait(msgbox('The noise parameters were successfully added to your user database. They will be added to the ''Camera calibration'' popup menu the next time you''ll start SpeckTackle. If you need to modify the entered record, click on ''Edit experiment database'' in fsmCenter and scroll down to the end to the file.','Info','modal'));
else
    error('The function fsmCalcNoiseParam returned an invalid status.');
end

function pushOptimCal_Callback(hObject, eventdata, handles)
fsmCenterAlphaBetaGui;

function pushCalReset_Callback(hObject, eventdata, handles)
set(handles.editCalY0,'String','');
set(handles.editCalY,'String','');
set(handles.editCalX0,'String','');
set(handles.editCalX,'String','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  IMKYMOANALYSIS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushImKymoAnalysis_Callback(hObject, eventdata, handles)

hCorrTrack=findall(0,'Tag','imKymoAnalysis','Name','imKymoAnalysis');
if isempty(hCorrTrack)
    imKymoAnalysis;
elseif ishandle(hCorrTrack)
    figure(hCorrTrack);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  POST PROCESSING
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushPostProc_Callback(hObject, eventdata, handles)
fsmPostProc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  TRANSITION ANALYSIS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushFsmTransition_Callback(hObject, eventdata, handles)

fsmTransition;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  EDGE TRACKER
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushEdgeTracker_Callback(hObject, eventdata, handles)

prPanel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  EXIT
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fsmCenter_DeleteFcn(hObject, eventdata, handles)

function fsmCenter_CloseRequestFcn(hObject, eventdata, handles)

fsmH = findall(0,'Tag','fsmCenter'); % Get the handle of fsmCenter
choice = questdlg('Are you sure you want to exit?','Exit request','Yes','No','No');

if isempty(choice) || strcmp(choice, 'No') == 1
    return;
end

% fsmPostProc
hFsmPostProc=findall(0,'Tag','fsmPostProc','Name','SpeckTackle - Post processing');
if ~isempty(hFsmPostProc)
    fsmPostProc('fsmPostProc_CloseRequestFcn',fsmH,[],guidata(hFsmPostProc)); % The calling GUI is fsmH -> fsmCenter
end

% fsmGuiMain (SpeckTackle)
hFsmGuiMain=findall(0,'Tag','fsmGuiMain','Name','SpeckTackle');
if ~isempty(hFsmGuiMain)
    fsmGuiMain('fsmGuiMain_CloseRequestFcn',fsmH,[],guidata(hFsmGuiMain)); % The calling GUI is fsmH -> fsmCenter
end

% imKymoAnalysis
hCorrTrack=findall(0,'Tag','imKymoAnalysis','Name','imKymoAnalysis');
if ~isempty(hCorrTrack)
    imKymoAnalysis('closeRequest',fsmH,[],guidata(hCorrTrack)); % The calling GUI is fsmH -> fsmCenter
end

% fsmTransition
hFsmTransition=findall(0,'Tag','fsmTransition','Name','fsmTransition');
if ~isempty(hFsmTransition)
    fsmTransition('fsmTransition_CloseRequestFcn',fsmH,[],guidata(hFsmTransition)); % The calling GUI is fsmH -> fsmCenter
end

% fsmCenterAlphaBetaGui
hFsmCenterABGui=findall(0,'Tag','NoiseParameterOptimizer','Name','Optimizer');
if ~isempty(hFsmCenterABGui)
    fsmCenterAlphaBetaGui('NoiseParameterOptimizer_CloseRequestFcn',fsmH,[],guidata(hFsmCenterABGui)); % The calling GUI is fsmH -> fsmCenter
end

% And now close fsmCenter
delete(fsmH);

% Menu PROJECT

function menuProject_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  PROJECT MANAGEMENT
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuSetupProject_Callback(hObject, eventdata, handles)

projDir = get(handles.textCurrentProject, 'String');

[projDir, imageDir, subProjects, imgDirList, firstImgList, physiParam] = ...
   projSetupGUI('a','b',projDir, handles.fsmCenter);

if isempty(projDir)
    % Nothing to do here - the user just canceled
    return
end

% Update project info in fsmCenter
set(handles.textCurrentProject,'String',projDir);
set(handles.textCurrentImage,'String',imageDir);

% Add projDir, imageDir and subProjects to the userData of fsmCenter
settings.projDir = projDir;
settings.imageDir = {imageDir}; % This is a cell (n-channel movies)
settings.subProjects = subProjects;
settings.imgDirList = imgDirList;
settings.firstImgList = firstImgList;

set(handles.fsmCenter,'UserData', settings);

% Enable all fsm software package.
handles = enableFsmPackages(handles, 'on');

% Enable menus
set(handles.menuCloseProject,'Enable','on');
set(handles.menuProjectDescription, 'Enable', 'on');
set(handles.menuProjectStatus, 'Enable', 'on');

if ~isempty(physiParam)
    handles.physiParam = physiParam{1};
    saveFsmPhysiParam(handles);

    handles = updateFsmGuiPhysiParamEdit(handles);
end

guidata(handles.fsmCenter, handles);

% Update other GUIs if they are already running

% fsmPostProc
hFsmPostProc=findall(0,'Tag','fsmPostProc','Name','SpeckTackle - Post processing');
if ~isempty(hFsmPostProc)
    fsmPostProc;
end

% fsmGuiMain (SpeckTackle)
hFsmGuiMain=findall(0,'Tag','fsmGuiMain','Name','SpeckTackle');
if ~isempty(hFsmGuiMain)
    fsmGuiMain;
end

% imKymoAnalysis
hCorrTrack=findall(0,'Tag','imKymoAnalysis','Name','imKymoAnalysis');
if ~isempty(hCorrTrack)
    imKymoAnalysis;
end

% fsmTransition
hFsmTransition=findall(0,'Tag','fsmTransition','Name','fsmTransition');
if ~isempty(hFsmTransition)
    fsmTransition;
end

function menuCloseProject_Callback(hObject, eventdata, handles)
% TODO: is there anything to be saved before closing?

% Update project info in fsmCenter
set(handles.textCurrentProject, 'String', '');
set(handles.textCurrentImage, 'String', '');

% How to clear settings from UserData?
% set(handles.fsmCenter, 'UserData', settings);

handles.physiParam = getDefFsmPhysiParam;
handles = updateFsmGuiPhysiParamEdit(handles);

% Disable all fsm software package.
handles = enableFsmPackages(handles,'off');

% Disable menus
set(handles.menuCloseProject,'Enable','off');
set(handles.menuProjectDescription, 'Enable', 'off');
set(handles.menuProjectStatus, 'Enable', 'off');

guidata(handles.fsmCenter, handles);

function menuProjectDescription_Callback(hObject, eventdata, handles)
% This function allows the user to create a text file where all notes, comments, and description are stored.

projDir = get(handles.textCurrentProject,'String');

disp('test description')

% This becomes useless since the menu is disabled if no project is open.
% if isempty(projDir)
%     % No project active
%     uiwait(errordlg('Plese create/open a project first.','Error','modal'));
%     return;
% end

% Check whether a file description.txt exists in the project directory
if exist([projDir,filesep,'description.txt'],'file') == 2
    % Open it
    edit([projDir,filesep,'description.txt']);
else
    % Create it
    fid = fopen([projDir,filesep,'description.txt'],'w');
    if fid==-1
        uiwait(errordlg('Cannot write to the project directory.','Error','modal'));
        return;
    else
        % Add an initial text
        fwrite(fid,'Please use this file to store your comments, notes, descriptions for the project.');
        fclose(fid);
        edit([projDir,filesep,'description.txt']);
    end
end

function menuProjectStatus_Callback(hObject, eventdata, handles)

fsmCenterProjectStatus;

% Menu HELP
function menuHelp_Callback(hObject, eventdata, handles)

function menuFsmCenter_Callback(hObject, eventdata, handles)

function menuSpeckTackle_Callback(hObject, eventdata, handles)

function menuTools_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Call back functions for physical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editBitDepth_Callback(hObject, eventdata, handles)

bitDepth = str2double(get(hObject,'String'));
if isempty(bitDepth)
   set(hObject,'String',num2str(handles.physiParam.bitDepth));
   errordlg('Not valid numerical input.','Error','modal');
   return;
end

handles.physiParam.bitDepth = bitDepth;
guidata(hObject,handles);
saveFsmPhysiParam(handles);

function editNA_Callback(hObject, eventdata, handles)
NA = str2double(get(hObject,'String'));
if isempty(NA)
   set(hObject,'String',num2str(handles.physiParam.NA));
   errordlg('Not valid numerical input.','Error','modal');
   return;
end

physiParam = handles.physiParam;

physiParam.NA      = NA;
physiParam.psfSigma   = round(0.21*physiParam.waveLen/physiParam.NA ...
                     /physiParam.pixelSize*1000)/1000;
set(handles.editCaliSigma,'String',sprintf('%1.3f',physiParam.psfSigma));

handles.physiParam = physiParam;
guidata(hObject,handles);
saveFsmPhysiParam(handles);

function editWaveLen_Callback(hObject, eventdata, handles)
waveLen = str2double(get(hObject,'String'));
if isempty(waveLen)
   set(hObject,'String',num2str(handles.physiParam.waveLen));
   errordlg('Not valid numerical input.','Error','modal');
   return;
end

physiParam = handles.physiParam;

physiParam.waveLen = waveLen;
physiParam.psfSigma   = round(0.21*physiParam.waveLen/physiParam.NA ...
                     /physiParam.pixelSize*1000)/1000;
set(handles.editCaliSigma,'String',sprintf('%1.3f',physiParam.psfSigma));

handles.physiParam = physiParam;
guidata(hObject,handles);
saveFsmPhysiParam(handles);

function editPixelSize_Callback(hObject, eventdata, handles)
pixelSize = str2double(get(hObject,'String'));
if isempty(pixelSize)
   set(hObject,'String',num2str(handles.physiParam.pixelSize));
   errordlg('Not valid numerical input.','Error','modal');
   return;
end

physiParam = handles.physiParam;

physiParam.pixelSize = pixelSize;
physiParam.psfSigma     = round(0.21*physiParam.waveLen/physiParam.NA ...
                       /physiParam.pixelSize*1000)/1000;
set(handles.editCaliSigma,'String',sprintf('%1.3f',physiParam.psfSigma));

handles.physiParam = physiParam;
guidata(hObject,handles);
saveFsmPhysiParam(handles);

function editFrameInterval_Callback(hObject, eventdata, handles)
frameInterval = str2double(get(hObject,'String'));
if isempty(frameInterval)
   set(hObject,'String',num2str(handles.physiParam.frameInterval));
   errordlg('Not valid numerical input.','Error','modal');
   return;
end

handles.physiParam.frameInterval = frameInterval;
guidata(hObject,handles);
saveFsmPhysiParam(handles);


function saveFsmPhysiParam(handles)

settings = get(handles.fsmCenter,'UserData');
if isempty(settings)
   return;
elseif ~isfield(settings,'projDir')
   return;
elseif isempty(settings.projDir)
   return;
end

physiParamFile = [settings.projDir filesep handles.physiParamFile];
fsmPhysiParam  = handles.physiParam;

save(physiParamFile,'fsmPhysiParam');

settingsMatFile = [settings.projDir filesep 'lastProjSettings.mat'];
if exist(settingsMatFile,'file') == 2
   s = load(settingsMatFile);
   projSettings = s.projSettings;
   projSettings.physiParam{1} = handles.physiParam;

   save(settingsMatFile,'projSettings');
end
