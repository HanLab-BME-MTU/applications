function varargout = projSetupGUI(varargin)
% FLOWSAVEGUI M-file for projSetupGUI.fig
%      FLOWSAVEGUI, by itself, creates a new FLOWSAVEGUI or raises the existing
%      singleton*.
%
%      H = FLOWSAVEGUI returns the handle to a new FLOWSAVEGUI or the handle to
%      the existing singleton*.
%
%      FLOWSAVEGUI('Property','Value',...) creates a new FLOWSAVEGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to projSetupGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FLOWSAVEGUI('CALLBACK') and FLOWSAVEGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FLOWSAVEGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projSetupGUI

% Last Modified by GUIDE v2.5 13-Aug-2004 13:17:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @projSetupGUI_OpeningFcn, ...
    'gui_OutputFcn',  @projSetupGUI_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before projSetupGUI is made visible.
function projSetupGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for projSetupGUI
%handles.output = hObject;

handles.projDir = '';
handles.imgDir  = '';
handles.tackDir = '';
handles.lplaDir = '';
handles.postDir = '';
handles.edgeDir = '';
handles.mergDir = '';
handles.fadhDir = '';
handles.corrDir = '';
projDir=[];
imgDir=[];
resDirList={};
    
% if nargin==3
%     projDir=[];
%     resDirList={};
% end

if nargin > 5
    projDir = varargin{3};
end
if nargin > 6
    imgDir = varargin{4};
end
if nargin > 7 
    resDirList = varargin{5};
end

if isempty(projDir)
    %Use the current directory.
    projDir = pwd;
end

handles.projDir = projDir;
if ~isempty(resDirList)
    handles.tackDir = resDirList{1};
    handles.lplaDir = resDirList{2};
    handles.postDir = resDirList{3};
    handles.edgeDir = resDirList{4};
    handles.mergDir = resDirList{5};
    handles.fadhDir = resDirList{6};
    handles.corrDir = resDirList{7};
end


%Get handles to GUI objects.
handles.projDirTFH = findobj('tag','projDir');
handles.imgDirTFH  = findobj('tag','imgDir');
handles.tackDirMH  = findobj('tag','tackSuffix');
handles.tackDirTFH = findobj('tag','tackSufNew');
handles.lplaDirMH  = findobj('tag','lplaSuffix');
handles.lplaDirTFH = findobj('tag','lplaSufNew');
handles.postDirMH  = findobj('tag','postSuffix');
handles.postDirTFH = findobj('tag','postSufNew');
handles.edgeDirMH  = findobj('tag','edgeSuffix');
handles.edgeDirTFH = findobj('tag','edgeSufNew');
handles.mergDirMH  = findobj('tag','mergSuffix');
handles.mergDirTFH = findobj('tag','mergSufNew');
handles.fadhDirMH  = findobj('tag','fadhSuffix');
handles.fadhDirTFH = findobj('tag','fadhSufNew');
handles.corrDirMH  = findobj('tag','corrSuffix');
handles.corrDirTFH = findobj('tag','corrSufNew');

handles = getProjSetting(handles);

updateGUI(hObject,handles);

if ~isempty(imgDir)
    set(handles.imgDirTFH,'string',imgDir);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projSetupGUI wait for user response (see UIRESUME)
uiwait(hObject);

function updateDirList(dirMH,dirTFH,dirList,selDir,dirName)

if isempty(dirList)
    set(dirMH,'string',{'New'});
else
    set(dirMH,'string',{'New',dirList{:}});
end
set(dirMH,'value',1);
set(dirTFH,'string','');

if ~isempty(selDir) 
    tmpInd = find(strcmp(dirList,selDir));
    if ~isempty(tmpInd)
        set(dirMH,'value',tmpInd+1);
        %else
        %   if length(selDir) < 4 | strcmp(selDir(1:4),dirName) == 0
        %      error(['The first 4 characters of ' dirName ' directiory has to ' ...
        %         'be ''' dirName '''.']);
        %   else
        %      if length(selDir) > 4
        %         set(dirTFH,'string',selDir(5:end));
        %      end
        %   end
    end
end

function updateGUI(hObject,handles)

%Search for available results directories for each package and update popup
% menu in the GUI.
projDir = handles.projDir;
imgDir  = handles.imgDir;
tackDir = handles.tackDir;
lplaDir = handles.lplaDir;
postDir = handles.postDir;
edgeDir = handles.edgeDir;
mergDir = handles.mergDir;
fadhDir = handles.fadhDir;
corrDir = handles.corrDir;

%tackDirList = {};
%lplaDirList = {};
%postDirList = {};
%edgeDirList = {};
%mergDirList = {};
%corrDirList = {};

projDirStruct=dir(projDir);
tackDirList = findProjSubDir(projDirStruct,'tack');
lplaDirList = findProjSubDir(projDirStruct,'lpla');
postDirList = findProjSubDir(projDirStruct,'post');
edgeDirList = findProjSubDir(projDirStruct,'edge');
mergDirList = findProjSubDir(projDirStruct,'merg');
fadhDirList = findProjSubDir(projDirStruct,'fadh');
corrDirList = findProjSubDir(projDirStruct,'corr');

set(handles.projDirTFH,'string',projDir);
%if ~isempty(projDir)
%if isdir(projDir)
%dirList = dir(projDir);
%dirNameList = {dirList(find([dirList.isdir])).name};
%for k = 1:length(dirNameList)
%   if length(dirNameList{k}) >= 4
%      if strcmp(dirNameList{k}(1:4),'tack') == 1
%         tackDirList = {tackDirList{:} dirNameList{k}};
%      elseif strcmp(dirNameList{k}(1:4),'lpla') == 1
%         lplaDirList = {lplaDirList{:} dirNameList{k}};
%      elseif strcmp(dirNameList{k}(1:4),'post') == 1
%         postDirList = {postDirList{:} dirNameList{k}};
%      elseif strcmp(dirNameList{k}(1:4),'edge') == 1
%         edgeDirList = {edgeDirList{:} dirNameList{k}};
%      elseif strcmp(dirNameList{k}(1:4),'merg') == 1
%         mergDirList = {mergDirList{:} dirNameList{k}};
%      elseif strcmp(dirNameList{k}(1:4),'corr') == 1
%         corrDirList = {corrDirList{:} dirNameList{k}};
%      end
%   end
%end
%   end
%   set(handles.projDirTFH,'string',projDir);
%end

set(handles.imgDirTFH,'string',imgDir);

updateDirList(handles.tackDirMH,handles.tackDirTFH,tackDirList,tackDir,'tack');
updateDirList(handles.lplaDirMH,handles.lplaDirTFH,lplaDirList,lplaDir,'lpla');
updateDirList(handles.postDirMH,handles.postDirTFH,postDirList,postDir,'post');
updateDirList(handles.edgeDirMH,handles.edgeDirTFH,edgeDirList,edgeDir,'edge');
updateDirList(handles.mergDirMH,handles.mergDirTFH,mergDirList,mergDir,'merg');
updateDirList(handles.fadhDirMH,handles.fadhDirTFH,fadhDirList,fadhDir,'fadh');
updateDirList(handles.corrDirMH,handles.corrDirTFH,corrDirList,corrDir,'corr');



% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = projSetupGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargout > 3
    error('Too many output arguments.');
end
if isempty(handles)
    for k = 1:nargout
        varargout{k} = '';
    end
    return;
end

if nargout >= 1
    varargout{1} = handles.projDir;
end
if nargout >= 2
    varargout{2} = handles.imgDir;
end
if nargout == 3
    varargout{3} = {handles.tackDir,handles.lplaDir,handles.postDir, ...
            handles.edgeDir,handles.mergDir,handles.fadhDir,handles.corrDir};
end

%Write image path to a file named 'lastProjSetting.txt'.
if isdir(handles.projDir) 
    if isunix==1
        settingsFileName='lastProjSettings_unix.txt';
    elseif ispc==1
        settingsFileName='lastProjSettings_win.txt';
    else
        error('Platform not supported.');
    end
    fid = fopen([handles.projDir filesep settingsFileName],'w');
    if fid ~= -1
        fprintf(fid,'%s\n',['    Image Path: ' handles.imgDir]);
        fprintf(fid,'%s\n',['   speckTackle: ' handles.tackDir]);
        fprintf(fid,'%s\n',['fsm Transition: ' handles.lplaDir]);
        fprintf(fid,'%s\n',['    fsm Center: ' handles.postDir]);
        fprintf(fid,'%s\n',['  Edge Tracker: ' handles.edgeDir]);
        fprintf(fid,'%s\n',['        Merger: ' handles.mergDir]);
        fprintf(fid,'%s\n',['Focal Adhesion: ' handles.fadhDir]);
        fprintf(fid,'%s\n',['   corrTracker: ' handles.corrDir]);
        fclose(fid);
    end
end

delete(hObject);



% --- Executes on button press on Cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.projDir = '';
handles.imgDir  = '';
handles.tackDir = '';
handles.corrDir = '';

% Update handles structure
guidata(hObject, handles);

uiresume;


% --- Executes on button press on Ok.
function Ok_Callback(hObject, eventdata, handles)
% hObject    handle to Ok(see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = get(handles.projDirTFH,'string');
imgDir  = get(handles.imgDirTFH,'string');

if isempty(projDir)
    warnH = warndlg('No project directory is set.','Warning','modal');
    return;
end

if isempty(imgDir)
    warnH = warndlg('No image directory is set.','Warning','modal');
    return;
end

item = get(handles.tackDirMH,'value');
menu = get(handles.tackDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    tackDir = menu{item};
else
    tackDir = ['tack' get(handles.tackDirTFH,'string')];
end

item = get(handles.lplaDirMH,'value');
menu = get(handles.lplaDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    lplaDir = menu{item};
else
    lplaDir = ['lpla' get(handles.lplaDirTFH,'string')];
end

item = get(handles.postDirMH,'value');
menu = get(handles.postDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    postDir = menu{item};
else
    postDir = ['post' get(handles.postDirTFH,'string')];
end

item = get(handles.edgeDirMH,'value');
menu = get(handles.edgeDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    edgeDir = menu{item};
else
    edgeDir = ['edge' get(handles.edgeDirTFH,'string')];
end

item = get(handles.mergDirMH,'value');
menu = get(handles.mergDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    mergDir = menu{item};
else
    mergDir = ['merg' get(handles.mergDirTFH,'string')];
end

item = get(handles.fadhDirMH,'value');
menu = get(handles.fadhDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    fadhDir = menu{item};
else
    fadhDir = ['fadh' get(handles.fadhDirTFH,'string')];
end

item = get(handles.corrDirMH,'value');
menu = get(handles.corrDirMH,'string');
if iscell(menu) & strcmp(menu{item},'New') == 0
    corrDir = menu{item};
else
    corrDir = ['corr' get(handles.corrDirTFH,'string')];
end

tackOK = 1;
lplaOK = 1;
postOK = 1;
edgeOK = 1;
mergOK = 1;
fadhOK = 1;
corrOK = 1;
if ~isdir([projDir filesep tackDir])
    [tackOK,msg,msgID] = mkdir(projDir,tackDir);
end
if ~isdir([projDir filesep lplaDir])
    [lplaOK,msg,msgID] = mkdir(projDir,lplaDir);
end
if ~isdir([projDir filesep postDir])
    [postOK,msg,msgID] = mkdir(projDir,postDir);
end
if ~isdir([projDir filesep edgeDir])
    [edgeOK,msg,msgID] = mkdir(projDir,edgeDir);
end
if ~isdir([projDir filesep mergDir])
    [mergOK,msg,msgID] = mkdir(projDir,mergDir);
end
if ~isdir([projDir filesep fadhDir])
    [fadhOK,msg,msgID] = mkdir(projDir,fadhDir);
end
if ~isdir([projDir filesep corrDir])
    [corrOK,msg,msgID] = mkdir(projDir,corrDir);
end

if tackOK == 1 & lplaOK == 1 & postOK == 1 & edgeOK == 1 & ...
   mergOK == 1 & fadhOK == 1 & corrOK == 1
    handles.projDir = projDir;
    handles.imgDir  = imgDir;
    handles.tackDir = tackDir;
    handles.lplaDir = lplaDir;
    handles.postDir = postDir;
    handles.edgeDir = edgeDir;
    handles.mergDir = mergDir;
    handles.fadhDir = fadhDir;
    handles.corrDir = corrDir;
else
    warning('Trouble making new directory.');
    handles.projDir = '';
    handles.imgDir  = '';
    handles.tackDir = '';
    handles.lplaDir = '';
    handles.postDir = '';
    handles.edgeDir = '';
    handles.mergDir = '';
    handles.fadhDir = '';
    handles.corrDir = '';
end

% Update handles structure
guidata(hObject, handles);

uiresume;


% --- Executes on button press on Browse.
function projDirBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to projDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = uigetdir('','Please select your project directory');

if isnumeric(projDir) & projDir == 0
    return;
end

set(handles.projDirTFH,'string',projDir);
handles.projDir = projDir;

handles = getProjSetting(handles);

updateGUI(hObject,handles);

% --- Executes on return projDir text field.
function projDir_Callback(hObject, eventdata, handles)
% hObject    handle to projDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = get(handles.projDirTFH,'string');

if ~samdir(handles.projDir,projDir)
    handles.projDir = projDir;
    handles = getProjSetting(handles);
    updateGUI(hObject,handles);
end

% --- Executes on button press on Browse.
function imgDirBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to imgDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imgDir = uigetdir('','Please select your image directory');

if isnumeric(imgDir) & imgDir == 0
    return;
end

set(handles.imgDirTFH,'string',imgDir);
handles.imgDir = imgDir;



function handles = getProjSetting(handles)

projDir = handles.projDir;
tackDir = '';
lplaDir = '';
postDir = '';
edgeDir = '';
mergDir = '';
fadhDir = '';
corrDir = '';

%Read last project setting in the selected project path.
if isdir(projDir)
    if isunix==1
        settingsFileName='lastProjSettings_unix.txt';
    elseif ispc==1
        settingsFileName='lastProjSettings_win.txt';
    else
        error('Platform not supported.');
    end
    fid = fopen([projDir filesep settingsFileName],'r');
    if fid ~= -1
        noProblem = 1;
        textL = fgetl(fid);
        while noProblem & ischar(textL)
            k = 1;
            while k <= length(textL) & ~strcmp(textL(k),':')
                k = k+1;
            end
            if k == 1 | k > length(textL)
                noProblem = 0;
            else
                %Use 'sscanf' to remove space.
                headStr = sscanf(textL(1:k-1),'%s');
                switch headStr
                    case 'ImagePath'
                        if k == length(textL)
                            imgDir = '';
                        else
                            imgDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'speckTackle'
                        if k == length(textL)
                            tackDir = '';
                        else
                            tackDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'fsmTransition'
                        if k == length(textL)
                            lplaDir = '';
                        else
                            lplaDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'fsmCenter'
                        if k == length(textL)
                            postDir = '';
                        else
                            postDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'EdgeTracker'
                        if k == length(textL)
                            edgeDir = '';
                        else
                            edgeDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'Merger'
                        if k == length(textL)
                            mergDir = '';
                        else
                            mergDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'FocalAdhesion'
                        if k == length(textL)
                            fadhDir = '';
                        else
                            fadhDir = sscanf(textL(k+1:end),'%s');
                        end
                    case 'corrTracker'
                        if k == length(textL)
                            corrDir = '';
                        else
                            corrDir = sscanf(textL(k+1:end),'%s');
                        end
                    otherwise
                end
            end
            
            textL = fgetl(fid);
        end
        if ~noProblem
            warnH = warndlg(['Project setting file is corrupted ' ...
                    'and will be ignored.'],'Warning','modal');
            
            handles.imgDir  = '';
            handles.tackDir = '';
            handles.lplaDir = '';
            handles.postDir= '';
            handles.edgeDir= '';
            handles.mergDir= '';
            handles.fadhDir= '';
            handles.corrDir= '';
        else
            handles.imgDir  = imgDir;
            handles.tackDir = tackDir;
            handles.lplaDir = lplaDir;
            handles.postDir = postDir;
            handles.edgeDir = edgeDir;
            handles.mergDir = mergDir;
            handles.fadhDir = fadhDir;
            handles.corrDir = corrDir;
        end
        fclose(fid);
    end
end


