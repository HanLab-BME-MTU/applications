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

subProjTags = {'tack', ...
               'lpla', ...
               'post', ...
               'edge', ...
               'merg', ...
               'fadh', ...
               'corr'};

subProjTitle = {'   speckTackle', ...
                'fsm Transition', ...
                '    fsm Center', ...
                '  Edge Tracker', ...
                '        Merger', ...
                'Focal Adhesion', ...
                '   corrTracker'};

numSubProj = length(subProjTags);

%Remove space in 'subProjTitle'.
subProjNames = cell(size(subProjTitle));
for k = 1:numSubProj
   subProjNames{k} = sscanf(subProjTitle{k},'%s');
end

projDir = '';
if nargin > 5
   projDir = varargin{3};
end

if ~isdir(projDir)
   projDir = pwd;
end

[imgDir,subProjDir] = getProjSetting(projDir,subProjNames);

%Get handles to GUI objects.
handles.projDirTFH = findobj('tag','projDir');
handles.imgDirTFH  = findobj('tag','imgDir');

%To get the handle to those 'subProj' GUI objects.
%Text Field Handle.
subProjTFH = cell(size(subProjTags));
%Menu Handle.
subProjMH = cell(size(subProjTags));
for k = 1:numSubProj
   handles.subProjMH{k}  = findobj('tag',[subProjTags{k} 'Suffix']);
   handles.subProjTFH{k} = findobj('tag',[subProjTags{k} 'SufNew']);
end

handles.subProjTags  = subProjTags;
handles.numSubProj   = numSubProj;
handles.subProjTitle = subProjTitle;
handles.subProjNames = subProjNames;


handles = updateGUI(handles,projDir,imgDir,subProjDir);

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
     else
        if length(selDir) > 4
           set(dirTFH,'String',selDir(5:end));
        end
    end
end

function handles = updateGUI(handles,projDir,imgDir,subProjDir)

numSubProj  = handles.numSubProj;
subProjTags = handles.subProjTags;
subProjMH   = handles.subProjMH;
subProjTFH  = handles.subProjTFH;

handles.projDir = projDir;
handles.imgDir  = imgDir;

%Search for available results directories for each package and update popup
% menu in the GUI.
set(handles.projDirTFH,'string',projDir);
set(handles.imgDirTFH,'string',imgDir);

projDirStruct=dir(projDir);
subProjDirList = cell(size(subProjDir));
for k = 1:numSubProj
   subProjDirList{k} = findProjSubDir(projDirStruct,subProjTags{k});
   updateDirList(subProjMH{k},subProjTFH{k},subProjDirList{k}, ...
      subProjDir{k},subProjTags{k});
end

handles.subProjDir = subProjDir;

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
    varargout{3} = handles.subProjDir;
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

       numSubProj   = handles.numSubProj;
       subProjDir   = handles.subProjDir;
       subProjTitle = handles.subProjTitle;
       for k = 1:numSubProj
          fprintf(fid,'%s\n',[subProjTitle{k} ': ' subProjDir{k}]);
       end
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

numSubProj  = handles.numSubProj;
subProjTags = handles.subProjTags;
subProjDir  = cell(size(subProjTags));

if isempty(projDir)
    warnH = warndlg('No project directory is set.','Warning','modal');
    return;
end

if isempty(imgDir)
    warnH = warndlg('No image directory is set.','Warning','modal');
    return;
end

for k = 1:numSubProj
   item = get(handles.subProjMH{k},'value');
   menu = get(handles.subProjMH{k},'string');
   if iscell(menu) & strcmp(menu{item},'New') == 0
      subProjDir{k} = menu{item};
   else
      subProjDir{k} = [subProjTags{k} get(handles.subProjTFH{k},'string')];
   end
end

subProjOK = 1; k =1;
while subProjOK & k <= numSubProj
   if ~isdir([projDir filesep subProjDir{k}])
      [subProjOK,msg,msgID] = mkdir(projDir,subProjDir{k});
   end
   k = k+1;
end

if subProjOK
   handles.projDir = projDir;
   handles.imgDir  = imgDir;
   handles.subProjDir = subProjDir;
else
   warning('Trouble making new directory.');
   handles.projDir = '';
   handles.imgDir  = '';
   for k = 1:numSubProj
      handles.subProjDir{k} = '';
   end
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

subProjNames = handles.subProjNames;
if ~samdir(handles.projDir,projDir)
    [imgDir subProjDir] = getProjSetting(projDir,subProjNames);
    handles = updateGUI(handles,projDir,imgDir,subProjDir);
end

guidata(hObject,handles);


% --- Executes on return projDir text field.
function projDir_Callback(hObject, eventdata, handles)
% hObject    handle to projDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = get(handles.projDirTFH,'string');

if ~isdir(projDir)
   errordlg('The entered path is not a directory.','Input error','modal');
   return;
end

subProjNames = handles.subProjNames;
if ~samdir(handles.projDir,projDir)
    [imgDir subProjDir] = getProjSetting(projDir,subProjNames);
    handles = updateGUI(handles,projDir,imgDir,subProjDir);
end

guidata(hObject,handles);



% --- Executes on button press on Browse.
function imgDirBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to imgDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isdir(handles.imgDir)
   imgDir = uigetdir([handles.imgDir filesep '..'], ...
      'Please select your image directory');
else
   imgDir = uigetdir('','Please select your image directory');
end

if isnumeric(imgDir) & imgDir == 0
    return;
end

set(handles.imgDirTFH,'string',imgDir);
handles.imgDir = imgDir;

guidata(hObject,handles);


function imgDir_Callback(hObject, eventdata, handles)
% hObject    handle to imgDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imgDir = get(handles.imgDirTFH,'string');

if ~isdir(imgDir)
   errordlg('The entered path is not a directory.','Input error','modal');
   return;
end

if ~samdir(handles.imgDir,imgDir)
    handles.imgDir = imgDir;
end

guidata(hObject,handles);


function [imgDir,subProjDir] = getProjSetting(projDir,subProjNames)

numSubProj = length(subProjNames);

subProjDir = cell(size(subProjNames));
for k = 1:numSubProj
   subProjDir{k} = '';
end
imgDir = '';

%Read last project setting in the selected project path.
noProblem = 0;
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
                if strcmp(headStr,'ImagePath')
                   if k ~= length(textL)
                      imgDir = sscanf(textL(k+1:end),'%s');
                   end
                else
                   j1 = 1;
                   while ~strcmp(headStr,subProjNames{j1}) & j1 <= numSubProj
                      j1 = j1+1;
                   end
                   if j1 <= numSubProj
                      if k ~= length(textL)
                         subProjDir{j1} = sscanf(textL(k+1:end),'%s');
                      end
                   end
                end
            end
            
            textL = fgetl(fid);
        end
        if ~noProblem
           warnH = warndlg(['Project setting file is corrupted ' ...
              'and will be ignored.'],'Warning','modal'); 
        end
        fclose(fid);
    end
end

if ~noProblem
   imgDir  = '';
   for k = 1:numSubProj
      subProjDir{k} = '';
   end
end

