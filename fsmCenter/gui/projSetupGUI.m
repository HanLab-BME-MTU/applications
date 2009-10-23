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
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before projSetupGUI is made visible.
function projSetupGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for projSetupGUI
%handles.output = hObject;

if nargin > 7
   error('Wrong number of input arguments.');
end

%First check if it is openned from 'fsmCenter'. If yes, always open it.
% Otherwise, check wether 'fsmCenter' if running. If yes, do not open.
callObjTag  = '';
callObjName = '';
if nargin == 7
   callObj = varargin{4};

   if ishandle(callObj)
      callObjTag  = get(callObj,'Tag');
      callObjName = get(callObj,'Name');
   end
end

if ~strcmp(callObjTag,'fsmCenter') || ~strcmp(callObjName,'fsmCenter')
   %It is not openned from 'fsmCenter'. Check whether 'fsmCenter' is running.
   hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');

   if ~isempty(hfsmC) && ishandle(hfsmC)
      handles.projDir     = '';
      handles.imgDirList  = {};
      handles.subProjDir  = {};
      guidata(hObject,handles);

      %'fsmCenter' is running.
      warndlg(['fsmCenter is running. In this case, you are ' ...
         'only allowed to change the project structure there.']);

      projSetupGUI('closeRequestFcn',hObject,[],guidata(hObject));
      return;
   end
end

subProjTags = {'tack', ...
               'lpla', ...
               'post', ...
               'edge', ...
               'merg', ...
               'fadh', ...
               'corr', ...
               'mech'};

subProjTitle = {'   speckTackle', ...
                'fsm Transition', ...
                '    fsm Center', ...
                '  Edge Tracker', ...
                '        Merger', ...
                'Focal Adhesion', ...
                '   corrTracker', ...
                '   Cont. Mech.'};

numSubProj = length(subProjTags);

unixMntRootMenu = {'/mnt/', ...
                   '/Volumes/', ...
                   '/'};

%The default unix mount root is '/mnt'.
handles.unixMntRootMI = 1;

%Remove space in 'subProjTitle'.
subProjNames = cell(numel(subProjTitle), 1);
for k = 1:numSubProj
   subProjNames{k} = sscanf(subProjTitle{k},'%s');
end

%Get handles to GUI objects.
handles.platformTH    = findobj(hObject,'tag','platformText');
handles.unixMntRootMH = findobj(hObject,'tag','unixMntRootMenu');
handles.projDirEH     = findobj(hObject,'tag','projDirEdit');
handles.imgDirMH      = findobj(hObject,'tag','imgDirMenu');
handles.imgDirEH      = findobj(hObject,'tag','imgDirEdit');
handles.firstImgTH    = findobj(hObject,'tag','firstImgText');

for k = 1:numSubProj
   handles.subProjMH{k}  = findobj(hObject,'tag',[subProjTags{k} 'Suffix']);
   handles.subProjTFH{k} = findobj(hObject,'tag',[subProjTags{k} 'SufNew']);
end

handles.unixMntRootMenu = unixMntRootMenu;
handles.subProjTags     = subProjTags;
handles.numSubProj      = numSubProj;
handles.subProjTitle    = subProjTitle;
handles.subProjNames    = subProjNames;

%Set the default unix mount root menu.
set(handles.unixMntRootMH,'String',handles.unixMntRootMenu);
set(handles.unixMntRootMH,'Value',handles.unixMntRootMI);
if isunix
   set(handles.unixMntRootMH,'Enable','on');
   set(handles.platformTH,'String','Unix Based');
else
   set(handles.unixMntRootMH,'Enable','off');
   set(handles.platformTH,'String','Windows PC Based');
end

projDir = '';
if nargin > 5
   projDir = varargin{3};
end

if ~isdir(projDir)
   projDir = pwd;
end

handles.selImgDir = 1;

handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
handles = updateGUI(handles);

handles.figH = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projSetupGUI wait for user response (see UIRESUME)
uiwait(handles.figH);

function updateDirList(dirMH,dirTFH,dirList,selDir)

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

function handles = updateGUI(handles)

numSubProj  = handles.numSubProj;
subProjTags = handles.subProjTags;
subProjMH   = handles.subProjMH;
subProjTFH  = handles.subProjTFH;
%platformTH  = handles.platformTH;

projDir       = handles.projDir;
imgDirList    = handles.imgDirList;
firstImgList  = handles.firstImgList;
%unix_imgDrive = handles.unix_imgDrive;
%win_imgDrive  = handles.win_imgDrive;
subProjDir    = handles.subProjDir;
selImgDir     = handles.selImgDir;

if isempty(firstImgList)
   firstImg = '';
else
   firstImg = firstImgList{selImgDir};
end

%Update the menu for unix mount root.
set(handles.unixMntRootMH,'Value',handles.unixMntRootMI);

%Update the menu for image directory list.
if selImgDir == 0
   imgDir = '';
   set(handles.imgDirMH,'String',{'--- Empty ---'});
   set(handles.imgDirMH,'Value',1);
   set(handles.imgDirMH,'Enable','off');
else
   imgDir = imgDirList{selImgDir};
   set(handles.imgDirMH,'String',imgDirList);
   set(handles.imgDirMH,'Value',selImgDir);
   set(handles.imgDirMH,'Enable','on');
end

set(handles.projDirEH,'string',projDir);
set(handles.firstImgTH,'string',firstImg);
set(handles.imgDirEH,'string',imgDir);

%Search for available results directories for each package.
projDirStruct=dir(projDir);
subProjDirList = cell(size(subProjDir));
for k = 1:numSubProj
   subProjDirList{k} = findProjSubDir(projDirStruct,subProjTags{k});
   updateDirList(subProjMH{k},subProjTFH{k},subProjDirList{k}, ...
       subProjDir{k});
end


% --- Outputs from this function are returned to the command line.
function varargout = projSetupGUI_OutputFcn(hObject, eventdata, handles) %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargout > 6
    error('Too many output arguments.');
end

if isempty(handles)
    for k = 1:nargout
        varargout{k} = '';
    end
    return;
end

unixMntRoot  = handles.unixMntRootMenu(handles.unixMntRootMI);
projDir      = handles.projDir;
imgDirList   = handles.imgDirList;
selImgDir    = handles.selImgDir;
firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjDir   = handles.subProjDir;
subProjTitle = handles.subProjTitle;
physiParam   = handles.physiParam;

unix_imgDrive     = handles.unix_imgDrive;
win_imgDrive      = handles.win_imgDrive;
lastUnix_imgDrive = handles.lastUnix_imgDrive;
lastWin_imgDrive  = handles.lastWin_imgDrive;

%Write image path to a file named 'lastProjSetting.mat'.
if isdir(projDir)
   if ~isempty(imgDirList)
      %Switch so that the first img dir is always the selected one.
      imgDir = imgDirList{selImgDir};
      imgDirList{selImgDir} = imgDirList{1};
      imgDirList{1}         = imgDir;
   
      if length(unix_imgDrive) == length(imgDirList)
         imgDrive = unix_imgDrive{selImgDir};
         unix_imgDrive{selImgDir} = unix_imgDrive{1};
         unix_imgDrive{1}         = imgDrive;
      end

      if length(win_imgDrive) == length(imgDirList)
         imgDrive = win_imgDrive{selImgDir};
         win_imgDrive{selImgDir} = win_imgDrive{1};
         win_imgDrive{1}         = imgDrive;
      end

      if length(lastUnix_imgDrive) == length(imgDirList)
         imgDrive = lastUnix_imgDrive{selImgDir};
         lastUnix_imgDrive{selImgDir} = lastUnix_imgDrive{1};
         lastUnix_imgDrive{1}         = imgDrive;
      end

      if length(lastWin_imgDrive) == length(imgDirList)
         imgDrive = lastWin_imgDrive{selImgDir};
         lastWin_imgDrive{selImgDir} = lastWin_imgDrive{1};
         lastWin_imgDrive{1}         = imgDrive;
      end

      if ~isempty(firstImgList)
         %Switch so that the first img dir is always the selected one.
         firstImg = firstImgList{selImgDir};
         firstImgList{selImgDir} = firstImgList{1};
         firstImgList{1}         = firstImg;
      end
   end

   projSettings.unixMntRoot   = unixMntRoot;
   projSettings.unix_imgDrive = unix_imgDrive;
   projSettings.win_imgDrive  = win_imgDrive;

   if ~isempty(physiParam)
      if length(physiParam) < length(firstImgList)
         defPhysiParam = physiParam{end};
         for k = length(physiParam)+1:length(firstImgList)
            physiParam{k} = defPhysiParam;
         end
      elseif length(physiParam) > length(firstImgList)
         physiParam(length(firstImgList)+1:length(physiParam)) = [];
      end
      selPhysiParam = physiParam{selImgDir};
      physiParam{selImgDir} = physiParam{1};
      physiParam{1}         = selPhysiParam;
   end


   projSettings.projDir    = projDir;
   projSettings.subProjDir = subProjDir;
   projSettings.physiParam = physiParam;

   settingsMatFile = [projDir filesep 'lastProjSettings.mat'];
   if isunix==1
      projSettings.unix_imgDirList = imgDirList;

      if ~iscell(lastUnix_imgDrive)
         samUnixImgDrive = 0;
      elseif length(lastUnix_imgDrive) ~= length(unix_imgDrive)
         samUnixImgDrive = 0;
      else
         samUnixImgDrive = 1;
      end

      k = 1;
      while k <= length(unix_imgDrive) && samUnixImgDrive
         samUnixImgDrive = samdir(lastUnix_imgDrive{k},unix_imgDrive{k});
         k = k+1;
      end

      if samUnixImgDrive
         if length(lastWin_imgDrive) == length(imgDirList)
            projSettings.win_imgDirList = ...
               dirUnix2PC(imgDirList,lastWin_imgDrive,unix_imgDrive);
            projSettings.win_imgDrive = lastWin_imgDrive;
         else
            projSettings.win_imgDrive   = {};
            projSettings.win_imgDirList = {};
         end
      else
         projSettings.win_imgDrive   = {};
         projSettings.win_imgDirList = {};
      end
   elseif ispc==1
      projSettings.win_imgDirList = imgDirList;
      if ~iscell(lastWin_imgDrive)
         samWinImgDrive = 0;
      elseif length(lastWin_imgDrive) ~= length(win_imgDrive)
         samWinImgDrive = 0;
      else
         samWinImgDrive = 1;
      end

      k = 1;
      while k <= length(win_imgDrive) && samWinImgDrive
         samWinImgDrive = samdir(lastWin_imgDrive{k},win_imgDrive{k});
         k = k+1;
      end

      if samWinImgDrive
         if length(lastUnix_imgDrive) == length(imgDirList)
            projSettings.unix_imgDirList = ...
               dirPC2Unix(imgDirList,lastUnix_imgDrive);
            projSettings.unix_imgDrive = lastUnix_imgDrive;
         else
            projSettings.unix_imgDrive   = {};
            projSettings.unix_imgDirList = {};
         end
      else
         projSettings.unix_imgDrive   = {};
         projSettings.unix_imgDirList = {};
      end
   else
      error('Platform not supported.');
   end
   projSettings.firstImgList = firstImgList; 

   try
       save(settingsMatFile,'projSettings');
   catch ME
       warndlg(ME.message);
   end
end

if isempty(imgDirList)
   imgDir = '';
else
   imgDir = imgDirList{1};
end

if nargout >= 1
    varargout{1} = projDir;
end
if nargout >= 2
    varargout{2} = imgDir;
end
if nargout >= 3
    varargout{3} = subProjDir;
end
if nargout >= 4
    varargout{4} = imgDirList;
end
if nargout >= 5
    varargout{5} = firstImgList;
end
if nargout == 6
    varargout{6} = physiParam;
end

delete(hObject);

function closeRequestFcn(hObject,eventdata,handles) %#ok<INUSL>

handles.projDir      = '';
handles.imgDirList   = {};
handles.firstImgList = {};
handles.subProjDir   = {};

% Update handles structure
guidata(hObject, handles);

delete(hObject);


% --- Executes on button press on Cancel.
function cancel_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.projDir      = '';
handles.imgDirList   = {};
handles.firstImgList = {};
handles.subProjDir   = {};

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figH);


% --- Executes on button press on Ok.
function Ok_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to Ok(see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir    = get(handles.projDirEH,'string');
if isempty(projDir)
    warndlg('No project directory is set.','Warning','modal');
    return;
end

selImgDir  = handles.selImgDir;

if selImgDir == 0;
   warndlg('No image directory is set. ' ,'Warning','modal');
   return;
end

imgDirList = handles.imgDirList;

%firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjTags  = handles.subProjTags;
subProjDir   = cell(size(subProjTags));

for k = 1:numSubProj
   item = get(handles.subProjMH{k},'value');
   menu = get(handles.subProjMH{k},'string');
   if iscell(menu) && strcmp(menu{item},'New') == 0
      subProjDir{k} = menu{item};
   else
      subProjDir{k} = [subProjTags{k} get(handles.subProjTFH{k},'string')];
   end
end

subProjOK = ones(1,numSubProj);
for k = 1:numSubProj
   if ~isdir([projDir filesep subProjDir{k}])
      subProjOK(k) = mkdir(projDir,subProjDir{k});
   end
   if subProjOK(k) ~= 1
      fprintf(1,'%s\n',['Trouble making ' subProjTags{k} ' directory: ' ...
         subProjDir{k} '.']);
      subProjDir{k} = '';
   end
end

handles.projDir    = projDir;
handles.imgDirList = imgDirList;
handles.selImgDir  = selImgDir;
handles.subProjDir = subProjDir;

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figH);


% --- Executes on button press on Browse.
function projDirBrowse_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to projDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = uigetdir('','Please select your project directory');

if isnumeric(projDir) && projDir == 0
    return;
end

subProjNames = handles.subProjNames;
subProjTags  = handles.subProjTags;

if ~samdir(handles.projDir,projDir)
   handles.selImgDir = 1;
   handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
   handles = updateGUI(handles);
end

guidata(hObject,handles);


% --- Executes on return projDir text field.
function projDir_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to projDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = get(handles.projDirEH,'string');

if ~isdir(projDir)
   errordlg('The entered path is not a directory.','Input error','modal');
   return;
end

subProjNames = handles.subProjNames;
subProjTags  = handles.subProjTags;
if ~samdir(handles.projDir,projDir)
    handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
    handles = updateGUI(handles);
end

guidata(hObject,handles);



% --- Executes on button press on Browse.
function delImgDir_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to imgDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selImgDir         = handles.selImgDir;
imgDirList        = handles.imgDirList;
unix_imgDrive     = handles.unix_imgDrive;
win_imgDrive      = handles.win_imgDrive;
lastUnix_imgDrive = handles.lastUnix_imgDrive;
lastWin_imgDrive  = handles.lastWin_imgDrive;

if selImgDir == 0
   return;
end

imgDirList(selImgDir)        = [];
if isunix
   if length(win_imgDrive) == length(unix_imgDrive)
      win_imgDrive(selImgDir) = [];
   end
   
   if length(lastWin_imgDrive) == length(unix_imgDrive)
      lastWin_imgDrive(selImgDir) = [];
   end
   
   unix_imgDrive(selImgDir)     = [];
   lastUnix_imgDrive(selImgDir) = [];   
elseif ispc
   if length(unix_imgDrive) == length(win_imgDrive)
      unix_imgDrive(selImgDir) = [];
   end
   
   if length(lastUnix_imgDrive) == length(win_imgDrive)
      lastUnix_imgDrive(selImgDir) = [];
   end
   
   win_imgDrive(selImgDir)      = [];
   lastWin_imgDrive(selImgDir)  = [];
end

firstImgList = handles.firstImgList;
if ~isempty(firstImgList)
   firstImgList(selImgDir) = [];
end

if selImgDir > length(imgDirList)
   selImgDir = selImgDir-1;
end

%Remove the corresponding physical parameter
if selImgDir > 0 && selImgDir <= length(handles.physiParam)
   handles.physiParam(selImgDir) = [];
end

if selImgDir == 0
   %All image directories are deleted.
   set(handles.imgDirMH,'Enable','off');
   set(handles.imgDirMH,'String',{'--- Empty ---'});
   set(handles.imgDirMH,'Value',1);
   set(handles.imgDirEH,'String','');
   set(handles.firstImgTH,'String','');
else
   if isunix
      [unixMntRoot,handles.unixMntRootMI] = ...
         autoExtractUnixMntRoot(handles.imgDirList{selImgDir},handles.unixMntRootMenu); %#ok<ASGLU>
      set(handles.unixMntRootMH,'Value',handles.unixMntRootMI);
   end
   %Update the menu list for image directory list.
   set(handles.imgDirMH,'Enable','on');
   set(handles.imgDirMH,'String',imgDirList);
   set(handles.imgDirMH,'Value',selImgDir);

   %Update the edit field for active (selected) image directory and the text
   %field for first image file name.
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

handles.selImgDir         = selImgDir;
handles.imgDirList        = imgDirList;
handles.firstImgList      = firstImgList;
handles.unix_imgDrive     = unix_imgDrive;
handles.win_imgDrive      = win_imgDrive;
handles.lastUnix_imgDrive = lastUnix_imgDrive;
handles.lastWin_imgDrive  = lastWin_imgDrive;

guidata(hObject,handles);


function imgDirMenu_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to imgDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

firstImgList = handles.firstImgList;
imgDirList   = handles.imgDirList;

selImgDir  = get(handles.imgDirMH,'Value');

%Update 'unixMntRoot'.
[unixMntRoot,handles.unixMntRootMI] = ...
   autoExtractUnixMntRoot(imgDirList{selImgDir},handles.unixMntRootMenu); %#ok<ASGLU>
set(handles.unixMntRootMH,'Value',handles.unixMntRootMI);

%Update the gui field for active (selected) image directory.
set(handles.imgDirEH,'String',imgDirList{selImgDir});

if isempty(firstImgList)
   set(handles.firstImgTH,'String','');
else
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

handles.selImgDir = selImgDir;
guidata(hObject,handles);

function imgDirEdit_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

selImgDir = handles.selImgDir;
imgDirList = handles.imgDirList;

if isempty(imgDirList)
   set(handles.imgDirEH,'String','');
else
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
end



function firstImgBrowse_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

oldDir = pwd;

selImgDir    = handles.selImgDir;
imgDirList   = handles.imgDirList;
firstImgList = handles.firstImgList;

if ~isempty(imgDirList)
   if ~isdir(imgDirList{selImgDir})
      errordlg(['Something is wrong with the current selected image ' ...
         'directory. If it does not exist any more, please delete it ' ...
         'using the delete button.'],'Image selection error','modal');
      return;
   end

   cd(imgDirList{selImgDir});
end

[firstImgFileName pathName filterIndex] = uigetfile( ...
   {'*.tif;*.gif;*.jpg;*.png', ...
   'Image Files (*.tif,*.gif,*.jpg,*.png)';
   '*.tif','TIFF files (*.tif)';
   '*.gif','GIF files (*.gif)';
   '*.jpg','JPEG files (*.jpg)';
   '*.png','PNG files (*.png)'},'Select the first image');

cd(oldDir);
if filterIndex == 0
   return;
end

%Check whether the selected new image directory exists in 'imgDirList'.
k = 1; imgDirExist = 0;
while ~imgDirExist && k <= length(imgDirList)
   if samdir(pathName,imgDirList{k})
      imgDirExist = 1;
   else
      k = k+1;
   end
end

if imgDirExist
   %First, make the matching image directory the active one.
   selImgDir = k; 
else
   %This is a new image directory. Add it to 'imgDirList'.
   imgDirList{end+1} = pathName;
   selImgDir = length(imgDirList);

   %Auto-extract drive name or drive letter.
   if isunix
      [unixMntRoot,handles.unixMntRootMI] = ...
         autoExtractUnixMntRoot(pathName,handles.unixMntRootMenu);
      unix_imgDrive = autoExtractUnixDriveName(pathName,unixMntRoot);
      handles.unix_imgDrive{end+1} = unix_imgDrive;
      handles.lastUnix_imgDrive{end+1} = '';
      
      handles.win_imgDrive   = {};
      handles.win_imgDirList = {};
      
      set(handles.unixMntRootMH,'Value',handles.unixMntRootMI);
   elseif ispc
      win_imgDrive = autoExtractWinDriveLetter(pathName);
      handles.win_imgDrive{end+1}  = win_imgDrive;
      handles.lastWin_imgDrive{end+1}  = '';
      
      handles.unix_imgDrive  = {};
      handles.unix_imgDirList = {};
   end
end
firstImgList{selImgDir} = firstImgFileName;

%Update the text field for first image file name.
set(handles.firstImgTH,'String',firstImgFileName);

if selImgDir ~= 0
   %Update the menu for image directory list.
   set(handles.imgDirMH,'String',imgDirList);
   set(handles.imgDirMH,'Value',selImgDir);
   set(handles.imgDirMH,'Enable','on');

   %Update the text field for active (selected) image directory.
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
else
   %Update the text field for active (selected) image directory.
   set(handles.imgDirEH,'String','');
end

handles.selImgDir    = selImgDir;
handles.imgDirList   = imgDirList;
handles.firstImgList = firstImgList;

guidata(hObject,handles);



function handles = getProjSetting(handles,projDir,subProjNames,subProjTags)

numSubProj = length(subProjNames);

subProjDir = cell(numSubProj, 1);
for k = 1:numSubProj
    subProjDir{k} = '';
end
imgDirList   = {};
firstImgList = {};

unix_imgDrive = {};
win_imgDrive  = {};

unixMntRoot = handles.unixMntRootMenu{handles.unixMntRootMI}; % unixMntRoot = '/mnt'

%Read last project setting in the selected project path.
noProblem = 0;
if isdir(projDir)
    title = 'Platform dependent';
    settingsMatFile = [projDir filesep 'lastProjSettings.mat'];
    if exist(settingsMatFile,'file') == 2
        s = load(settingsMatFile);
        projSettings = s.projSettings;
        subProjDir   = projSettings.subProjDir;
        firstImgList = projSettings.firstImgList;

        %Get saved physical parameter from last project.
        if isfield(projSettings,'physiParam')
            physiParam = projSettings.physiParam;
        else
            physiParam = cell(length(firstImgList), 1);
            for k = 1:length(firstImgList)
                physiParam{k} = getDefFsmPhysiParam;
            end
        end
        handles.subProjDir = subProjDir;
        handles.physiParam = physiParam;

        %First, try to get saved name of mounted unix drive.
        if isfield(projSettings,'unix_imgDrive')
            unix_imgDrive = projSettings.unix_imgDrive;
        end

        %Get saved unix mount root from last project.
        if isfield(projSettings,'unixMntRoot')
            unixMntRoot = projSettings.unixMntRoot;
            matchMI = strmatch(unixMntRoot,handles.unixMntRootMenu,'exact');

            if isempty(matchMI)
                warndlg(['The unix mount root stored in the last ' ...
                    'project setting is not recogonized. Will try ''/''.'], ...
                    'Warning','modal');
                unixMntRoot = '/';
                handles.unixMntRootMI = strmatch(unixMntRoot,handles.unixMntRootMenu,'exact');
            else
                handles.unixMntRootMI = matchMI(1);
            end
        end

        %Get saved image directory list and extract drive letter (PC) or drive
        % name (unix)
        %First, Unix based platform. Since in unix, the mounting root of a
        % disk drive is quite flexible, we will only try to extract the drive
        % mount root and drive name according to conventional rule (e.g. /mnt,
        % /Volumes (mac)).
        if isfield(projSettings,'unix_imgDirList') && ...
                ~isempty(projSettings.unix_imgDirList)
            unix_imgDirList = projSettings.unix_imgDirList;
            if ~isfield(projSettings,'unixMntRoot')
                %When 'unix_imgDirList' is a field of 'projSettings', 'unixMntRoot'
                % should also be a field of 'projSettings. However, in old project
                % settings, we did not store 'unixMntRoot'. For backward
                % compatibility, we try to extract mount root by checking if one of
                % the item in 'unixMntRootMenu' exist as a head string 'unix_imgDirList'.
                [unixMntRoot,handles.unixMntRootMI] = ...
                    autoExtractUnixMntRoot(unix_imgDirList{1},handles.unixMntRootMenu); %#ok<ASGLU>
            end

            if length(unix_imgDrive) ~= length(unix_imgDirList)
                unix_imgDrive = cell(size(unix_imgDirList));
            end
            for k = 1:length(unix_imgDirList)
                if ~isdir(unix_imgDrive{k})
                    %Either because it is not saved last time or the saved one is
                    % auto-detected and is wrong. Try auto-detect again.
                    mntRoot = autoExtractUnixMntRoot(unix_imgDirList{k}, ...
                        handles.unixMntRootMenu);
                    unix_imgDrive{k} = autoExtractUnixDriveName(unix_imgDirList{k},mntRoot);
                end
            end
        else
            unix_imgDirList = {};
        end

        %Then, PC based platform.
        if isfield(projSettings,'win_imgDirList') && ...
                ~isempty(projSettings.win_imgDirList)
            win_imgDirList = projSettings.win_imgDirList;
            for k = 1:length(win_imgDirList)
                win_imgDrive{k} = autoExtractWinDriveLetter(win_imgDirList{k});
            end
        else
            win_imgDirList = {};
        end

        %Check platform and convert image directory if neccessary.
        if isunix == 1
            if ~isempty(unix_imgDirList)
                noProblem = 1;
            else
                noProblem = 0;
            end

            for k = 1:length(unix_imgDirList)
                noProblem = noProblem & isdir(unix_imgDirList{k});
            end
            
            if noProblem
                imgDirList = unix_imgDirList;
            elseif ~isempty(win_imgDrive)
                tryAgain = 'Yes';
                winImgDriveNameList = win_imgDrive{1};
                for k = 2:length(win_imgDrive)
                    winImgDriveNameList = [winImgDriveNameList ', ' win_imgDrive{k}];
                end
                prompt = sprintf(['Last project is set up or modified in Windows. \n' ...
                    'The multi-channel image drive letters are \n   ' ...
                    winImgDriveNameList  '.\n\n' ...
                    'Please enter the corresponding image drive names \n' ...
                    'under Unix based platform: \n' ...
                    '(Enter one name if they are the same. ' ...
                    'Otherwise, separate by comma.)']);

                pat = '(/(\w+/*)+),*\s*';
                while strcmp(tryAgain,'Yes')
                    answer = inputdlg(prompt,'title',1,{''});
                    if isempty(answer)
                        tryAgain  = 'No';
                        noProblem = 0;
                        break;
                    end
                    unix_imgDrive = regexp(answer{1},pat,'tokens');
                    unix_imgDrive = [unix_imgDrive{:}];

                    if length(unix_imgDrive) ~= 1 && ...
                            length(unix_imgDrive) ~= length(win_imgDirList)
                        noProblem = 0;
                    else
                        %Extract unix mount root from 'unix_imgDrive'.
                        [unixMntRoot,handles.unixMntRootMI] = ...
                            autoExtractUnixMntRoot(unix_imgDrive{1},handles.unixMntRootMenu);

                        %Convert image directories to Unix format.
                        imgDirList = dirPC2Unix(win_imgDirList,unix_imgDrive);

                        %Check conversion.
                        if isempty(imgDirList)
                            noProblem = 0;
                        else
                            noProblem = 1;
                        end

                        for k = 1:length(unix_imgDirList)
                            noProblem = noProblem & isdir(unix_imgDirList{k});
                        end
                    end

                    if ~noProblem
                        question = sprintf(['Invalid or insufficient unix drive names or \n' ...
                            'the saved project setting file is corrupted or \n' ...
                            'image directories have changed since the last ' ...
                            'project is saved.\n' ...
                            'Please check the cause before you answer the following question.\n\n' ...
                            'Do you want to try again? \n' ...
                            'If no, the old image directories will be removed. ' ...
                            'Be cautious!!!']);
                        tryAgain = questdlg(question,'Alert','Yes','No','Yes');
                    else
                        projSettings.unix_imgDirList = imgDirList;
                        if length(unix_imgDrive) == 1
                            imgDrive = unix_imgDrive{1};
                            unix_imgDrive = cell(size(imgDirList));
                            for k = 1:length(imgDirList)
                                unix_imgDrive{k} = imgDrive;
                            end
                        end
                        noProblem = 1;
                        tryAgain  = 'No';
                    end
                end
            else
                warndlg(['The project setting file is corrupted. ' ...
                    'Project needs to be reset.'], 'Warning','modal');
                return;
            end
        elseif ispc == 1
            if isempty(win_imgDirList)
                noProblem = 0;
            else
                noProblem = 1;
            end
            
            for k = 1:length(win_imgDirList)
                noProblem = noProblem & isdir(win_imgDirList{k});
            end

            if noProblem
                imgDirList = win_imgDirList;
            elseif ~isempty(unix_imgDirList)
                tryAgain = 'Yes';
                unixImgDriveNameList = unix_imgDrive{1};
                for k = 2:length(unix_imgDrive)
                    unixImgDriveNameList = [unixImgDriveNameList ', ' unix_imgDrive{k}];
                end
                prompt = {sprintf(['Last project is set up or modified in Unix based platform. \n' ...
                    'The auto-detected image drive names are \n   ' ...
                    unixImgDriveNameList  '.\n\n' ...
                    'If it is not correct, please enter the correct image drive names: ']),
                    sprintf(['Please also enter the corresponding image drive letters ' ...
                    'under PC platform: \n' ...
                    '(Enter one name if they are the same. ' ...
                    'Otherwise, separate by comma.)'])};

                unixPat = '(/(\w+/*)+),*\s*';
                winPat  = '(\w+:?),*\s*';
                while strcmp(tryAgain,'Yes')
                    answer = inputdlg(prompt,'title',2,{unixImgDriveNameList,''});
                    if isempty(answer)
                        tryAgain  = 'No';
                        noProblem = 0;
                        retrun;
                    end
                    unix_imgDrive = regexp(answer{1},unixPat,'tokens');
                    win_imgDrive  = regexp(answer{2},winPat,'tokens');
                    unix_imgDrive = [unix_imgDrive{:}];
                    win_imgDrive  = [win_imgDrive{:}];

                    if length(win_imgDrive) ~= 1 && ...
                            length(win_imgDrive) ~= length(unix_imgDirList)
                        noProblem = 0;
                    elseif length(unix_imgDrive) ~= 1 && ...
                            length(unix_imgDrive) ~= length(unix_imgDirList)
                        noProblem = 0;
                    else
                        for k = 1:length(win_imgDrive)
                            colonInd = findstr(':',win_imgDrive{k});
                            if isempty(colonInd)
                                win_imgDrive{k}(end+1) = ':';
                            end
                        end
                        %Convert image directories to PC format.
                        imgDirList = dirUnix2PC(unix_imgDirList,win_imgDrive,unix_imgDrive);

                        %Check conversion.
                        if isempty(imgDirList)
                            noProblem = 0;
                        else
                            noProblem = 1;
                        end

                        for k = 1:length(imgDirList)
                            noProblem = noProblem & isdir(imgDirList{k});
                        end
                    end

                    if ~noProblem
                        question = sprintf(['Invalid or insufficient PC drive letters or ' ...
                            'unix drive names or \n' ...
                            'the saved project setting file is corrupted or \n'...
                            'image directories have changed since the last project. \n' ...
                            'Please check the cause before answer the following question.\n\n' ...
                            'Do you want to try again? \n' ...
                            'If no, the old image directories will be removed. ' ...
                            'Be cautious!!!']);
                        tryAgain = questdlg(question,'Alert','Yes','No','Yes');
                    else
                        projSettings.win_imgDirList = imgDirList;
                        if length(unix_imgDrive) == 1
                            imgDrive = unix_imgDrive{1};
                            unix_imgDrive = cell(size(imgDirList));
                            for k = 1:length(imgDirList)
                                unix_imgDrive{k} = imgDrive;
                            end
                        end
                        if length(win_imgDrive) == 1
                            imgDrive = win_imgDrive{1};
                            win_imgDrive = cell(size(imgDirList));
                            for k = 1:length(imgDirList)
                                win_imgDrive{k} = imgDrive;
                            end
                        end
                        noProblem = 1;
                        tryAgain  = 'No';
                    end
                end
            else
                warndlg(['The project setting file is corrupted. ' ...
                    'The project needs to be reset.'],'Warning','modal');
            end
        else
            error('Platform not supported.');
        end
    else
        % lastProjSettings.mat doesn't exist.
        noProblem = 0;
    end
else
    projDir = '';
end

if noProblem == 1 && length(firstImgList) ~= length(imgDirList)
    if ~isempty(firstImgList)
        noProblem = 0;
        warndlg(['Number of image directories and first images ' ...
            'do not match. Project setting file is likely corrupted ' ...
            'and therefore will be ignored.'],'Warning','modal');
    end
end

if ~noProblem
    imgDirList   = {};
    firstImgList = {};
    for k = 1:numSubProj
        subProjDir{k} = '';
    end

    win_imgDrive  = {};
    unix_imgDrive = {};
    physiParam    = [];
else
    %Check whether we can find all the subprojects in the loaded last
    %project settings. We do this check because there might be new
    %subprojects added.
    lastSubProjDir = subProjDir;
    for k = 1:numSubProj
        ind = strmatch(subProjTags{k},lastSubProjDir);
        if isempty(ind)
            subProjDir{k} = '';
        elseif length(ind) == 1
            subProjDir{k} = lastSubProjDir{ind};
        else
            error(['There are multiple directories for the same subproject. ' ...
                'The last project setting file is likely corruppted.']);
        end
    end
end

handles.projDir           = projDir;
handles.imgDirList        = imgDirList;
handles.firstImgList      = firstImgList;
handles.unix_imgDrive     = unix_imgDrive;
handles.win_imgDrive      = win_imgDrive;
handles.lastUnix_imgDrive = unix_imgDrive;
handles.lastWin_imgDrive  = win_imgDrive;
handles.subProjDir        = subProjDir;
handles.physiParam        = physiParam;

if isempty(imgDirList)
    handles.selImgDir = 0;
else
    handles.selImgDir = 1;
end

function [mntRoot,mntMI] = autoExtractUnixMntRoot(inDir,mntRootMenu)
%Extract the root of unix mounting point by checking through 'mntRootMenu'
%which stores possible mount root.

%First remove extra filesep in 'inDir'.
inDir = rmExtraFilesep(inDir);

for k = 1:length(mntRootMenu)
   if ~strcmp(mntRootMenu{k}(end),'/')
      mntRootMenu{k}(end+1) = '/';
   end
end

%First check whether '/' is in 'mntRootMenu'.
rootMI = strmatch('/',mntRootMenu,'exact');
if isempty(rootMI)
   rootMI = 0;
end

mntMI  = 1;
mntInd = [];
while mntMI ~= rootMI && mntMI < length(mntRootMenu) && isempty(mntInd)
   mntInd = findstr(mntRootMenu{mntMI}, inDir);
   mntMI  = mntMI+1;
end

if isempty(mntInd) || mntInd(1) ~= 1
   mntRoot = '/';
   mntMI = rootMI;
else
   mntMI = mntMI-1;
   mntRoot = mntRootMenu{mntMI};
end

function unix_drive = autoExtractUnixDriveName(inDir,mntRoot)

mntRoot = rmExtraFilesep(mntRoot);

if ~strcmp(mntRoot(end),'/')
   mntRoot = [mntRoot '/'];
end

if ~strcmp(mntRoot,'/')
   mntInd = findstr(mntRoot,inDir);
else
   mntInd = [];
end

inDir      = rmExtraFilesep(inDir);
filesepInd = findstr('/',inDir);
if isempty(mntInd)
   unix_drive = inDir(1:filesepInd(2)-1);
else
   unix_drive = inDir(1:filesepInd(3)-1);
end

function win_drive = autoExtractWinDriveLetter(inDir)

colonInd = findstr(':',inDir);
if ~isempty(colonInd)
   win_drive   = inDir(1:colonInd(1));
end
