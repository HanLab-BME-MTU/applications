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

if ~strcmp(callObjTag,'fsmCenter') | ~strcmp(callObjName,'fsmCenter')
   %It is not openned from 'fsmCenter'. Check whether 'fsmCenter' is running.
   hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');

   if ~isempty(hfsmC) & ishandle(hfsmC)
      handles.projDir     = '';
      handles.imgDirList  = {};
      handles.subProjDir  = {};
      guidata(hObject,handles);

      %'fsmCenter' is running.
      warnH = warndlg(['fsmCenter is running. In this case, you are ' ...
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

%Get handles to GUI objects.
handles.projDirEH  = findobj(hObject,'tag','projDirEdit');
handles.imgDirMH   = findobj(hObject,'tag','imgDirMenu');
handles.imgDirEH   = findobj(hObject,'tag','imgDirEdit');
handles.firstImgTH = findobj(hObject,'tag','firstImgText');

%To get the handle to those 'subProj' GUI objects.
%Text Field Handle.
subProjTFH = cell(size(subProjTags));
%Menu Handle.
subProjMH = cell(size(subProjTags));
for k = 1:numSubProj
   handles.subProjMH{k}  = findobj(hObject,'tag',[subProjTags{k} 'Suffix']);
   handles.subProjTFH{k} = findobj(hObject,'tag',[subProjTags{k} 'SufNew']);
end

handles.subProjTags  = subProjTags;
handles.numSubProj   = numSubProj;
handles.subProjTitle = subProjTitle;
handles.subProjNames = subProjNames;

projDir = '';
if nargin > 5
   projDir = varargin{3};
end

if ~isdir(projDir)
   projDir = pwd;
end

selImgDir = 1;
[imgDirList,firstImgList,subProjDir] = getProjSetting(projDir,subProjNames);
handles = updateGUI(handles,projDir,imgDirList,selImgDir, ...
   firstImgList,subProjDir);

handles.figH = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projSetupGUI wait for user response (see UIRESUME)
uiwait(handles.figH);

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

function handles = updateGUI(handles,projDir,imgDirList,selImgDir, ...
   firstImgList,subProjDir)

numSubProj  = handles.numSubProj;
subProjTags = handles.subProjTags;
subProjMH   = handles.subProjMH;
subProjTFH  = handles.subProjTFH;

handles.projDir      = projDir;
handles.imgDirList   = imgDirList;
handles.firstImgList = firstImgList;

handles.selImgDir = selImgDir;

if selImgDir > length(imgDirList) | isempty(firstImgList)
   firstImg = '';
else
   firstImg = firstImgList{selImgDir};
end

if selImgDir > length(imgDirList)
   imgDir = '';
else
   imgDir = imgDirList{selImgDir};
end

if isempty(imgDirList)
   imgDirList = {'-- New Image Directory (by selecting first image) --'};
else
   imgDirList{end+1} = '-- New Image Directory (by selecting first image) --';
end

%Update popup menu in the GUI.
set(handles.projDirEH,'string',projDir);
set(handles.firstImgTH,'string',firstImg);
set(handles.imgDirMH,'Value',1);
set(handles.imgDirMH,'String',imgDirList);
set(handles.imgDirMH,'Value',selImgDir);
set(handles.imgDirEH,'string',imgDir);

%Search for available results directories for each package.
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

if nargout > 5 
    error('Too many output arguments.');
end

if isempty(handles)
    for k = 1:nargout
        varargout{k} = '';
    end
    return;
end

projDir      = handles.projDir;
imgDirList   = handles.imgDirList;
selImgDir    = handles.selImgDir;
firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjDir   = handles.subProjDir;
subProjTitle = handles.subProjTitle;

%Write image path to a file named 'lastProjSetting.txt'.
if isdir(projDir) 
    if isunix==1
        settingsFileName='lastProjSettings_unix.txt';
    elseif ispc==1
        settingsFileName='lastProjSettings_win.txt';
    else
        error('Platform not supported.');
    end
    fid = fopen([handles.projDir filesep settingsFileName],'w');
    if fid ~= -1
       %Write image path.
       if isempty(imgDirList)
          fprintf(fid,'%s\n',['    Image Path: ']);
       else
          imgDir = imgDirList{selImgDir};

          %Switch so that the first img dir is always the selected one.
          imgDirList{selImgDir} = imgDirList{1};
          imgDirList{1}         = imgDir;

          fprintf(fid,'%s\n',['    Image Path: ' imgDir]);
          for k = 2:length(imgDirList)
             fprintf(fid,'%s\n',['              : ' imgDirList{k}]);
          end
       end
       %Write first image name.
       if isempty(firstImgList)
          fprintf(fid,'%s\n',['   First Image: ']);
       else
          firstImg = firstImgList{selImgDir};
          %Switch so that the first img dir is always the selected one.
          firstImgList{selImgDir} = firstImgList{1};
          firstImgList{1}         = firstImg;

          fprintf(fid,'%s\n',['   First Image: ' firstImg]);
          for k = 2:length(imgDirList)
             fprintf(fid,'%s\n',['              : ' firstImgList{k}]);
          end
       end

       for k = 1:numSubProj
          fprintf(fid,'%s\n',[subProjTitle{k} ': ' subProjDir{k}]);
       end
       fclose(fid);
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
if nargout == 5
    varargout{5} = firstImgList;
end

delete(hObject);

function closeRequestFcn(hObject,eventdata,handles)

handles.projDir      = '';
handles.imgDirList   = {};
handles.firstImgList = {};
handles.subProjDir   = {};

% Update handles structure
guidata(hObject, handles);

delete(hObject);


% --- Executes on button press on Cancel.
function cancel_Callback(hObject, eventdata, handles)
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
function Ok_Callback(hObject, eventdata, handles)
% hObject    handle to Ok(see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir    = get(handles.projDirEH,'string');
imgDirList = get(handles.imgDirMH,'string');
selImgDir  = get(handles.imgDirMH,'Value');

if length(imgDirList) == 1
   imgDirList = {};
else
   imgDirList(end) = [];
end

firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjTags  = handles.subProjTags;
subProjDir   = cell(size(subProjTags));

if isempty(projDir)
    warnH = warndlg('No project directory is set.','Warning','modal');
    return;
end

if isempty(imgDirList)
    warnH = warndlg('No image directory is set. ' ,'Warning','modal');
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

subProjOK = ones(1,numSubProj);
for k = 1:numSubProj
   if ~isdir([projDir filesep subProjDir{k}])
      [subProjOK(k),msg,msgID] = mkdir(projDir,subProjDir{k});
   end
   if subProjOK(k) ~= 1
      fprintf(1,'%s\n',['Trouble making ' subProjTag{k} ' directory: ' ...
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
   selImgDir = 1;
   [imgDirList firstImgList subProjDir] = ...
      getProjSetting(projDir,subProjNames);
   handles = updateGUI(handles,projDir,imgDirList,selImgDir, ...
      firstImgList,subProjDir);
end

guidata(hObject,handles);


% --- Executes on return projDir text field.
function projDir_Callback(hObject, eventdata, handles)
% hObject    handle to projDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projDir = get(handles.projDirEH,'string');

if ~isdir(projDir)
   errordlg('The entered path is not a directory.','Input error','modal');
   return;
end

subProjNames = handles.subProjNames;
if ~samdir(handles.projDir,projDir)
    [imgDirList firstImgList subProjDir] = ...
       getProjSetting(projDir,subProjNames);
    handles = updateGUI(handles,projDir,imgDirList,selImgDir, ...
       firstImgList,subProjDir);
end

guidata(hObject,handles);



% --- Executes on button press on Browse.
function delImgDir_Callback(hObject, eventdata, handles)
% hObject    handle to imgDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selImgDir  = get(handles.imgDirMH,'Value');
imgDirList = get(handles.imgDirMH,'String');
if selImgDir == length(imgDirList)
   %' New Image Directory ' is selected. Do nothing.
   return;
end

imgDirList(selImgDir) = [];

firstImgList = handles.firstImgList;
if ~isempty(firstImgList)
   firstImgList(selImgDir) = [];
end

if selImgDir == length(imgDirList) & selImgDir > 1
   selImgDir = selImgDir-1;
end

if selImgDir == length(imgDirList) | isempty(firstImgList)
   set(handles.firstImgTH,'String','');
else
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

set(handles.imgDirMH,'Value',selImgDir);
set(handles.imgDirMH,'String',imgDirList);

if selImgDir == length(imgDirList)
   set(handles.imgDirEH,'String','');
else
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
end

handles.selImgDir    = selImgDir;
handles.imgDirList   = imgDirList(1:end-1);
handles.firstImgList = firstImgList;

guidata(hObject,handles);


function imgDirMenu_Callback(hObject, eventdata, handles)
% hObject    handle to imgDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

firstImgList = handles.firstImgList;

selImgDir  = get(handles.imgDirMH,'Value');
imgDirList = get(handles.imgDirMH,'String');
if selImgDir == length(imgDirList)
   oldDir = pwd;

   if isdir(imgDirList{handles.selImgDir})
      %Go to the parent of last selected image directory.
      cd([imgDirList{handles.selImgDir} filesep '..']);
   end

   %' New Image Directory ' is selected. Choose new image directory and 
   % first image.
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

   %Add the new image dir to 'imgDirList'.
   imgDirList{end+1} = imgDirList{end}; %'New Image Directory'.
   imgDirList{end-1} = pathName;

   set(handles.imgDirMH,'String',imgDirList);
   handles.imgDirList = imgDirList(1:end-1);

   if isempty(firstImgList)
      for k = 1:length(imgDirList)-1
         firstImgList{k} = '';
      end
   end
   firstImgList{selImgDir} = firstImgFileName;
   handles.firstImgList    = firstImgList;
end

handles.selImgDir = selImgDir;
if length(imgDirList) == 1
   set(handles.imgDirEH,'String','');
else
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
end

if isempty(firstImgList)
   set(handles.firstImgTH,'String','');
else
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

guidata(hObject,handles);

function imgDirEdit_Callback(hObject, eventdata, handles)

selImgDir = handles.selImgDir;
imgDirList = handles.imgDirList;

if isempty(imgDirList)
   set(handles.imgDirEH,'String','');
else
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
end



function firstImgBrowse_Callback(hObject, eventdata, handles)

selImgDir    = handles.selImgDir;
imgDirList   = handles.imgDirList;
firstImgList = handles.firstImgList;

if selImgDir == length(imgDirList)
   %' New Image Directory ' is selected in 'imgDir' menu. Do nothing.
   return;
end

oldDir = pwd;

if ~isdir(imgDirList{selImgDir})
   errordlg(['Something is wrong with the current selected image ' ...
      'directory. If it does not exist any more, please delete it ' ...
      'using the delete button.'],'Image selection error','modal');
   return;
end

cd(imgDirList{selImgDir});
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

if ~samdir(pathName,imgDirList{selImgDir})
   errordlg(['Please select image from the current selected ' ...
      'image directory. The selection will be discarded. ' ...
      'Please reselect.'],'Image selection error','modal');
   return;
end

if isempty(firstImgList)
   for k = 1:length(imgDirList)-1
      firstImgList{k} = '';
   end
end

firstImgList{selImgDir} = firstImgFileName;
set(handles.firstImgTH,'String',firstImgFileName);
handles.firstImgList = firstImgList;

guidata(hObject,handles);



function [imgDirList,firstImgList,subProjDir] =  ...
   getProjSetting(projDir,subProjNames)

numSubProj = length(subProjNames);

subProjDir = cell(size(subProjNames));
for k = 1:numSubProj
   subProjDir{k} = '';
end
imgDirList   = {};
firstImgList = {};

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
        textL     = fgetl(fid);
        lineNo    = 1;
        imgLineNo = -1;
        while noProblem & ischar(textL)
            k = findstr(':',textL);
            if isempty(k) | k == 1
                noProblem = 0;
            else
                %Use 'sscanf' to remove space.
                headStr = sscanf(textL(1:k-1),'%s');
                if strcmp(headStr,'ImagePath')
                   imgLineNo = lineNo;
                   if k ~= length(textL)
                      imgDir = sscanf(textL(k+1:end),'%s');
                      if isempty(imgDir)
                         numImgDirs = 0;
                         imgDirList = {};
                      else
                         numImgDirs = 1;
                         imgDirList{numImgDirs} = imgDir;
                      end
                   end
                   lastHead = headStr; 
                elseif strcmp(headStr,'FirstImage')
                   if strcmp(lastHead,'ImagePath') 
                      %The first image has to be written next to image path.
                      if k ~= length(textL)
                         imgName = sscanf(textL(k+1:end),'%s');
                      else
                         imgName = '';
                      end

                      if isempty(imgDirList) 
                         if ~isempty(imgName)
                            noProblem = 0;
                         else
                            firstImgList = {};
                         end
                      else
                         numImgDirs = 1;
                         firstImgList{numImgDirs} = imgName;
                      end
                   else
                      noProblem = 0;
                   end
                   lastHead = headStr; 
                elseif isempty(headStr) 
                   if strcmp(lastHead,'ImagePath') 
                      if isempty(imgDirList)
                         noProblem = 0;
                      elseif k ~= length(textL)
                         %Get next image path
                         imgDir = sscanf(textL(k+1:end),'%s');
                      else
                         imgDir = '';
                      end

                      if isempty(imgDir)
                         noProblem = 0;
                      else
                         numImgDirs = numImgDirs+1;
                         imgDirList{numImgDirs} = imgDir;
                      end
                   elseif strcmp(lastHead,'FirstImage') 
                      if k ~= length(textL)
                         %Get next first image name.
                         imgName = sscanf(textL(k+1:end),'%s');
                      else
                         imgName = '';
                      end
                      numImgDirs = numImgDirs+1;

                      if numImgDirs > length(imgDirList)
                         noProblem = 0;
                      else
                         firstImgList{numImgDirs} = imgName;
                      end
                   else
                      noProblem = 0;
                   end
                else
                   j1 = strmatch(headStr,subProjNames,'exact');
                   if isempty(j1)
                      noProblem = 0;
                   elseif k ~= length(textL)
                      subProjDir{j1} = sscanf(textL(k+1:end),'%s');
                   end
                   lastHead = headStr; 
                end
            end
            
            textL  = fgetl(fid);
            lineNo = lineNo+1;
        end
        if ~noProblem
           warnH = warndlg(['Project setting file is corrupted ' ...
              'and will be ignored.'],'Warning','modal'); 
        end
        fclose(fid);
    end
end

if noProblem == 1 & length(firstImgList) ~= length(imgDirList)
   if ~isempty(firstImgList)
      noProblem = 0;
      warnH = warndlg(['Number of image directories and first images ' ...
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
end

