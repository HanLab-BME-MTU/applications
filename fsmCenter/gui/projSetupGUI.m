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

handles.selImgDir = 1;
handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
handles = updateGUI(handles);

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

function handles = updateGUI(handles)

numSubProj  = handles.numSubProj;
subProjTags = handles.subProjTags;
subProjMH   = handles.subProjMH;
subProjTFH  = handles.subProjTFH;

projDir       = handles.projDir;
imgDirList    = handles.imgDirList;
firstImgList  = handles.firstImgList;
unix_imgDrive = handles.unix_imgDrive;
win_imgDrive  = handles.win_imgDrive;
subProjDir    = handles.subProjDir;
selImgDir     = handles.selImgDir;

if isempty(firstImgList)
   firstImg = '';
else
   firstImg = firstImgList{selImgDir};
end

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
      subProjDir{k},subProjTags{k});
end


% --- Outputs from this function are returned to the command line.
function varargout = projSetupGUI_OutputFcn(hObject, eventdata, handles)
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

projDir      = handles.projDir;
imgDirList   = handles.imgDirList;
selImgDir    = handles.selImgDir;
firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjDir   = handles.subProjDir;
subProjTitle = handles.subProjTitle;
physiParam   = handles.physiParam;

unix_imgDrive = [];
win_imgDrive  = [];
%Write image path to a file named 'lastProjSetting.mat'.
if isdir(projDir)
   if ~isempty(imgDirList)
      %Switch so that the first img dir is always the selected one.
      imgDir = imgDirList{selImgDir};
      imgDirList{selImgDir} = imgDirList{1};
      imgDirList{1}         = imgDir;

      if isunix == 1
         filesepInd = findstr('/',imgDir);
         mntInd = findstr('/mnt/',imgDir);
         if isempty(mntInd)
            unix_imgDrive = imgDir(1:filesepInd(2)-1);
         else
            unix_imgDrive = imgDir(1:filesepInd(3)-1);
         end
      elseif ispc == 1
         colonInd = findstr(':',imgDir);
         win_imgDrive = imgDir(1:colonInd(1));
      end
   end

   if ~isempty(firstImgList)
      %Switch so that the first img dir is always the selected one.
      firstImg = firstImgList{selImgDir};
      firstImgList{selImgDir} = firstImgList{1};
      firstImgList{1}         = firstImg;
   end

   if ~isempty(physiParam)
      if length(physiParam) < length(firstImgList)
          defPhysiParam = physiParam{end};
          for k = length(physiParam)+1:length(firstImgList)
             physiParam{k} = defPhysiParam;
          end
      end
      selPhysiParam = physiParam{selImgDir};
      physiParam{selImgDir} = physiParam{1};
      physiParam{1}         = selPhysiParam;
   end


    projSettings.projDir         = projDir;

    settingsMatFile = [projDir filesep 'lastProjSettings.mat'];
    if isunix==1
        projSettings.unix_imgDirList = imgDirList;

        if samdir(handles.unix_imgDrive,unix_imgDrive)
           if ~isempty(handles.win_imgDrive)
              projSettings.win_imgDirList = ...
                 dirUnix2PC(imgDirList,handles.win_imgDrive);
           end
        end
        settingsFileName='lastProjSettings_unix.txt';
    elseif ispc==1
        projSettings.win_imgDirList = imgDirList;
        if samdir(handles.win_imgDrive,win_imgDrive)
           if ~isempty(handles.unix_imgDrive)
              projSettings.unix_imgDirList = ...
                 dirPC2Unix(imgDirList,handles.unix_imgDrive);
           end
        end
        settingsFileName='lastProjSettings_win.txt';
    else
        error('Platform not supported.');
    end
    projSettings.firstImgList    = firstImgList;
    projSettings.subProjDir      = subProjDir;
    projSettings.physiParam      = physiParam;
    
    save(settingsMatFile,'projSettings');
    
    %We also save a text file of the project settings.
    fid = fopen([handles.projDir filesep settingsFileName],'w');
    if fid ~= -1
       %Write image path.
       if isempty(imgDirList)
          fprintf(fid,'%s\n',['    Image Path: ']);
       else
          fprintf(fid,'%s\n',['    Image Path: ' imgDirList{1}]);
          for k = 2:length(imgDirList)
             fprintf(fid,'%s\n',['              : ' imgDirList{k}]);
          end
       end
       %Write first image name.
       if isempty(firstImgList)
          fprintf(fid,'%s\n',['   First Image: ']);
       else
          fprintf(fid,'%s\n',['   First Image: ' firstImgList{1}]);
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
if nargout >= 5
    varargout{5} = firstImgList;
end
if nargout == 6
    varargout{6} = physiParam;
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
if isempty(projDir)
    warnH = warndlg('No project directory is set.','Warning','modal');
    return;
end

selImgDir  = handles.selImgDir;

if selImgDir == 0;
   warnH = warndlg('No image directory is set. ' ,'Warning','modal');
   return;
end

imgDirList = handles.imgDirList;

firstImgList = handles.firstImgList;
numSubProj   = handles.numSubProj;
subProjTags  = handles.subProjTags;
subProjDir   = cell(size(subProjTags));

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
subProjTags  = handles.subProjTags;
if ~samdir(handles.projDir,projDir)
   handles.selImgDir = 1;
   handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
   handles = updateGUI(handles);
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
subProjTags  = handles.subProjTags;
if ~samdir(handles.projDir,projDir)
    handles = getProjSetting(handles,projDir,subProjNames,subProjTags);
    handles = updateGUI(handles);
end

guidata(hObject,handles);



% --- Executes on button press on Browse.
function delImgDir_Callback(hObject, eventdata, handles)
% hObject    handle to imgDirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selImgDir  = handles.selImgDir;
imgDirList = handles.imgDirList;

if selImgDir == 0
   return;
end

imgDirList(selImgDir) = [];

firstImgList = handles.firstImgList;
if ~isempty(firstImgList)
   firstImgList(selImgDir) = [];
end

if selImgDir > length(imgDirList)
   selImgDir = selImgDir-1;
end

if selImgDir == 0
   %All image directories are deleted.
   set(handles.imgDirMH,'Enable','off');
   set(handles.imgDirMH,'String',{'--- Empty ---'});
   set(handles.imgDirMH,'Value',1);
   set(handles.imgDirEH,'String','');
   set(handles.firstImgTH,'String','');
else
   %Update the menu list for image directory list.
   set(handles.imgDirMH,'Enable','on');
   set(handles.imgDirMH,'String',imgDirList);
   set(handles.imgDirMH,'Value',selImgDir);

   %Update the edit field for active (selected) image directory and the text
   %field for first image file name.
   set(handles.imgDirEH,'String',imgDirList{selImgDir});
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

handles.selImgDir    = selImgDir;
handles.imgDirList   = imgDirList;
handles.firstImgList = firstImgList;

guidata(hObject,handles);


function imgDirMenu_Callback(hObject, eventdata, handles)
% hObject    handle to imgDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

firstImgList = handles.firstImgList;
imgDirList   = handles.imgDirList;

selImgDir  = get(handles.imgDirMH,'Value');

%Update the gui field for active (selected) image directory.
set(handles.imgDirEH,'String',imgDirList{selImgDir});

if isempty(firstImgList)
   set(handles.firstImgTH,'String','');
else
   set(handles.firstImgTH,'String',firstImgList{selImgDir});
end

handles.selImgDir = selImgDir;
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
while ~imgDirExist & k <= length(imgDirList)
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

subProjDir = cell(size(subProjNames));
for k = 1:numSubProj
   subProjDir{k} = '';
end
imgDirList   = {};
firstImgList = {};

unix_imgDrive = [];
win_imgDrive  = [];

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
        
        if isfield(projSettings,'physiParam')
           physiParam = projSettings.physiParam;
        else
           for k = 1:length(firstImgList)
              physiParam{k} = getDefFsmPhysiParam;
           end
        end
        if isfield(projSettings,'unix_imgDirList')
            unix_imgDirList = projSettings.unix_imgDirList;
            filesepInd = findstr('/',unix_imgDirList{1});
            mntInd = findstr('/mnt/',unix_imgDirList{1});
            if isempty(mntInd)
                unix_imgDrive = unix_imgDirList{1}(1:filesepInd(2)-1);
            else
                unix_imgDrive = unix_imgDirList{1}(1:filesepInd(3)-1);
            end
        end
        if isfield(projSettings,'win_imgDirList')
            win_imgDirList = projSettings.win_imgDirList;
            colonInd       = findstr(':',win_imgDirList{1});
            win_imgDrive   = win_imgDirList{1}(1:colonInd(1));
        end
        if isunix == 1
            if ~isempty(unix_imgDrive)
                imgDirList   = unix_imgDirList;
                
                noProblem = 1;
            elseif ~isempty(win_imgDrive)                
                tryAgain = 'Yes';
                prompt = sprintf(['Last project is set up in Windows. \n' ...
                    'The image drive letter is ' win_imgDrive  '.\n' ...
                    'Please enter image drive name in Unix format:']);
                while strcmp(tryAgain,'Yes')
                    answer = inputdlg(prompt,'title',1,{''});
                    unix_imgDrive = answer{1};

                    %Convert image directories to Unix format.
                    imgDirList = dirPC2Unix(win_imgDirList,unix_imgDrive);
                    if ~isdir(imgDirList{1})
                        question = ['Invalid drive name. Do you want to try again? ' ...
                            'If no, the old image directories will be removed. ' ...
                            'Be cautious!!!'];
                        tryAgain = questdlg(question,'Alert','Yes','No','Yes');
                    else
                        projSettings.unix_imgDirList = imgDirList;
                        noProblem = 1;
                        tryAgain  = 'No';
                    end
                end
            else
                error('The project setting file is corrupted. You need to reset project.');
            end
        elseif ispc == 1
            if ~isempty(win_imgDrive)
                imgDirList = win_imgDirList;
                
                noProblem = 1;
            elseif ~isempty(unix_imgDrive)
                tryAgain = 'Yes';
                prompt = sprintf(['Last project is set up in Unix platform. \n' ...
                    'The image drive name is ' unix_imgDrive  '.\n' ...
                    'Please enter image drive letter in PC format:']);
                while strcmp(tryAgain,'Yes')
                    answer = inputdlg(prompt,'title',1,{''});
                    win_imgDrive = answer{1};

                    %Convert image directories to PC format.
                    imgDirList = dirUnix2PC(unix_imgDirList,win_imgDrive);
                    if ~isdir(imgDirList{1})
                        question = ['Invalid drive letter. Do you want to try again? ' ...
                            'If no, the old image directories will be removed. ' ...
                            'Be cautious!!!'];
                        tryAgain = questdlg(question,'Alert','Yes','No','Yes');
                    else
                        projSettings.win_imgDirList = imgDirList;
                        noProblem = 1;
                        tryAgain  = 'No';
                    end
                end
            else
                error('The project setting file is corrupted. You need to reset project.');
            end
        else
             error('Platform not supported.');
        end       
    else
        %The project setting used to be saved in text files and it is not
        %convinient for coding. For backward compatibility, we transfer all
        %the old setting file to mat file.
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

            if noProblem
               for k = 1:length(imgDirList)
                  physiParam{k} = getDefFsmPhysiParam;
               end
            end
        end
    end
else
   projDir = '';
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

   win_imgDrive  = [];
   unix_imgDrive = [];
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

handles.projDir       = projDir;
handles.imgDirList    = imgDirList;
handles.firstImgList  = firstImgList;
handles.unix_imgDrive = unix_imgDrive;
handles.win_imgDrive  = win_imgDrive;
handles.subProjDir    = subProjDir;
handles.physiParam    = physiParam;

if isempty(imgDirList)
   handles.selImgDir = 0;
else
   handles.selImgDir = 1;
end
