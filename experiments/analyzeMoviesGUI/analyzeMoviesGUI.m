function varargout = analyzeMoviesGUI(varargin)
% ANALYZEMOVIESGUI M-file for analyzeMoviesGUI.fig
%      ANALYZEMOVIESGUI, by itself, creates a new ANALYZEMOVIESGUI or raises the existing
%      singleton*.
%
%      H = ANALYZEMOVIESGUI returns the handle to a new ANALYZEMOVIESGUI or the handle to
%      the existing singleton*.
%
%      ANALYZEMOVIESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZEMOVIESGUI.M with the given input arguments.
%
%      ANALYZEMOVIESGUI('Property','Value',...) creates a new ANALYZEMOVIESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyzeMoviesGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyzeMoviesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyzeMoviesGUI

% Last Modified by GUIDE v2.5 04-Feb-2003 17:22:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzeMoviesGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzeMoviesGUI_OutputFcn, ...
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


% --- Executes just before analyzeMoviesGUI is made visible.
function analyzeMoviesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyzeMoviesGUI (see VARARGIN)

% Choose default command line output for analyzeMoviesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analyzeMoviesGUI wait for user response (see UIRESUME)
% uiwait(handles.AMG);


% --- Outputs from this function are returned to the command line.
function varargout = analyzeMoviesGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function dirListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dirListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in dirListBox.
function dirListBox_Callback(hObject, eventdata, handles)
% hObject    handle to dirListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dirListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dirListBox


% --- Executes on button press in amg_load_PB.
function amg_load_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_load_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check if default biodata-dir exists
mainDir = getenv('BIODATA');
if isempty(mainDir)
    mainDir = getenv('HOME');
end

%get project directory
oldPath = pwd;
cd(mainDir);
[jobFile,jobPath] = uigetfile('job*','select a job file');


%return if user has selected nothing
if jobFile==0
    cd(oldPath);
    return
end
cd(jobPath);
load(jobFile);
cd(oldPath);

if isfield(job,'jobs2run')
    
    %ensure portability of jobfiles between windows and linux
    nJobs = length(job);
        for i = 1:nJobs
            fileSepIdx = [findstr(job(i).projProperties.dataPath,'\'),...
                findstr(job(i).projProperties.dataPath,'/')];
            job(i).projProperties.dataPath(fileSepIdx) = filesep;
        end
    
    %find strings for dirList
    jobCell = struct2cell(job);
    jobNames = fieldnames(job);
    nameIdx = strmatch('projName',jobNames);
    nameList = jobCell(nameIdx,:)';
    set(handles.dirListBox,'Value',1);
    
    dirList = get(handles.dirListBox,'String');
    if ~iscell(dirList)
        set(handles.dirListBox,'String',nameList);
        handles.job = job;
        dirListSelect = 1;
    else
        keep = questdlg('Do you want to keep the other projects?','','yes','no','no');
        if strcmp(keep,'no')
            set(handles.dirListBox,'String',nameList);
            handles.job = job;
            dirListSelect = 1;
        else
            dirListSelect = 1+size(dirList,1);
            dirList = [dirList;nameList];
            set(handles.dirListBox,'String',dirList);
            handles.job = [handles.job,job];
        end
    end
    set(handles.dirListBox,'Value',dirListSelect);
    dirListLength = size(get(handles.dirListBox,'String'),1);
    guidata(hObject,handles);
    
    ans = questdlg('Do you want to update properties?','','yes','no','yes');
    if strcmp(ans,'yes')
        editProperties(handles,[dirListSelect:dirListLength]);
    end
    
end

% --- Executes on button press in amg_editProp_PB.
function amg_editProp_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_editProp_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

editProperties(handles);

% --- Executes on button press in amg_guiProp_PB.
function amg_guiProp_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_guiProp_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ans = questdlg('Do you really think this program needs more options?','What the ...!?','yes','no, of course not!','no, of course not!');
if strcmp('yes',ans)
    h = helpdlg('Great! Why don''t you do some work on it yourself?');
    uiwait(h);
    return
else
    h = helpdlg('But why click on this button?','Correct...');
    uiwait(h);
end

% --- Executes on button press in amg_cancel_PB.
function amg_cancel_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_cancel_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%confirm user input
really = questdlg('Do you really want to quit without running/saving jobs?','Quit?','yes','no','yes');
if strcmp(really,'no')
    return %end evaluation here
end

delete(handles.AMG);

% --- Executes on button press in amg_deletePB.
function amg_deletePB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_deletePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

dirList = get(handles.dirListBox,'String');
activeJobNum = get(handles.dirListBox,'Value');
if ~iscell(dirList)
    return
    %no cell, no job, no delete
end

if length(dirList)==1
    dirList = char('No project loaded');
else
    dirList(activeJobNum) = [];
end
%store new dirList
set(handles.dirListBox,'Value',1);
set(handles.dirListBox,'String',dirList);
%store job data
handles.job(activeJobNum) = [];
guidata(hObject,handles);

% --- Executes on button press in amg_add_PB.
function amg_add_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_add_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%check if default biodata-dir exists

mainDir = cdBiodata(2);

%guihandle
guiH = handles.AMG;

%get project name. If the name of the movie and the name of the directory
%are not the same, create subdirectory with projectname = moviename and move
%movie and logfile
movieNames = searchFiles('(.r3[cd]|[dr]3d.dv)','(log|part|bad)','ask');
selectedIdx = listSelectGUI(movieNames(:,1));
nSelected = length(selectedIdx);
movieNames = movieNames(selectedIdx,:);

% if there are any movieNames, we loop, otherwise, we automatically return
for iMovie = 1:nSelected
    fileName = movieNames{iMovie,1};
    pathName = movieNames{iMovie,2};
if strcmpi(fileName(end-5:end),'3D.dv')
    % change fileName to *.r3d
    destFileName = [fileName(1:end-7),'.r3d'];
elseif strcmp(fileName(end-3:end-1),'.r3')
    % we're happy
    destFileName = fileName;
else
    h=errordlg('This file extension is not recognized. Please notify the authorities')
    uiwait(h);
    return
end

% %make sure user has not inadvertedly selected the r3d instead of the r3c
% if strcmp(fileName(end),'d')&exist([pathName,fileName(1:end-4),'_corr.r3c'])
%     ans = questdlg('Are you sure you want to load the uncorrected movie?','','load uncorrected','load corrected','cancel','load corrected');
%     switch ans
%         case 'load uncorrected'
%             %leave everything as is
%         case 'load corrected'
%             %change fileName
%             fileName = [fileName(1:end-4),'_corr.r3c'];
%             destFileName = fileName;
%         otherwise
%             %cancel or []
%             return
%     end
% end

%if last part of pathname differs from moviename: create project directory
%and move movie and logfile; else read projectName
projName = destFileName(1:end-4);

%get extension of movie file
%fileExtension = destFileName(end-3:end); %.r3d/.r3c

%get name of directory
listSep = findstr(pathName,filesep);
% lastSep = listSep(end-1); %with uigetfile, the path ends with a filesep!
% dirName = pathName(lastSep+1:end-1); %ditto
% with the new search/list, there are no trailing fileseps anymore
dirName = pathName(listSep(end)+1:end);

if strcmpi(dirName(end-2:end),'bad')
    % then don't do any calculations
    disp(sprintf('%s is labelled ''bad'' - will not be analyzed',fileName))
else
if strcmp(projName,dirName)|strcmp(projName,[dirName,'_corr'])
    %projName==dirName or projName==dirName_corr
    fullPathName = pathName; 
else
    %projName~=dirName
    cd(pathName);
    % add filesep to pathName
    if ~strcmp(pathName(end),filesep)
        pathName = [pathName,filesep];
    end
    %create new directory
    mkdir(projName);
    %store directory path
    fullPathName = [pathName,projName,filesep]; 
    %move movie and rename if necessary
    if ispc
        movefile([pathName,fileName],[fullPathName,destFileName]);
    else %there is a bug in linux with copying
        system(['mv ' [pathName,fileName] ' ' [fullPathName,destFileName]]);
    end
    
    %move logfile
    logfileName = [pathName,fileName,'.log'];
    if exist(logfileName)
        if ispc
            movefile([pathName,fileName,'.log'],fullPathName);
        else
            system(['mv ' [pathName,fileName,'.log'] ' ' fullPathName]);
        end
    end
end
cd(fullPathName);

%see if project is in mainDir, do it case-unsensitive
if strcmpi(fullPathName(1:length(mainDir)),[mainDir])==1
    relPathName = fullPathName(length(mainDir)+2:end); %begins without filesep, ends with none
else
    ans = questdlg('moviefile is not in main file structure or you have no env-var BIODATA. There will be problems sharing your data',...
        'WARNING','Continue','Cancel','Cancel');
    if strcmp(ans,'Cancel')
        return %end evaluation here
    else
        relPathName = fullPathName(1:end-1); %eliminate last filesep
    end
end


%write project name into listbox
dirList = get(handles.dirListBox,'String');
if ~iscell(dirList)
    dirList = cellstr(projName);
else
    dirList(end+1) = cellstr(projName);
end
set(handles.dirListBox,'String',dirList);

%number of already selected projects
projNum = length(dirList); 

%select current project
set(handles.dirListBox,'Value',projNum);


%write job data
%handles.job(projNum).dirName = directoryName; 
handles.job(projNum).projName = projName;
handles.job(projNum).projProperties.dataPath = relPathName;
handles.job(projNum).projProperties.datafileName = [];
handles.job(projNum).projProperties.status = 0;
handles.job(projNum).dataProperties = []; %to be defined in editPropertiesGUI
handles.job(projNum).createNew = 1; %default
handles.job(projNum).lastName = [];
handles.job(projNum).mainDir = [];
handles.job(projNum).mainSaveDir = [];
%handles.job(projNum).fileExtension = fileExtension;
handles.job(projNum).correctBackground = [];

%save data
guidata(gcbo,handles);

%launch property window
editProperties(handles);
end % check 'bad'
end % end loop movies

%-----------------------------------------------------------

% --- Executes on button press in amg_run_PB.
function amg_run_PB_Callback(hObject, eventdata, handles)
% hObject    handle to amg_run_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

%check if there is a jobfile at all
if ~isfield(handles,'job')
    h = errordlg('please specify a job first');
    uiwait(h);
    return
end

%get jobfile
fileName = ['job-',nowString];
job = handles.job;
job(1).lastName = fileName;

%make sure there are dataProperties defined for all jobs
for i = 1:length(job);
    if isempty(job(i).dataProperties)
        h = errordlg(['There are no properties defined for job no. ',num2str(i),'!']);
        uiwait(h);
        return
    end
end

%check if default biodata-dir exists
mainDir = getenv('BIODATA');
if isempty(mainDir)
    mainDir = getenv('HOME');
end
job(1).mainDir = mainDir;
job(1).mainSaveDir = [];

%ask for confirmation
ans = questdlg('run and/or save job?','choose wisely!','run only','save&run','save only','run only');
switch ans
    case ''
       %do nothing
       return
   case 'run only'
       %no problemo
   case 'save&run'
       cd(mainDir);
       [fileName,pathName] = uiputfile('job-*','save jobfile'); %careful: pathName includes filesep
       if ~strcmp(pathName,[mainDir,filesep])
           ans2 = questdlg('do you want to store your jobdata in this folder, too?','','yes','no','no');
           if strcmp(ans2, 'yes')
               job(1).mainSaveDir = pathName(1:end-1);
           end
       end
       if pathName == 0
           return
       end
       save([pathName,fileName],'job');
       
   case 'save only'
       cd(mainDir);
       [fileName,pathName] = uiputfile('job-*','save jobfile'); %careful: pathName includes filesep
       if ~strcmp(pathName,[mainDir,filesep])
           ans2 = questdlg('do you want to store your jobdata in this folder, too?','','yes','no','no');
           if strcmp(ans2, 'yes')
               job(1).mainSaveDir = pathName(1:end-1);
               mainDir = pathName(1:end-1);
           end
       end
       save([pathName,fileName],'job');
       return %do not run job
end


%save to mainDir
if isempty(job(1).mainSaveDir)
    save([mainDir,filesep,fileName],'job');
else
    save([job(1).mainSaveDir,filesep,fileName],'job');
end

%close figure
delete(handles.AMG);

runCtBatch(job);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function editProperties(handles,dirListVector)
%wrapper for calling editPropertiesGUI

%if there is no job loaded, do not open editPropertiesGUI
dirList = get(handles.dirListBox,'String');
if ~iscell(dirList)
    return
end

%if another GUI is open, close it
h = findobj('Tag','EPGUI');
if ~isempty(h)
    close(h);
end

%get number of projects in dirListBox
dirListLength = size(get(handles.dirListBox,'String'),1);
if dirListLength==1
    handles.apply2all = -1;
else
    handles.apply2all = 0;
end

switch nargin
    case 1 %edit current project
        handles.acceptAll = -1;
        guidata(handles.AMG,handles);
        set(handles.AMG,'Visible','off')
        epH = editPropertiesGUI;
        s = warning('query', 'all');
        warning off all;
        uiwait(epH);
        set(handles.AMG,'Visible','on')
        warning(s); %turn back on the warnings
        handles = guidata(handles.AMG);

    case 2
        i = dirListVector(1);
        done = 0;
        if length(dirListVector)==1
            handles.acceptAll = -1;
        else
            handles.acceptAll = 0;
        end
        
        while i<=dirListVector(end)&done~=1
            %set current job to update
            set(handles.dirListBox,'Value',i);
            %erase job-properties, but remember cropInfo
            if isfield(handles.job(i).dataProperties,'crop')
                cropInfo = handles.job(i).dataProperties.crop;
            else
                cropInfo = [];
            end
            handles.job(i).dataProperties = [];
            handles.job(i).dataProperties.crop = cropInfo;
            handles.job(i).jobs2run = 0;
            handles.job(i).createNew = 0;
            
            guidata(handles.AMG,handles);
            epH = editPropertiesGUI;
            switch handles.acceptAll %if acceptAll: do not need to show GUI
                case 0
                    set(handles.AMG,'Visible','off')
                    uiwait(epH);
                    set(handles.AMG,'Visible','on')
                case -1
                    set(handles.AMG,'Visible','off')
                    uiwait(epH);
                    set(handles.AMG,'Visible','on')
                case 1
                    editPropertiesGUI('edit_OK_PB_Callback',epH,[],guidata(epH));
            end
            
            handles = guidata(handles.AMG);
            
            %check if quit loop to set current for all
            if handles.apply2all>0
                handles.acceptAll = 1;
            end
            i = i+1;
        end
end

handles = guidata(handles.AMG);

if handles.apply2all>0 %could be 1 or 2
    
    % read the properties of myJob, and apply it to all the other projects. If
    % they are not as far as myJob yet, they will have to do the catching
    % up with the default parameters. If they are further advanced, the
    % concerned jobs are rerun. (It all depends on what job settings the user wants
    % to apply to all!)
    
    myJobNum = get(handles.dirListBox,'Value');
    myJob = handles.job(myJobNum);
    myJobs2run = bsum2bvec(myJob.jobs2run);
    
    % get individual job properties of myJob
    myFieldNames = '';
    if any(myJobs2run == 1) %filter movie. Copy correct background
        myFieldNames = [myFieldNames;{'correctBackground'}];
    end
    if any(myJobs2run==2)
        myFieldNames = [myFieldNames;{'dataProperties.CH_MAXSLOPE';'dataProperties.F_TEST_PROB'}];
    end
    if any(myJobs2run==4)
        myFieldNames = [myFieldNames;{'dataProperties.IDopt'}];
    end
%     if any(myJobs2run==16)
%         myFieldNames = [myFieldNames];
%     end
    if handles.apply2all==2 %update movie properties, too
        myFieldNames = [myFieldNames;{'dataProperties.cellCycle';...
                                      'dataProperties.strains';...
                                      'dataProperties.drugs';...
                                      'dataProperties.temperature'}];
    end
    
    %set job properties
    for i = 1:dirListLength
        if i~=myJobNum
            
            % run the same jobs as in myJob; if there has been any previous
            % analysis, discard those results!
            handles.job(i).jobs2run = myJob.jobs2run;
            
            % check if project is "ready" to run jobs
            % or whether we have to run some jobs first, before we get to
            % the point where the current project would start
            
            statVec = bsum2bvec(handles.job(i).projProperties.status);
            if isempty(statVec) 
                statVec = 0; % I do not want to risk changing bsum2bvec for this
            end
            if ~isempty(myJobs2run) %only look at jobs2run if there is any
                jobi = myJobs2run(1);
                while statVec(end)<jobi/2&jobi>2
                    jobi = jobi/2;
                    handles.job(i).jobs2run = handles.job(i).jobs2run+jobi;
                end
                % check whether we have to force createNew
                if myJob.createNew & myJobs2run(1)<=2 
                    handles.job(i).createNew = 1;
                end
            end
            
            % set all the necessary job properties
            handles.job(i).eraseAllPrev = myJob.eraseAllPrev;
            for j = 1:size(myFieldNames,1)
                eval(['handles.job(i).',myFieldNames{j},' = myJob.',myFieldNames{j},';']);
            end
        end
    end
end

%update guidata
guidata(handles.AMG,handles);
