function varargout = plusTipGroupAnalysis(varargin)
% PLUSTIPGROUPANALYSIS M-file for plusTipGroupAnalysis.fig
%      PLUSTIPGROUPANALYSIS, by itself, creates a new PLUSTIPGROUPANALYSIS or raises the
%      existing
%      singleton*.
%
%      H = PLUSTIPGROUPANALYSIS returns the handle to a new PLUSTIPGROUPANALYSIS or the handle to
%      the existing singleton*.
%
%      PLUSTIPGROUPANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPGROUPANALYSIS.M with the given input arguments.
%
%      PLUSTIPGROUPANALYSIS('Property','Value',...) creates a new PLUSTIPGROUPANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipGroupAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipGroupAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipGroupAnalysis

% Last Modified by GUIDE v2.5 30-Jun-2011 10:59:03


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plusTipGroupAnalysis_OpeningFcn, ...
    'gui_OutputFcn',  @plusTipGroupAnalysis_OutputFcn, ...
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


% --- Executes just before plusTipGroupAnalysis is made visible.
function plusTipGroupAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipGroupAnalysis (see VARARGIN)

% Choose default command line output for plusTipGroupAnalysis
handles.output = hObject;


% for "select projects" pushbutton
handles.projList=[]; % select projects pushbutton
handles.loadProjList = 0; % load projList checkbox
handles.getStr = 0; % narrow down list checkbox
handles.projData=[]; % if one project is selected, projData will be retrieved

% for "create groups" pushbutton
userData=get(handles.figure1,'UserData');
userData.groupList=[]; % also select groups pushbutton
set(handles.figure1,'UserData',userData);

%place image onto the axes, remove tick marks
pic=imread('pTT_logo_sm.png');
set(handles.figure1,'CurrentAxes',handles.logoAxes);
image(pic);
axis off

set(handles.getHelpPush,'CData',imread('help_icon.png'),'Callback',...
    @(hObject,eventdata)open('plusTipGroupAnalysis.pdf'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plusTipGroupAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipGroupAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plusTipGuiSwitch(hObject,eventdata,handles,'resetButton');


% --- Executes on button press in getProjPush.
function getProjPush_Callback(hObject, eventdata, handles)
% hObject    handle to getProjPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plusTipGuiSwitch(hObject,eventdata,handles,'getProjPush');

if ~isempty(handles.projList)
    
    % here we do NOT filter out any sub-directories
    a=projList2Cell(handles.projList);
    a=a(:,1);

    % allow multiple projects to be selected
    if isempty(a)
        selection=[];
    else
        [selection,selectionList]=listSelectGUI(a,[],'move',1);
    end

    handles.dataDir=[];
    % if only one project was selected, save projData info and get data
    if isempty(selection)
        msgbox('No projects selected or tracking has not been completed.')
        handles.projList=[];
    elseif size(selection,1)==1
        handles.projList=handles.projList(selection,1);
        handles.dataDir=selectionList{1,1};
        p=load([handles.dataDir filesep 'meta' filesep 'projData.mat']);
        handles.projData=p.projData;
        set(handles.anDirEdit,'String',handles.projData.anDir);
        
        assignin('base','projData',handles.projData)
    else
        handles.projList=handles.projList(selection,1);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);
guidata(hObject, handles);


% --- Executes on button press in getQueryStr_Check.
function getQueryStr_Check_Callback(hObject, eventdata, handles)
% hObject    handle to getQueryStr_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getQueryStr_Check
handles.getStr=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in getProjListFile_check.
function getProjListFile_check_Callback(hObject, eventdata, handles)
% hObject    handle to getProjListFile_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getProjListFile_check
handles.loadProjList=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in selectOutputDirPush.
function selectOutputDirPush_Callback(hObject, eventdata, handles)
% hObject    handle to selectOutputDirPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1,'UserData');

if isfield(handles,'dataDir')
    dirStart=handles.dataDir;
else
    dirStart=pwd;
end

outDir=uigetdir(dirStart,'Please select OUTPUT directory for group analyis');

if isequal(outDir,0), return; end

userData = get(handles.figure1,'UserData');
userData.saveDir=outDir;
set(handles.figure1,'UserData',userData);

% --- Executes on button press in pickGroupsPush.
function pickGroupsPush_Callback(hObject, eventdata, handles)

autoGrp =get(handles.checkbox_autoGrp,'Value');

[handles.groupList]=plusTipPickGroups(autoGrp,[],...
    handles.projList,1);
assignin('base','groupList',handles.groupList);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_poolData.
function pushbutton_poolData_Callback(hObject, eventdata, handles)

userData = get(handles.figure1,'UserData');
doBtw=get(handles.checkbox_doBtw,'Value');
doWtn=get(handles.checkbox_doWtn,'Value');
doPlot=get(handles.checkbox_doPlot,'Value');
remBegEnd=get(handles.checkbox_remBegEnd,'Value');
[groupData]=plusTipPoolGroupData(handles.groupList,...
    userData.saveDir,doBtw,doWtn,doPlot,remBegEnd);
