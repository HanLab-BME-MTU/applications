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

% Last Modified by GUIDE v2.5 23-Sep-2011 15:14:48


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

% Set-up popup-menus
testList = {'t-test of the means';'Wilcoxon ranksum test';'Kolmogorov-Smirnov test (K-S test)';...
    'Mean substracted K-S test';'Median substracted K-S test';...
    'Permutation t-test of the means';'Calibrated mean subtracted K-S test'};
testValues=[1 2 10 11 12 20 21];
set(handles.popupmenu_testID1,'String',testList,'UserData',testValues);
set(handles.popupmenu_testID2,'String',testList,'UserData',testValues);
uipanel_analysisMode_SelectionChangeFcn(hObject, eventdata, handles)

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

% two-sided p-value: proportion of abs(delta) values greater than deltaPop
%pValue = sum(abs(delta)>deltaPop)/nReps;

% calculate the one-sided p-value
pValue = 1-normcdf(deltaPop,mean(delta),std(delta));

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
        assignin('base','projData',handles.projData)
    else
        handles.projList=handles.projList(selection,1);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
    return
end
temp=projList2Cell(handles.projList);
assignin('base','selectedProjects',temp);

% Launch the group creation routine
autoGrp =get(handles.checkbox_autoGrp,'Value');
userData.groupList=plusTipPickGroups(autoGrp,[],handles.projList,1);
if isempty(userData.groupList),return; end
userData.groupData=plusTipExtractGroupData(userData.groupList);
set(handles.figure1,'UserData',userData);
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

if isfield(handles,'dataDir')
    dirStart=handles.dataDir;
else
    dirStart=pwd;
end

outDir=uigetdir(dirStart,'Please select OUTPUT directory for group analyis');

if isequal(outDir,0), return; end

userData = get(handles.figure1,'UserData');
set(handles.edit_outDir,'String',outDir);
set(handles.figure1,'UserData',userData);


% --- Executes on button press in pushbutton_loadGroup.
function pushbutton_loadGroup_Callback(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat','Select the group list to open',pwd);
if ~any([file,path]), return; end
try
    userData=get(handles.figure1,'UserData');    
    s=load([path file]);
    userData.groupList = s.groupList;
    userData.groupData=plusTipExtractGroupData(userData.groupList);
    set(handles.figure1,'UserData',userData);   
catch ME
    throw(ME)
end

% --- Executes on button press in checkbox_doPlot.
function checkbox_doPlot_Callback(hObject, eventdata, handles)

if get(hObject,'Value'); enable='on'; else enable='off'; end
set(get(handles.uipanel_histogram,'Children'),'Enable',enable);


% --- Executes on button press in pushbutton_analyzeGroups.
function pushbutton_analyzeGroups_Callback(hObject, eventdata, handles)

userData = get(handles.figure1,'UserData');
saveDir = get(handles.edit_outDir,'String');
% Test the group setup (groupList and output directory)
if isempty(userData.groupList) || isempty(saveDir)
    warndlg('Select a group and an output directory first');
    return
end

% Read common value for statistical tests
stringency=str2double(get(handles.edit_stringency,'String'));
testValues = get(handles.popupmenu_testID1,'UserData');
testID1=testValues(get(handles.popupmenu_testID1,'Value'));
testID2=testValues(get(handles.popupmenu_testID2,'Value'));

doPlot=get(handles.checkbox_doPlot,'Value');
if get(handles.radiobutton_poolData,'Value')
    doWtn=get(handles.checkbox_doWtn,'Value');
    plusTipPoolGroupData(userData.groupData,...
        [saveDir filesep 'pooledData'],doWtn,doPlot);
    plusTipTestDistrib(userData.groupData,[saveDir filesep 'pooledData'],...
    stringency,testID1,testID2);
end

if get(handles.radiobutton_perCell,'Value')
    plusTipGetHits(userData.groupData,[saveDir filesep 'perCell'],...
    stringency,testID1,testID2);
end

% --- Executes when selected object is changed in uipanel_analysisMode.
function uipanel_analysisMode_SelectionChangeFcn(hObject, eventdata, handles)

if get(handles.radiobutton_poolData,'Value')
    set(handles.checkbox_doWtn,'Enable','on');
    set(handles.popupmenu_testID1,'Value',6);
    set(handles.popupmenu_testID2,'Value',7);
else
    set(handles.checkbox_doWtn,'Enable','off');
    set(handles.popupmenu_testID1,'Value',1);
    set(handles.popupmenu_testID2,'Value',6);    

end
