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

% Last Modified by GUIDE v2.5 27-Oct-2011 14:42:59


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
handles.projData=[]; % if one project is selected, projData will be retrieved
handles.strList='';

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

% for "sub-ROIs" panel
set(handles.subroiDistUnitPop,'String',{'Microns','Fraction'},'Value',1);
set(handles.subroiTimeUnitPop,'String',{'Fraction','Seconds'},'Value',1);

% for "quadrant scatter plots" panel
scatterData = {'growthSpeed','growthLifetime','growthDisp','bgapSpeed',...
    'bgapLifetime','bgapDisp''fgapSpeed','fgapLifetime','fgapDisp'};
scatterString={'Growth speed (microns/min)','Growth lifetime (s)',...
    'Growth displacement (microns)','Bgap speed (microns/min)','Bgap lifetime (s)',...
    'Bgap displacement (microns)','Fgap speed (microns/min)','Fgap lifetime (s)',...
    'Fgap displacement (microns)'};
set(handles.xaxisScatterDrop,'String',scatterString','UserData',scatterData,'Value',1);
set(handles.yaxisScatterDrop,'String',scatterString','UserData',scatterData,'Value',2);
set([handles.xParamDrop handles.yParamDrop],'String',{'Value','Percentile'},'Value',1);

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
        assignin('base','projData',handles.projData)
    else
        handles.projList=handles.projList(selection,1);
    end
else
    msgbox('No projects selected.')
    handles.dataDir=[];
    return
end

% Allow the user to save the project list if not from a single projList
isSingleProjList = handles.nProjLists == 1 && ~get(handles.getQueryStr_Check,'Value');
if ~isempty(handles.projList) && ~isSingleProjList
    if ~isempty(handles.strList)
        defaultListName = ['projList', sprintf('_%s',handles.strList{:}) '.mat'];
    else
        defaultListName='projList.mat';
    end
    [file,path] = uiputfile('projList.mat','Find a location to save your project list',...
        [pwd filesep defaultListName]);
    if ~isequal(file,0), 
        save([path file],'-struct','handles','projList');
    end
end

% If checked, create groups from selected projects
if get(handles.checkbox_createGroups,'Value')
    autoGrp =get(handles.checkbox_autoGrp,'Value');
    userData.groupList=plusTipPickGroups(autoGrp,[],handles.projList,1);
    set(handles.figure1,'UserData',userData);
end
guidata(hObject, handles);

% --- Executes on button press in selectOutputDirPush.
function selectOutputDirPush_Callback(hObject, eventdata, handles)

if isfield(handles,'dataDir')
    dirStart=handles.dataDir;
else
    dirStart=pwd;
end

outputDir=uigetdir(dirStart,'Please select output directory for group analyis');

if isequal(outputDir,0), return; end

userData = get(handles.figure1,'UserData');
set(handles.edit_outputDir,'String',outputDir);
set(handles.figure1,'UserData',userData);


% --- Executes on button press in pushbutton_loadGroup.
function pushbutton_loadGroup_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
[userData.groupList]=combineGroupListFiles(0);
set(handles.figure1,'UserData',userData);

% --- Executes on button press in pushbutton_analyzeGroups.
function pushbutton_analyzeGroups_Callback(hObject, eventdata, handles)

userData = get(handles.figure1,'UserData');
saveDir = get(handles.edit_outputDir,'String');
% Test the group setup (groupList and output directory)
if isempty(userData.groupList) || isempty(saveDir)
    warndlg('Select a group and an output directory first');
    return
end

% Load group data
remBegEnd = get(handles.checkbox_remBegEnd,'Value');
userData.groupData=plusTipExtractGroupData(userData.groupList,remBegEnd);

% Read common value for statistical tests
stringency=str2double(get(handles.edit_stringency,'String'));
testValues = get(handles.popupmenu_testID1,'UserData');
testID1=testValues(get(handles.popupmenu_testID1,'Value'));
testID2=testValues(get(handles.popupmenu_testID2,'Value'));

doPlot=1;
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

% --- Executes when selected object is changed in subroiRadioPanel.
function subroiRadioPanel_SelectionChangeFcn(hObject, eventdata, handles)

manualSelection = strcmp(get(hObject,'Tag'),'subroiManualRadio');
if manualSelection, state = 'off'; else state='on'; end
set(handles.subroiDistValEdit,'Enable',state);
set(handles.subroiDistUnitPop,'Enable',state);
set(handles.subroiAutoDivPeriphCheck,'Enable',state);
if manualSelection,set(handles.subroiDistValEdit,'String',''); end

% --- Executes on button press in subRoiPush.
function subRoiPush_Callback(hObject, eventdata, handles)

% Read selection type
if get(handles.subroiManualRadio,'Value'),
    subroiSelectType =0;
elseif get(handles.subroiAutoDivPeriphCheck,'Value'),
    subroiSelectType =2;
else
    subroiSelectType =1;
end

% Read distance units and value
if subroiSelectType>0
    props = get(handles.subroiDistUnitPop,{'String','Value'});
    subroiDistUnit=props{1}{props{2}};
    subroiDistVal=str2double(get(handles.subroiDistValEdit,'String'));
    if strcmpi(subroiDistUnit,'fraction');
        check=@(x) x>0 && x<=1 && ~isnan(x);
    else
        check= @(x) x>0 && ~isnan(x);
    end
    if ~check(subroiDistVal)
        errordlg('Please enter a valid value for the distance from cell edge');
        return;
    end
else
    subroiDistVal=[];
    subroiDistUnit='';
end

% Read time units and value
props = get(handles.subroiTimeUnitPop,{'String','Value'});
subroiTimeUnit=props{1}{props{2}};
subroiTimeVal=str2double(get(handles.subroiTimeValEdit,'String'));
if strcmpi(subroiTimeUnit,'fraction');
    check=@(x) x>0 && x<=1 && ~isnan(x);
else
    check= @(x) x>0 && ~isnan(x);
end

if ~check(subroiTimeVal)
    errordlg('Please enter a valid value for the time in the sub-region');
    return;
end


subroiExcludeRegion = get(handles.subroiExcludeCheck,'Value');
useSegMask = get(handles.checkbox_useSegMask,'Value');

plusTipSubRoiTool(handles.projList,subroiSelectType,...
    subroiDistUnit,subroiDistVal,subroiTimeUnit,subroiTimeVal,...
    [],subroiExcludeRegion,useSegMask);


% --- Executes on button press in quadScatterPlotPush.
function quadScatterPlotPush_Callback(hObject, eventdata, handles)

% Read groupList and output directory
userData = get(handles.figure1,'UserData');
saveDir=get(handles.edit_outputDir,'String');
if isempty(userData.groupList) || isempty(saveDir)
    warndlg('Select a group and an output directory first');
    return
end

% Read scatter plot parameters
xAxisInfo=getAxisParameters(handles,'x');
yAxisInfo=getAxisParameters(handles,'y');
timeRange=getTimeRange(handles);
doPlotQuad=~get(handles.summQuadPlotOnlyCheck,'Value');
remBegEnd = get(handles.remTrackBegEnd,'Value');

plusTipQuadScatter(xAxisInfo,yAxisInfo,userData.groupList,remBegEnd,...
    timeRange,doPlotQuad,saveDir);


function timeRange = getTimeRange(handles)
% Read post-processing time range
sFVal=get(handles.startFrame,'String');
if strcmpi(sFVal,'min')
    timeRange(1)=1;
else
    timeRange(1)=str2double(sFVal);
end
eFVal=get(handles.endFrame,'String');
if strcmpi(eFVal,'max')
    timeRange(2)=Inf;
else
    timeRange(2)=str2double(eFVal);
end
if ~all(isposint(timeRange)),
    errordlg('Invalid time range. Please check again','Missing Input');
    return
end


function axisInfo = getAxisParameters(handles,type)
% Read axis scatter name
props=get(handles.([type 'axisScatterDrop']),{'UserData','Value'});
axisInfo.name=props{1}{props{2}};

% Read axis scatter input
value = str2double(get(handles.([type 'AxisScatterInput']),'String'));
if isnan(value),
    errordlg(['Please enter a valid value for the '  type '-axis scatter']);
    return;
end
props = get(handles.([type 'ParamDrop']),{'String','Value'});
if strcmpi(props{1}{props{2}},'percentile')
    axisInfo.splitPercentile=value;
    axisInfo.splitValue=[];
else
    axisInfo.splitPercentile=[];
    axisInfo.splitValue=value;
end

% Read axis limits
value=get(handles.([type 'AxisMin']),'String');
if strcmpi(value,'min'), axisInfo.minMax(1) = -Inf;
else axisInfo.minMax(1) = str2double(value);
end

value=get(handles.([type 'AxisMax']),'String');
if strcmpi(value,'max'), axisInfo.minMax(2) = Inf;
else axisInfo.minMax(2) = str2double(value);
end