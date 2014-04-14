function varargout = plusTipGroupAnalysisGUI(varargin)
% plusTipGroupAnalysisGUI M-file for plusTipGroupAnalysisGUI.fig
%      plusTipGroupAnalysisGUI, by itself, creates a new plusTipGroupAnalysisGUI or raises the existing
%      singleton*.
%
%      H = plusTipGroupAnalysisGUI returns the handle to a new plusTipGroupAnalysisGUI or the handle to
%      the existing singleton*.
%
%      plusTipGroupAnalysisGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in plusTipGroupAnalysisGUI.M with the given input arguments.
%
%      plusTipGroupAnalysisGUI('Property','Value',...) creates a new plusTipGroupAnalysisGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipGroupAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipGroupAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipGroupAnalysisGUI

% Last Modified by GUIDE v2.5 20-Mar-2014 21:02:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plusTipGroupAnalysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @plusTipGroupAnalysisGUI_OutputFcn, ...
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


% --- Executes just before plusTipGroupAnalysisGUI is made visible.
function plusTipGroupAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.text_copyright, 'String', getLCCBCopyright())

userData = get(handles.figure1, 'UserData');
% Get main figure handle and process id
ip = inputParser;
ip.addParamValue('mainFig', [], @isscalar)
ip.parse(varargin{:});
userData.mainFig = ip.Results.mainFig;
userData.handles_main = guidata(userData.mainFig);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.ML = userData_main.ML;

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% Set test values
testList = {'t-test of the means';'Wilcoxon ranksum test';'Kolmogorov-Smirnov test (K-S test)';...
    'Mean substracted K-S test';'Median substracted K-S test';...
    'Permutation t-test of the means';'Calibrated mean subtracted K-S test'};
testValues=[1 2 10 11 12 20 21];
set(handles.popupmenu_testID1, 'String', testList, 'UserData', testValues);
set(handles.popupmenu_testID2, 'String', testList, 'UserData', testValues);
uipanel_analysisMode_SelectionChangeFcn(hObject, eventdata, handles)

handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = plusTipGroupAnalysisGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

userData = get(handles.figure1, 'UserData');
if isempty(userData.ML)
    warndlg('At least one movie list is required to perform group analysis.',...
        'Input error', 'modal');
    close(handles.figure1);
    return
end

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
    delete(userData.helpFig)
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)


userData = get(handles.figure1,'UserData');

% Load group data
remBegEnd = get(handles.checkbox_remBegEnd,'Value');
userData.groupData = plusTipExtractGroupData(userData.ML, remBegEnd);


% Read common value for statistical tests
alpha =str2double(get(handles.edit_alpha, 'String'));
testValues = get(handles.popupmenu_testID1,'UserData');
testID1 = testValues(get(handles.popupmenu_testID1,'Value'));
testID2 = testValues(get(handles.popupmenu_testID2,'Value'));

% Run within group comparison
if get(handles.checkbox_doWtn,'Value')
    for i = 1 : numel(userData.ML)
        outputDir = fullfile(userData.ML(i).outputDirectory_,...
            'withinGroupComparison');
        plusTipWithinGroupComparison(userData.groupData, i, outputDir, 1);
    end
end

% Run per-cell group analysis
if numel(userData.ML) > 1 && get(handles.radiobutton_poolData,'Value')
    outputDir = fullfile(userData.ML(1).outputDirectory_, 'pooledData');
    plusTipPoolGroupData(userData.groupData, outputDir, 0, 1);
    
    % Perform statistical tests if more than one list is passed
    plusTipTestDistrib(userData.groupData, outputDir,...
        alpha, testID1, testID2);
end

% Run pooled group analysis
if get(handles.radiobutton_perCell,'Value')
    outputDir = fullfile(userData.ML(1).outputDirectory_, 'perCell');
    plusTipGetHits(userData.groupData, outputDir, alpha, testID1, testID2);
end

arrayfun(@save, userData.ML)