function varargout = fsmCenterProjectStatus(varargin)
% FSMCENTERPROJECTSTATUS M-file for fsmCenterProjectStatus.fig
%      FSMCENTERPROJECTSTATUS, by itself, creates a new FSMCENTERPROJECTSTATUS or raises the existing
%      singleton*.
%
%      H = FSMCENTERPROJECTSTATUS returns the handle to a new FSMCENTERPROJECTSTATUS or the handle to
%      the existing singleton*.
%
%      FSMCENTERPROJECTSTATUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMCENTERPROJECTSTATUS.M with the given input arguments.
%
%      FSMCENTERPROJECTSTATUS('Property','Value',...) creates a new FSMCENTERPROJECTSTATUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fsmCenterProjectStatus_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fsmCenterProjectStatus_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fsmCenterProjectStatus

% Last Modified by GUIDE v2.5 19-Mar-2009 18:27:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fsmCenterProjectStatus_OpeningFcn, ...
                   'gui_OutputFcn',  @fsmCenterProjectStatus_OutputFcn, ...
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


% --- Executes just before fsmCenterProjectStatus is made visible.
function fsmCenterProjectStatus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fsmCenterProjectStatus (see VARARGIN)

% Choose default command line output for fsmCenterProjectStatus
handles.output = hObject;

%
% Check whether the active project is valid
%

fsmH = findall(0, 'Tag', 'fsmCenter');

if ~ishandle(fsmH)
    error('Unable to find fsmCenter window.');
end

settings = get(fsmH,'UserData');

% Reminder how sub projects names are stored:
% subProjects: {'tack'  'lpla'  'post'  'edge'  'merg'  'fadh'  'corr'  'mech'}

projDir = settings.projDir;
tackDir = [projDir filesep settings.subProjects{1}];
postDir = [projDir filesep settings.subProjects{3}];
edgeDir = [projDir filesep settings.subProjects{4}];
corrDir = [projDir filesep settings.subProjects{7}];

if ~exist(projDir, 'dir')
    error('Unable to find the project root directory.');
end

if ~exist(tackDir, 'dir')
    error('Unable to find the ''tack'' directory.');
end

if ~exist(postDir, 'dir')
    error('Unable to find the ''post'' directory.');
end

if ~exist(edgeDir, 'dir')
    error('Unable to find the ''edge'' directory.');
end

if ~exist(corrDir, 'dir')
    error('Unable to find the ''corr'' directory.');
end

%
% Set status
%

% Edge Tracker
path = [edgeDir filesep 'cell_mask'];
status = 2;
if ~(exist(path, 'dir') && ...
     size(dir([path, filesep, '*.tif']), 2))
 status = 0;
end
setElementStatus(handles.textEdgeTracker, status);

% Kymograph Analysis
status = 2;
flowPath = [corrDir, filesep, 'flow'];
if ~(size(dir([corrDir, filesep, 'flowTrack*.mat']), 1) && ...
     exist(flowPath, 'dir') && ...
     size(dir(flowPath, 'flow*.mat'), 1))
 status  = 0;
end
setElementStatus(handles.textKymographAnalysis, status);

% Preprocessing Module
status = 2;
candsPath = [tackDir filesep 'cands'];
locMaxPath = [tackDir filesep 'locMax'];
if ~(exist(candsPath, 'dir') && ...
     size(dir([candsPath, filesep, 'cands*.mat']), 1) && ...
     exist(locMaxPath, 'dir') && ...
     size(dir([locMaxPath, filesep, 'locMax*.mat']), 1))
 status = 0;
end
setElementStatus(handles.textPreprocessingModule, status);

% Tracker Module
status = 2;
gapListPath = [tackDir filesep 'gapList'];     
if ~(exist('mpm.mat', 'file') && ...
     exist(gapListPath, 'dir') && ...
     size(dir([gapListPath, filesep, 'gapList*.mat']), 1))
 status = 0;
end
setElementStatus(handles.textTrackerModule, status);

% Builder Module
status = 2;
if ~exist('speckleArray.mat', 'file')
    status = 0;
end
setElementStatus(handles.textBuilderModule, status);

% Kinetics Analysis Module
status = 2;
kinScorPath = [tackDir filesep 'kinScor'];
if ~(exist('SCORE.mat', 'file') && ...
     exist(kinScorePath, 'dir') && ...
     size(dir([kinScorPath, filesep, 'kinScore*.mat']), 1))
 status = 0;
end
setElementStatus(handles.textKineticsAnalysisModule, status);

% Result Display Module
status = 2;
dispPath = [tackDir filesep 'movies'];
if ~(exist(dispPath, 'dir') && ...
     size(dir([dispPath, filesep, 'score*.jpeg']), 1))
 status = 0;
end
setElementStatus(handles.textResultDisplayModule, status);

% Post process Flow Map
status = 2;
if ~size(dir([postDir, filesep, 'flowMap_d0*.tif']), 1)
    status = 0;
end
setElementStatus(handles.textFlowMap, status);

postTifPath = [postDir, filesep, 'tif'];
postMatPath = [postDir, filesep, 'mat'];
postEpsPath = [postDir, filesep, 'eps'];

if ~(exist(postTifPath, 'dir') && ...
     exist(postMatPath, 'dir') && ...
     exist(postEpsPath, 'dir'))
 setElementStatus(handles.textSpeedMap, 0);
 setElementStatus(handles.textKineticsMap, 0);
else
    % Post speed map
    status = 2;
    if ~(size(dir([postTifPath, filesep, 'speedMap_d0*.tif']), 1) && ...
         size(dir([postMatPath, filesep, 'speedMap_d0*.mat']), 1) && ...
         size(dir([postEpsPath, filesep, 'speedMap_d0*.eps']), 1))
     status = 0;
    end
    setElementStatus(handles.textSpeedMap, status);
    
    % Kinetics map
    status = 2;
    if ~(size(dir([postTifPath, filesep, 'kinMap2C*.tif']), 1) && ...
         size(dir([postMatPath, filesep, 'kinMap2C*.mat']), 1) && ...
         size(dir([postMatPath, filesep, 'depolyMap*.mat']), 1) && ...
         size(dir([postMatPath, filesep, 'polyMap*.mat']), 1))
     status = 0;
    end
    setElementStatus(handles.textKineticsMap, status);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fsmCenterProjectStatus wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fsmCenterProjectStatus_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushOK.
function pushOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close the status window.
close(get(hObject, 'Parent'));

function setElementStatus(hElement, status)
colorStatus = zeros(1, 3);
stringStatus = '';

switch status
    case 0, colorStatus = [0.847, 0.161, 0.0]; stringStatus = 'Not done';
    case 1, colorStatus = [1.0, 1.0, 0]; stringStatus = 'Problem';
    case 2, colorStatus = [0, 0.49, 0]; stringStatus = 'Done';
end

set(hElement, 'String', stringStatus);
set(hElement, 'ForegroundColor', colorStatus);

