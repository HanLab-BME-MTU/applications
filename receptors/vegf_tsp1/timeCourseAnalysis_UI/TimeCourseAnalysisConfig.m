function varargout = TimeCourseAnalysisConfig(varargin)
% TIMECOURSEANALYSISCONFIG MATLAB code for TimeCourseAnalysisConfig.fig
%      TIMECOURSEANALYSISCONFIG, by itself, creates a new TIMECOURSEANALYSISCONFIG or raises the existing
%      singleton*.
%
%      H = TIMECOURSEANALYSISCONFIG returns the handle to a new TIMECOURSEANALYSISCONFIG or the handle to
%      the existing singleton*.
%
%      TIMECOURSEANALYSISCONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIMECOURSEANALYSISCONFIG.M with the given input arguments.
%
%      TIMECOURSEANALYSISCONFIG('Property','Value',...) creates a new TIMECOURSEANALYSISCONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TimeCourseAnalysisConfig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TimeCourseAnalysisConfig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TimeCourseAnalysisConfig

% Last Modified by GUIDE v2.5 26-Jan-2016 11:42:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TimeCourseAnalysisConfig_OpeningFcn, ...
                   'gui_OutputFcn',  @TimeCourseAnalysisConfig_OutputFcn, ...
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


% --- Executes just before TimeCourseAnalysisConfig is made visible.
function TimeCourseAnalysisConfig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TimeCourseAnalysisConfig (see VARARGIN)

% Choose default command line output for TimeCourseAnalysisConfig
handles.output = hObject;
handles.p = struct();
skipUI = false;
nIn = nargin;

if(nIn < 4 || isempty(varargin{1}))
    %% Prompt user
    %prompt user to select a folder where all figures and data will be stored
    handles.p.outputDir = uigetdir('', 'Select output folder');
elseif(isstruct(varargin{1}))
    handles.p = varargin{1};
    varargin{2} = handles.p.CML_FullPath;
    varargin{3} = handles.p.channelTable;
    nIn = 6;
    skipUI = true;
else
    handles.p.outputDir = varargin{1}{1};
end

if(nIn < 5 || isempty(varargin{2}))
    %prompt user to select Combined Movie Data objects
    %until they press cancel.
    handles.p.CML_FullPath = {};
    [fileName, filePath] = uigetfile('*.mat', ...
        'Select CombinedMovieLists', 'MultiSelect', 'on');
    while ~isnumeric(fileName)
        if iscell(fileName)
            handles.p.CML_FullPath = [handles.p.CML_FullPath ...
                cellfun(@(x) [filePath x], fileName, ...
                'UniformOutput', false)]; %#ok<*AGROW>
        else
            handles.p.CML_FullPath{end+1} = [filePath fileName];
        end
        [fileName, filePath] = uigetfile('*.mat', ...
            'Select CombinedMovieLists', 'MultiSelect', 'on');
    end
else
    handles.p.CML_FullPath = varargin{2};
end

if(nIn < 6 || isempty(varargin{3}))
    % Try to get some basic channel data
    try
        S = load(handles.p.CML_FullPath{1});
        S = load(S.CML.movieListDirectory_{1});
        S = load(S.ML.movieDataFile_{1});
        c = S.MD.channels_;
        data(:,3:4) = [{c.name_}' {c.emissionWavelength_}'];
        data(:,2) = strcat({'Channel '},cellstr(num2str((1:size(data,1))')));
        data(:,1) = num2cell(true(size(data,1),1));
        set(handles.channelTable,'data',data);
    catch err
        disp('Could not obtain channel data');
        disp(getReport(err));
    end
else
    set(handles.channelTable,'data',varargin{3});
end

% Retrieve list of available cluster profiles and append to menu
profiles = parallel.clusterProfiles;
set(handles.batchmenu,'String',...
    [get(handles.batchmenu,'String') ; profiles']);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes TimeCourseAnalysisConfig wait for user response (see UIRESUME)
if(~skipUI)
    uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = TimeCourseAnalysisConfig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if(isstruct(handles))
    varargout{1} = handles.p;
    close(handles.figure1);
    pause(1);
else
    % Canceled
    warning('TimeCourseAnalysisConfig canceled');
    varargout{1} = [];
end


% --- Executes on button press in start2zero.
function start2zero_Callback(hObject, eventdata, handles)
% hObject    handle to start2zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of start2zero


% --- Executes on button press in doPartition.
function doPartition_Callback(hObject, eventdata, handles)
% hObject    handle to doPartition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doPartition


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%closes the dialogue box
handles.p.doNewAnalysis = get(handles.doNewAnalysis,'Value');
handles.p.doPartition = get(handles.doPartition,'Value');
handles.p.start2zero = get(handles.start2zero,'Value');
handles.p.shiftPlotPositive = get(handles.shiftPlotPositive,'Value');
handles.p.nBootstrp = get(handles.nBootstrp,'Value');
handles.p.channelTable = get(handles.channelTable,'Data');
handles.p.detectOutliers_k_sigma = get(handles.detectOutliers,'Value') ...
                                   .*str2double(get(handles.k_sigma,'String'));
handles.p.ignoreIsolatedPts = get(handles.isolatedPts,'Value');
handles.p.channels = find([handles.p.channelTable{:,1}]);
handles.p.channelNames = handles.p.channelTable(:,2);
handles.p.showPlots = get(handles.showPlots,'Value');
menuOptions = get(handles.batchmenu,'String');
handles.p.batchClusterName = menuOptions{get(handles.batchmenu,'Value')};
guidata(handles.figure1, handles);
uiresume(handles.figure1);
% Go to TimeCourseAnalysisConfig_OutputFcn


% --- Executes on button press in shiftPlotPositive.
function shiftPlotPositive_Callback(hObject, eventdata, handles)
% hObject    handle to shiftPlotPositive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shiftPlotPositive


% --- Executes on button press in nBootstrp.
function nBootstrp_Callback(hObject, eventdata, handles)
% hObject    handle to nBootstrp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nBootstrp


% --- Executes on button press in detectOutliers.
function detectOutliers_Callback(hObject, eventdata, handles)
% hObject    handle to detectOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detectOutliers
detectOutliers = get(hObject, 'Value') == get(hObject,'Max');
visibleStates = {'off', 'on'};
set([handles.k_sigma_label handles.k_sigma], ...
    'Visible',visibleStates{detectOutliers + 1});


function k_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to k_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_sigma as text
%        str2double(get(hObject,'String')) returns contents of k_sigma as a double


% --- Executes on button press in isolatedPts.
function isolatedPts_Callback(hObject, eventdata, handles)
% hObject    handle to isolatedPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isolatedPts


% --- Executes during object creation, after setting all properties.
function k_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showPlots.
function showPlots_Callback(hObject, eventdata, handles)
% hObject    handle to showPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPlots


% --- Executes on selection change in batchmenu.
function batchmenu_Callback(hObject, eventdata, handles)
% hObject    handle to batchmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns batchmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from batchmenu


% --- Executes during object creation, after setting all properties.
function batchmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batchmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
