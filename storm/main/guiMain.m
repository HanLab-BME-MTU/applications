function varargout = guiMain(varargin)
% GUIMAIN MATLAB code for guiMain.fig
%      GUIMAIN, by itself, creates a new GUIMAIN or raises the existing
%      singleton*.
%
%      H = GUIMAIN returns the handle to a new GUIMAIN or the handle to
%      the existing singleton*.
%
%      GUIMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIMAIN.M with the given input arguments.
%
%      GUIMAIN('Property','Value',...) creates a new GUIMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUIMAIN before guiMain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiMain_OpeningFcn via varargin.
%
%      *See GUIMAIN Options on GUIDE's Tools menu.  Choose "GUIMAIN allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiMain

% Last Modified by GUIDE v2.5 10-Feb-2012 11:24:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @guiMain_OpeningFcn, ...
    'gui_OutputFcn',  @guiMain_OutputFcn, ...
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

% --- Executes just before guiMain is made visible.
function guiMain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiMain (see VARARGIN)

% Choose default command line output for guiMain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiMain wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Set the close request function
set(gcf,'CloseRequestFcn',{@closeFcn,handles})
% Load the last GUIMAIN state
loadState(handles);
% Init the data path
data = guidata(gcf);
data.rootPath = getStormPath();
guidata(gcf,data);
% Detect the data sets
updateDataSets(handles);
% Detect the presets
updatePresets(handles);
% Detect queue items
updateQueue(handles);
% Autoload the current configuration
data = guidata(gcf);
if exist(data.presets{get(handles.listbox2,'Value'),1},'file')
    updateCurrentConfig(data.presets{get(handles.listbox2,'Value'),1},handles);
end


% --- Outputs from this function are returned to the command line.
function varargout = guiMain_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function closeFcn(src,event,handles)
saveState(handles);
clc;
delete(src);


function saveState(handles)
state.listbox1Value = get(handles.listbox1,'value');
state.listbox2Value = get(handles.listbox2,'value');
state.listbox3Value = get(handles.listbox3,'value');
fileName = [mfilename('fullpath') '.mat'];
save(fileName,'state');


function loadState(handles)
fileName = [mfilename('fullpath') '.mat'];
if exist(fileName,'file')
    load(fileName);
    set(handles.listbox1,'value',state.listbox1Value);
    set(handles.listbox2,'value',state.listbox2Value);
    set(handles.listbox3,'value',state.listbox3Value);
    delete(fileName);
end


function updateDataSets(handles)
data = guidata(gcf);
path = [data.rootPath '_data\'];
list = dir(path);
dataSetIdx = 0;
data.dataSets = [];
for i=3:numel(list)
    if list(i).isdir
        sublist = dir([path list(i).name '\']);
        for k=3:numel(sublist)
            if strcmp(sublist(k).name(end-3:end),'.bin')
                dataSetIdx = dataSetIdx + 1;
                data.dataSets{dataSetIdx,1} = [path list(i).name '\' sublist(k).name];
                data.dataSets{dataSetIdx,2} = list(i).name;
                break;
            end
        end
    end
end
% Check if the listbox value is valid
if get(handles.listbox1,'Value') > size(data.dataSets,1)
    set(handles.listbox1,'Value',1);
end
% Setup data set listbox
set(handles.listbox1,'String',data.dataSets(:,2));
guidata(gcf,data);


function updatePresets(handles)
data = guidata(gcf);
path = [data.rootPath '_presets\'];
list = dir(path);
presetIdx = 0;
data.presets = [];
for i=3:numel(list)
    if strcmp(list(i).name(end-3:end),'.cfg')
        presetIdx = presetIdx + 1;
        data.presets{presetIdx,1} = [path list(i).name];
        data.presets{presetIdx,2} = list(i).name(1:end-4);
    end
end
% Check if the listbox value is valid
if get(handles.listbox2,'Value') > size(data.presets,1)
    set(handles.listbox2,'Value',1);
end
% Setup presets listbox
set(handles.listbox2,'String',data.presets(:,2));
guidata(gcf,data);


function updateQueue(handles)
data = guidata(gcf);
path = [data.rootPath '_queue\'];
list = dir(path);
itemIdx = 0;
data.queue = [];
for i=3:numel(list)
    if list(i).isdir
        sublist = dir([path list(i).name '\']);
        sublist = sublist(~vertcat(sublist.isdir));
        for k=1:numel(sublist)
            if strcmp(sublist(k).name(end-3:end),'.cfg')
                itemIdx = itemIdx + 1;
                data.queue{itemIdx,1} = [path list(i).name '\'];
                data.queue{itemIdx,2} = list(i).name;
                data.queue{itemIdx,3} = [sublist(k).name];
                data.queue{itemIdx,4} = 0;
                for j=1:numel(sublist)
                    if numel(sublist(j).name) > 5
                        switch (sublist(j).name(end-5:end))
                            case '.d.dat'
                                data.queue{itemIdx,4} = data.queue{itemIdx,4} + bin2dec('100');
                                data.queue{itemIdx,5} = [path list(i).name '\' sublist(j).name];
                            case '.i.dat'
                                data.queue{itemIdx,4} = data.queue{itemIdx,4} + bin2dec('010');
                                data.queue{itemIdx,5} = [path list(i).name '\' sublist(j).name];
                            case '.p.dat'
                                data.queue{itemIdx,4} = data.queue{itemIdx,4} + bin2dec('001');
                                data.queue{itemIdx,5} = [path list(i).name '\' sublist(j).name];
                        end
                    end
                end
                break;
            end
        end
    end
end
% Check if the listbox value is valid
if get(handles.listbox3,'Value') > size(data.queue,1)
    set(handles.listbox3,'Value',1);
end
% Setup queue listbox
for i=1:size(data.queue,1)
    status = dec2bin(data.queue{i,4},3);
    status = strrep(status,'1','/');
    status = strrep(status,'0',' ');
    if status(2) == '/'
    status(2) = '\';
    end
    dispNames{i} = [status '   ' data.queue{i,2}];
end
set(handles.listbox3,'String',dispNames);
guidata(gcf,data);


function updateCurrentConfig(fullPath,handles)
data = guidata(gcf);
data.currentConfig = Config.load(fullPath);
data.configLoadedFrom = fullPath;
set(handles.edit2,'String',data.currentConfig.configName);
guidata(gcf,data);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb1');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb2');


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb4');
data = guidata(gcf);
fullPathName = data.dataSets{get(handles.listbox1,'Value'),1};
delimiterPos = strfind(fullPathName,'\');
pathName = fullPathName(1:delimiterPos(end));
fileName = fullPathName(delimiterPos(end)+1:end);
myROISelector = ROISelector(pathName,fileName);
myROISelector.roiPos = data.currentConfig.roiPosition;
myROISelector.roiSize = data.currentConfig.roiSize;
myROISelector.run();
data.currentConfig.roiPosition = myROISelector.roiPos;
data.currentConfig.roiSize = myROISelector.roiSize;
guidata(gcf,data);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb5');
data = guidata(gcf);
data.currentConfig.save([data.rootPath '_presets\' data.currentConfig.configName '.cfg']);
% Detect the presets
updatePresets(handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb6');
[fileName,pathName] = uigetfile('*.cfg','Select the configuration file');
if pathName ~= 0
    updateCurrentConfig([pathName fileName],handles);
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb9');
data = guidata(gcf);
updateCurrentConfig(data.presets{get(handles.listbox2,'Value'),1},handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb10');
runAlgorithm();


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb11');
data = guidata(gcf);
dataDirPath = data.queue{get(handles.listbox3,'Value'),1};
configName = data.queue{get(handles.listbox3,'Value'),3};
dataPath = [dataDirPath strtok(configName,'.') '.p.dat'];
if exist(dataPath,'file')
    configPath = [dataDirPath configName];
    assignin('base','configPath',configPath);
    assignin('base','dataPath',dataPath);
    evalin('base','cfg = Config.load(configPath);');
    evalin('base','data = Data.load(dataPath);');
else
    disp('Main: No data found!')
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb13');
data = guidata(gcf);
delete(data.presets{get(handles.listbox2,'Value'),1});
updatePresets(handles);


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb14');
updateQueue(handles);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb15');
dos('C:\Users\PB93\Desktop\Software\putty.exe');


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb16');


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb19'); % Add to Queue
data = guidata(gcf);
fullPath = data.dataSets{get(handles.listbox1,'Value'),1};
delimiterPos = strfind(fullPath,'\');
data.currentConfig.fileName = fullPath(delimiterPos(end)+1:end);
data.currentConfig.path = fullPath(1:delimiterPos(end));
dataSetName = data.dataSets{get(handles.listbox1,'Value'),2};
dirName = [dataSetName '--' datestr(now,'yymmdd--HHMM--SS--') strtok(data.currentConfig.configName,'.')];
mkdir([data.rootPath '_queue\'],dirName);
data.currentConfig.save([data.rootPath '_queue\' dirName '\' data.currentConfig.configName '.cfg']);
guidata(gcf,data);
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb20');


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb21');
updatePresets(handles);


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb22');
updateDataSets(handles);


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb23');
data = guidata(gcf);
prompt = {...
    'dataReductionEnabled',...
    'nReductionRun',...
    'reductionEdgeRadius',...
    'densityFilteringEnabled',...
    'nNeighborsThreshold',...
    'neighborBallRadius',...
    'nIterEM',...
    'maxIterMerge',...
    'nIterGeomMatching',...
    'modelLength',...
    'angleThreshold',...
    'snapshotsPath',...
    'subsampleFraction',...
    };
name = 'Edit Advanced Parameters';
defaultanswer = {...
    num2str(data.currentConfig.dataReductionEnabled),...
    num2str(data.currentConfig.nReductionRun),...
    num2str(data.currentConfig.reductionEdgeRadius),...
    num2str(data.currentConfig.densityFilteringEnabled),...
    num2str(data.currentConfig.nNeighborsThreshold),...
    num2str(data.currentConfig.neighborBallRadius),...
    num2str(data.currentConfig.nIterEM),...
    num2str(data.currentConfig.maxIterMerge),...
    num2str(data.currentConfig.nIterGeomMatching),...
    num2str(data.currentConfig.modelLength),...
    num2str(data.currentConfig.angleThreshold),...
    data.currentConfig.snapshotsPath,...
    num2str(data.currentConfig.subsampleFraction),...
    };
options.Resize = 'on';
options.WindowStyle = 'normal';
answer = inputdlg(prompt,name,[ones(13,1),ones(13,1)*50],defaultanswer,options);
if ~isempty(answer)
    data.currentConfig.dataReductionEnabled  = str2double(answer{1});
    data.currentConfig.nReductionRun  = str2double(answer{2});
    data.currentConfig.reductionEdgeRadius  = str2double(answer{3});
    data.currentConfig.densityFilteringEnabled  = str2double(answer{4});
    data.currentConfig.nNeighborsThreshold  = str2double(answer{5});
    data.currentConfig.neighborBallRadius  = str2double(answer{6});
    data.currentConfig.nIterEM  = str2double(answer{7});
    data.currentConfig.maxIterMerge  = str2double(answer{8});
    data.currentConfig.nIterGeomMatching  = str2double(answer{9});
    data.currentConfig.modelLength  = str2double(answer{10});
    data.currentConfig.angleThreshold  = str2double(answer{11});
    data.currentConfig.snapshotsPath = answer{12};
    data.currentConfig.subsampleFraction = str2double(answer{13});
end
guidata(gcf,data);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb8');
data = guidata(gcf);
prompt = {...
    'Configuration name',...
    'maxDegreeBezier',...
    'maxCurvature',...
    'errorX,errorY,errorZ',...
    'nSigmaThreshold',...
    'roiPosition(1),roiPosition(2),roiPosition(3)',...
    'roiSize(1),roiSize(2),roiSize(3)',...
    'initialEdgeRadiusGeom',...
    'initialEdgeRadius',...
    'filterLength',...
    'angularSampling',...
    'displayEnabled',...
    'snapshotsEnabled',...
    'fitMethod',...
    'edgeWidthInitFree',...
    'betaVar',...
    'modeVar'...
    };
name = 'Edit Parameters';
defaultanswer = {...
    data.currentConfig.configName,...
    num2str(data.currentConfig.maxDegreeBezier),...
    num2str(data.currentConfig.maxCurvature),...
    [num2str(data.currentConfig.errorX) ',' num2str(data.currentConfig.errorY) ',' num2str(data.currentConfig.errorZ)],...
    num2str(data.currentConfig.nSigmaThreshold),...
    [num2str(data.currentConfig.roiPosition(1)) ',' num2str(data.currentConfig.roiPosition(2)) ',' num2str(data.currentConfig.roiPosition(3))],...
    [num2str(data.currentConfig.roiSize(1)) ',' num2str(data.currentConfig.roiSize(2)) ',' num2str(data.currentConfig.roiSize(3))],...
    num2str(data.currentConfig.initialEdgeRadiusGeom),...
    num2str(data.currentConfig.initialEdgeRadius),...
    num2str(data.currentConfig.filterLength),...
    num2str(data.currentConfig.angularSampling),...
    num2str(data.currentConfig.displayEnabled),...
    num2str(data.currentConfig.snapshotsEnabled),...
    num2str(data.currentConfig.fitMethod),...
    num2str(data.currentConfig.edgeWidthInitFree),...
    num2str(data.currentConfig.betaVar),...
    num2str(data.currentConfig.modeVar)...
    };
options.Resize = 'on';
options.WindowStyle = 'normal';
answer = inputdlg(prompt,name,[ones(17,1),ones(17,1)*50],defaultanswer,options);
if ~isempty(answer)
    data.currentConfig.configName = answer{1};
    set(handles.edit2,'String',data.currentConfig.configName);
    data.currentConfig.maxDegreeBezier = str2double(answer{2});
    data.currentConfig.maxCurvature = str2double(answer{3});
    delimiterPos = strfind(answer{4},',');
    data.currentConfig.errorX = str2double(answer{4}(1:delimiterPos(1)-1));
    data.currentConfig.errorY = str2double(answer{4}(delimiterPos(1)+1:delimiterPos(2)-1));
    data.currentConfig.errorZ = str2double(answer{4}(delimiterPos(2)+1:end));
    data.currentConfig.nSigmaThreshold = str2double(answer{5});
    delimiterPos = strfind(answer{6},',');
    data.currentConfig.roiPosition(1) = str2double(answer{6}(1:delimiterPos(1)-1));
    data.currentConfig.roiPosition(2) = str2double(answer{6}(delimiterPos(1)+1:delimiterPos(2)-1));
    data.currentConfig.roiPosition(3) = str2double(answer{6}(delimiterPos(2)+1:end));
    delimiterPos = strfind(answer{7},',');
    data.currentConfig.roiSize(1) = str2double(answer{7}(1:delimiterPos(1)-1));
    data.currentConfig.roiSize(2) = str2double(answer{7}(delimiterPos(1)+1:delimiterPos(2)-1));
    data.currentConfig.roiSize(3) = str2double(answer{7}(delimiterPos(2)+1:end));
    data.currentConfig.initialEdgeRadiusGeom = str2double(answer{8});
    data.currentConfig.initialEdgeRadius = str2double(answer{9});
    data.currentConfig.filterLength = str2double(answer{10});
    data.currentConfig.angularSampling = str2double(answer{11});
    data.currentConfig.displayEnabled = str2double(answer{12});
    data.currentConfig.snapshotsEnabled = str2double(answer{13});
    data.currentConfig.fitMethod = str2double(answer{14});
    data.currentConfig.edgeWidthInitFree = str2double(answer{15});
    data.currentConfig.betaVar = str2double(answer{16});
    data.currentConfig.modeVar = str2double(answer{17});
end
guidata(gcf,data);


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb25');
data = guidata(gcf);
% dirPath = data.queue{get(handles.listbox3,'Value'),1};
% rmdir(dirPath,'s');
dirPathSrc = data.queue{get(handles.listbox3,'Value'),1};
dirPathDest = strrep(dirPathSrc,'_queue','_deleted');
movefile(dirPathSrc,dirPathDest,'f');
updateQueue(handles);



% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb26');
data = guidata(gcf);
src = data.queue{get(handles.listbox3,'Value'),1};
if strcmp(src(end-1),'+')
    dest = [data.queue{get(handles.listbox3,'Value'),1}(1:end-2) '-\'];
elseif strcmp(src(end-1),'-')
    dest = [data.queue{get(handles.listbox3,'Value'),1}(1:end-2) '\'];
else
    dest = [data.queue{get(handles.listbox3,'Value'),1}(1:end-1) '-\'];
end
movefile(src,dest,'f');
currentPos = get(handles.listbox3,'Value');
if currentPos < size(data.queue,1)
    set(handles.listbox3,'Value',currentPos+1);
end
updateQueue(handles);


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb27');
data = guidata(gcf);
updateCurrentConfig([data.queue{get(handles.listbox3,'Value'),1} data.queue{get(handles.listbox3,'Value'),3}],handles);


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb28');
data = guidata(gcf);
src = data.queue{get(handles.listbox3,'Value'),1};
if ~strcmp(src(end-1),'+')
    for i=1:size(data.queue,1)
        if strcmp(data.queue{i,2}(end),'+')
            src = data.queue{i,1};
            dest = [data.queue{i,1}(1:end-2) '\'];
            movefile(src,dest,'f')
        end
    end
    updateQueue(handles);
    data = guidata(gcf);
    src = data.queue{get(handles.listbox3,'Value'),1};
    if strcmp(src(end-1),'-')
        dest = [data.queue{get(handles.listbox3,'Value'),1}(1:end-2) '+\'];
    else
        dest = [data.queue{get(handles.listbox3,'Value'),1}(1:end-1) '+\'];
    end
    movefile(src,dest,'f');
    updateQueue(handles);
else
    for i=1:size(data.queue,1)
        if strcmp(data.queue{i,2}(end),'+')
            src = data.queue{i,1};
            dest = [data.queue{i,1}(1:end-2) '\'];
            movefile(src,dest,'f')
        end
    end
    updateQueue(handles);
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb30');
data = guidata(gcf);
% Load configuration
fileName = data.queue{get(handles.listbox3,'Value'),3};
path = data.queue{get(handles.listbox3,'Value'),1};
fullPath = [path fileName];
assignin('base', 'fullPath', fullPath);
evalin('base','cfg = Config.load(fullPath);');
if data.queue{get(handles.listbox3,'Value'),4}
    dataPath = data.queue{get(handles.listbox3,'Value'),5};
    assignin('base', 'dataPath', dataPath);
    evalin('base','data = Data.load(dataPath);');
else
    % Read data
    evalin('base','data = Data.read([cfg.path cfg.fileName]);');
    evalin('base','pro = Processor(data);');
    evalin('base','pro.cropData(cfg.roiPosition,cfg.roiSize);');
end
% Launch Imaris
evalin('base','dis = Show(data);');
% Display the data
evalin('base','dis.points();');
disp('Main: Show Data in Imaris: Done!')


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb29');
data = guidata(gcf);
dirPathSrc = data.queue{get(handles.listbox3,'Value'),1};
dirPathDest = strrep(dirPathSrc,'_queue','_archive');
movefile(dirPathSrc,dirPathDest,'f');
updateQueue(handles);


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb31');
data = guidata(gcf);
fullPathName = data.dataSets{get(handles.listbox1,'Value'),1};
delimiterPos = strfind(fullPathName,'\');
pathName = fullPathName(1:delimiterPos(end));
fileName = fullPathName(delimiterPos(end)+1:end);
myROISelector = ROISelector(pathName,fileName);
myROISelector.roiPos = data.currentConfig.roiPosition;
myROISelector.roiSize = data.currentConfig.roiSize;
myROISelector.runWindowing();
% data.currentConfig.roiPosition = myROISelector.roiPos;
% data.currentConfig.roiSize = myROISelector.roiSize;
% guidata(gcf,data);


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb32');
data = guidata(gcf);
fileName = data.queue{get(handles.listbox3,'Value'),3};
path = data.queue{get(handles.listbox3,'Value'),1};
cfg = Config.load([path fileName]);
cfg.display();



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
disp('e2');
data = guidata(gcf);
data.currentConfig.configName = get(hObject,'String');
guidata(gcf,data);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb33'); % Repetition
data = guidata(gcf);
fullPathName = data.dataSets{get(handles.listbox1,'Value'),1};
delimiterPos = strfind(fullPathName,'\');
pathName = fullPathName(1:delimiterPos(end));
fileName = fullPathName(delimiterPos(end)+1:end);
myROISelector = ROISelector(pathName,fileName);
myROISelector = myROISelector.load([pathName fileName]);
myROISelector.roiPos = data.currentConfig.roiPosition;
myROISelector.roiSize = data.currentConfig.roiSize;
state = guiRepetition(myROISelector,data.currentConfig);

if state
    for x=1:myROISelector.repetitions(1)
        for y=1:myROISelector.repetitions(2)
            pos = data.currentConfig.roiPosition(1:2) + ([x y]-1).*myROISelector.offset;
            posBackup = data.currentConfig.roiPosition(1:2);
            data.currentConfig.roiPosition(1:2) = pos;
            
            configNameBackup = data.currentConfig.configName;
            data.currentConfig.configName = sprintf([data.currentConfig.configName '_%02u_%02u'],x,y);
            
            fullPath = data.dataSets{get(handles.listbox1,'Value'),1};
            delimiterPos = strfind(fullPath,'\');
            data.currentConfig.fileName = fullPath(delimiterPos(end)+1:end);
            data.currentConfig.path = fullPath(1:delimiterPos(end));
            dataSetName = data.dataSets{get(handles.listbox1,'Value'),2};
            dirName = [dataSetName '--' datestr(now,'yymmdd--HHMM--SS--') strtok(data.currentConfig.configName,'.') '-'];
            mkdir([data.rootPath '_queue\'],dirName);
            data.currentConfig.save([data.rootPath '_queue\' dirName '\' data.currentConfig.configName '.cfg']);
            
            data.currentConfig.roiPosition(1:2) = posBackup;
            data.currentConfig.configName = configNameBackup;
        end
    end
end

guidata(gcf,data);

% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb34');
data = guidata(gcf);
% Load configuration
fileName = data.queue{get(handles.listbox3,'Value'),3};
path = data.queue{get(handles.listbox3,'Value'),1};
fullPath = [path fileName];
assignin('base', 'fullPath', fullPath);
evalin('base','cfg = Config.load(fullPath);');
if data.queue{get(handles.listbox3,'Value'),4}
    dataPath = [data.queue{get(handles.listbox3,'Value'),5}(1:end-5) 'p.dat'];
    assignin('base', 'dataPath', dataPath);
    evalin('base','data = Data.load(dataPath);');
    % Launch Imaris
    evalin('base','dis = Show(data);');
    % Display the data
    evalin('base','dis.points();');
    evalin('base','dis.models();');
    disp('Main: Show Models in Imaris: Done!')
else
    disp('Main: Show Models in Imaris: No data found!');
end


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb37'); % Overwrite
data = guidata(gcf);
% Check if configuration names match
configNameDest = data.queue{get(handles.listbox3,'Value'),3};
configNameSrc = [data.currentConfig.configName '.cfg'];
configFolderDest = data.queue{get(handles.listbox3,'Value'),1};
if ~strcmp([configFolderDest configNameSrc],data.configLoadedFrom)
    % Ask if config file should be replaced
    choice = questdlg('Do you want to overwrite the file anyway?','Configuration names do not match!','Yes','No, thank you!','Yes');
    switch choice
        case 'No, thank you!'
            % If no, return
            return;
        case 'Yes'
            configFolderDestOld = configFolderDest;
            [~,endPos] = regexp(configFolderDestOld,'--......--....--..--');
            configFolderDest = [configFolderDest(1:endPos) data.currentConfig.configName '\'];
            movefile(configFolderDestOld,configFolderDest,'f');
    end
end
% If the backup subfolder does not exist, create it
if ~exist([configFolderDest '\bak'],'file')
    mkdir([configFolderDest '\bak']);
end
% Create backup in subfolder of the old config file
try
    movefile([configFolderDest configNameDest],[configFolderDest 'bak\' configNameDest],'f');
catch
    fprintf('from = %s',[configFolderDest configNameDest])
    fprintf('to = %s',[configFolderDest 'bak\' configNameDest])
end
% Write the new config file
data.currentConfig.save([configFolderDest configNameSrc]);
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb38'); % Delete Data (d)
data = guidata(gcf);
dirPath = data.queue{get(handles.listbox3,'Value'),1};
list = dir(dirPath);
list = list(~vertcat(list.isdir));
for i=1:numel(list)
    if strcmp(list(i).name(end-5:end),'.d.dat')
        filePath = [dirPath list(i).name];
        delete(filePath);
    end
end
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb39'); % Delete Data (i)
data = guidata(gcf);
dirPath = data.queue{get(handles.listbox3,'Value'),1};
list = dir(dirPath);
list = list(~vertcat(list.isdir));
for i=1:numel(list)
    if strcmp(list(i).name(end-5:end),'.i.dat')
        filePath = [dirPath list(i).name];
        delete(filePath);
    end
end
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb40'); % Delete Data (p)
data = guidata(gcf);
dirPath = data.queue{get(handles.listbox3,'Value'),1};
list = dir(dirPath);
list = list(~vertcat(list.isdir));
for i=1:numel(list)
    if strcmp(list(i).name(end-5:end),'.p.dat')
        filePath = [dirPath list(i).name];
        delete(filePath);
    end
end
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb41'); % Duplicate data
data = guidata(gcf);
% Get source directory
dataFolderSrc = data.queue{get(handles.listbox3,'Value'),1};
% Update directory name
dataFolderDest = dataFolderSrc;
patternPos = strfind(dataFolderSrc,'--');
dataFolderDest(patternPos(1)+2:patternPos(1)+2+15) = datestr(now,'yymmdd--HHMM--SS');
% Copy directory
copyfile(dataFolderSrc,dataFolderDest);
% Detect the queued items
updateQueue(handles);


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb42');
data = guidata(gcf);
% Load timing object
fullPath = [data.queue{get(handles.listbox3,'Value'),5}(1:end-5) 'tim'];
assignin('base', 'fullPath', fullPath);
if exist(fullPath,'file')
    evalin('base','tim = Timing.load(fullPath);');
else
    disp('Main: Timer data not found!');
end


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pb43');
data = guidata(gcf);
dirPath = data.queue{get(handles.listbox3,'Value'),1};
system(['c:\Windows\explorer.exe ' '"' dirPath  '"']);
