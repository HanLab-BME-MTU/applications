function varargout = CellCycleStateIdentificationDataSelectionGUI(varargin)
% CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI MATLAB code for CellCycleStateIdentificationDataSelectionGUI.fig
%      CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI, by itself, creates a new CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI or raises the existing
%      singleton*.
%
%      H = CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI returns the handle to a new CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI or the handle to
%      the existing singleton*.
%
%      CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI.M with the given input arguments.
%
%      CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI('Property','Value',...) creates a new CELLCYCLESTATEIDENTIFICATIONDATASELECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellCycleStateIdentificationDataSelectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellCycleStateIdentificationDataSelectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellCycleStateIdentificationDataSelectionGUI

% Last Modified by GUIDE v2.5 11-Jul-2014 15:18:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellCycleStateIdentificationDataSelectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellCycleStateIdentificationDataSelectionGUI_OutputFcn, ...
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


% --- Executes just before CellCycleStateIdentificationDataSelectionGUI is made visible.
function CellCycleStateIdentificationDataSelectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellCycleStateIdentificationDataSelectionGUI (see VARARGIN)

% Choose default command line output for CellCycleStateIdentificationDataSelectionGUI
handles.output = [];
handles.lastVisitedDir = [];

handles.data.flagAlignFucciDataToNuclearMarker = get(handles.checkboxAlignFucciToNuclearMarker, 'Value');
handles.data.nuclearMarker = [];
handles.data.fucciGemininMarker = [];
handles.data.fucciCdt1Marker = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellCycleStateIdentificationDataSelectionGUI wait for user response (see UIRESUME)
uiwait(handles.figCellCycleStateIdentificationDataSelectionGUI);


% --- Outputs from this function are returned to the command line.
function varargout = CellCycleStateIdentificationDataSelectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % close set parameters dialog
    close(handles.figCellCycleStateIdentificationDataSelectionGUI);

% --- Executes on button press in btnOk.
function btnOk_Callback(hObject, eventdata, handles)
% hObject    handle to btnOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % make sure files for all three channels were set
    if isempty(handles.data.nuclearMarker) || ...
       isempty(handles.data.fucciGemininMarker) || ...
       isempty(handles.data.fucciCdt1Marker)
        errordlg('Files for all three markers must be specified.' )
        return;
    end

    % make sure that the metadata is consistent accross the channels
    checkFieldList = {'imageSize', 'pixelSize', 'bitsPerPixel', 'pixelType', 'timepointId'};
    
    for i = 1:numel(checkFieldList)

        fname = checkFieldList{i};
        if any(handles.data.nuclearMarker.metadata.(fname) ~= handles.data.fucciGemininMarker.metadata.(fname)) || ...
           any(handles.data.nuclearMarker.metadata.(fname) ~= handles.data.fucciCdt1Marker.metadata.(fname))
            errordlg( sprintf('%s of all three markers must be the same.', fname) );
            return;
        end
        
    end

    % prepare output
    handles.output = [];
    
        % should channels be registered
        handles.output.flagAlignFucciDataToNuclearMarker = handles.data.flagAlignFucciDataToNuclearMarker;    

        % metadata
        metadata = handles.data.nuclearMarker.metadata;

        metadata.dataFilePath = {handles.data.nuclearMarker.metadata.dataFilePath, ...
                                 handles.data.fucciGemininMarker.metadata.dataFilePath, ...
                                 handles.data.fucciCdt1Marker.metadata.dataFilePath};

        metadata.format = {handles.data.nuclearMarker.metadata.format, ...
                           handles.data.fucciGemininMarker.metadata.format, ...
                           handles.data.fucciCdt1Marker.metadata.format};

        metadata.numSeries = [handles.data.nuclearMarker.metadata.numSeries, ...
                              handles.data.fucciGemininMarker.metadata.numSeries, ...
                              handles.data.fucciCdt1Marker.metadata.numSeries];

        metadata.seriedId = [handles.data.nuclearMarker.metadata.seriesId, ...
                             handles.data.fucciGemininMarker.metadata.seriesId, ...
                             handles.data.fucciCdt1Marker.metadata.seriesId];

        metadata.channelId = [handles.data.nuclearMarker.metadata.channelId, ...
                             handles.data.fucciGemininMarker.metadata.channelId, ...
                             handles.data.fucciCdt1Marker.metadata.channelId];
                         
        metadata.numChannels = 3;
        metadata.channelNames = {'Nuclei', 'FUCCI Geminin ', 'FUCCI Cdt1'};
        
        handles.output.metadata = metadata;
        
        % metadata
        handles.output.imageData = [handles.data.nuclearMarker.imageData, ...
                                    handles.data.fucciGemininMarker.imageData, ...
                                    handles.data.fucciCdt1Marker.imageData];
    
    % Update handles structure
    guidata(hObject, handles);

    % close set parameters dialog
    close(handles.figCellCycleStateIdentificationDataSelectionGUI);
    

% --- Executes on button press in btnSelectFucciCdt1File.
function btnSelectFucciCdt1File_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectFucciCdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the user to select the oif file for the nucleus channel
    status = bfCheckJavaPath(1);
    [fileName,pathName] = uigetfile(bfGetFileExtensions, ...
                                    'Select the file contaning the FUCCI Cdt1 marker for G1 Phase', ...
                                    handles.lastVisitedDir);
                            
    if ~fileName 
        return;
    end
                                
    handles.lastVisitedDir = pathName;
    dataFilePath = fullfile(pathName, fileName);
    
    % load image data
    [imageData, metadata] = uiGetBioImageData(dataFilePath, 'channelDescriptionsAndNamesList', {{'FUCCI Cdt1 marker for G1', 'channelId'}});
    
    if isempty(imageData)
        return;
    end
    
    handles.data.fucciCdt1Marker.metadata = metadata;
    handles.data.fucciCdt1Marker.imageData = imageData(metadata.channelId);

    % update ui controls
    set(handles.editFUCCICdt1File, 'String', dataFilePath);
    set(handles.editFUCCICdt1File, 'TooltipString', dataFilePath);
    
    % Update handles structure
    guidata(hObject, handles);

function editFUCCICdt1File_Callback(hObject, eventdata, handles)
% hObject    handle to editFUCCICdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFUCCICdt1File as text
%        str2double(get(hObject,'String')) returns contents of editFUCCICdt1File as a double


% --- Executes during object creation, after setting all properties.
function editFUCCICdt1File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFUCCICdt1File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSelectFucciGemininFile.
function btnSelectFucciGemininFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the user to select the oif file for the nucleus channel
    status = bfCheckJavaPath(1);
    [fileName,pathName] = uigetfile(bfGetFileExtensions, ...
                                    'Select the file contaning FUCCI Geminin marker for S/G2/M Phases', ...
                                    handles.lastVisitedDir);
                                
    if ~fileName 
        return;
    end
                                
    handles.lastVisitedDir = pathName;
    dataFilePath = fullfile(pathName, fileName);
    
    % load image data
    [imageData, metadata] = uiGetBioImageData(dataFilePath, 'channelDescriptionsAndNamesList', {{'FUCCI Geminin marker for S/G2/M', 'channelId'}});
    
    if isempty(imageData)
        return;
    end
    
    handles.data.fucciGemininMarker.metadata = metadata;
    handles.data.fucciGemininMarker.imageData = imageData(metadata.channelId);

    % update ui controls
    set(handles.editFucciGemininFile, 'String', dataFilePath);
    set(handles.editFucciGemininFile, 'TooltipString', dataFilePath);
    
    % Update handles structure
    guidata(hObject, handles);

function editFucciGemininFile_Callback(hObject, eventdata, handles)
% hObject    handle to editFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFucciGemininFile as text
%        str2double(get(hObject,'String')) returns contents of editFucciGemininFile as a double


% --- Executes during object creation, after setting all properties.
function editFucciGemininFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFucciGemininFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSelectNuclearMarkerFile.
function btnSelectNuclearMarkerFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the user to select the oif file for the nucleus channel
    status = bfCheckJavaPath(1);
    [fileName,pathName] = uigetfile(bfGetFileExtensions, ...
                                    'Select the file contaning the nuclear marker - Histone-2B', ...
                                    handles.lastVisitedDir);
                   
    if ~fileName 
        return;
    end
                                
    dataFilePath = fullfile(pathName, fileName);
    handles.lastVisitedDir = pathName;
    
    % load image data
    [imageData, metadata] = uiGetBioImageData(dataFilePath, 'channelDescriptionsAndNamesList', {{'Nuclear Marker Channel', 'channelId'}});
    
    if isempty(imageData)
        return;
    end
    
    handles.data.nuclearMarker.metadata = metadata;
    handles.data.nuclearMarker.imageData = imageData(metadata.channelId);

    % update ui controls
    set(handles.editNuclearMarkerFile, 'String', dataFilePath);
    set(handles.editNuclearMarkerFile, 'TooltipString', dataFilePath);
    
    % Update handles structure
    guidata(hObject, handles);

function editNuclearMarkerFile_Callback(hObject, eventdata, handles)
% hObject    handle to editNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNuclearMarkerFile as text
%        str2double(get(hObject,'String')) returns contents of editNuclearMarkerFile as a double



% --- Executes during object creation, after setting all properties.
function editNuclearMarkerFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNuclearMarkerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxAlignFucciToNuclearMarker.
function checkboxAlignFucciToNuclearMarker_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAlignFucciToNuclearMarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAlignFucciToNuclearMarker

    handles.data.flagAlignFucciDataToNuclearMarker = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);


% --- Executes when user attempts to close figCellCycleStateIdentificationDataSelectionGUI.
function figCellCycleStateIdentificationDataSelectionGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figCellCycleStateIdentificationDataSelectionGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
