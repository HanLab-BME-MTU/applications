function [varargout] = CellSegmentationQualityAnnotatorParametersGUI( varargin )
% CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI MATLAB code for CellSegmentationQualityAnnotatorParametersGUI.fig
%      CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI, by itself, creates a new CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI returns the handle to a new CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI.M with the given input arguments.
%
%      CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI('Property','Value',...) creates a new CELLSEGMENTATIONQUALITYANNOTATORPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellSegmentationQualityAnnotatorParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellSegmentationQualityAnnotatorParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellSegmentationQualityAnnotatorParametersGUI

% Last Modified by GUIDE v2.5 06-Mar-2014 12:43:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellSegmentationQualityAnnotatorParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellSegmentationQualityAnnotatorParametersGUI_OutputFcn, ...
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


% --- Executes just before CellSegmentationQualityAnnotatorParametersGUI is made visible.
function CellSegmentationQualityAnnotatorParametersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellSegmentationQualityAnnotatorParametersGUI (see VARARGIN)

% get initial parameters
initialParameters = varargin{1};

% set default output
handles.output = initialParameters;

% initialize uicontrols of all parameters
handles.data.initialParameters = initialParameters;
handles.data.curParameters = initialParameters;

    % cell diameter range
    set( handles.editMinCellDiameter, 'String', num2str(initialParameters.cellDiameterRange(1)) );
    set( handles.editMaxCellDiameter, 'String', num2str(initialParameters.cellDiameterRange(2)) );
    
    % minimum cell volume
    set( handles.editMinCellVolume, 'String', num2str(initialParameters.minCellVolume) );
    
    minCellVolumeSuggested = 0.75 * (4/3) * pi * (0.5 * initialParameters.cellDiameterRange(1))^3;    
    set( handles.labelSuggestedMinVolume, 'String', sprintf( 'suggested value - %f', minCellVolumeSuggested ) );
    
    % cell seed point detection method
    handles.data.seedPointDetectionAlgorithmList = { 'IntensityMaxima', ...
                                                     'MultiscaleLoG', ...
                                                     'AdaptiveMultiscaleLoG', ...
                                                     'MultiscaleLoBG' };

    initMethodId = find( strcmpi( handles.data.seedPointDetectionAlgorithmList, ...
                                  initialParameters.seedPointDetectionAlgorithm ) );
                            
    if isempty( initMethodId )
        error( 'ERROR: unknown seed point detection algorithm' );
    end
    
    set( handles.popupSeedPointDetectionMethod, 'String', handles.data.seedPointDetectionAlgorithmList );    
    set( handles.popupSeedPointDetectionMethod, 'Value', initMethodId );    
    
    % post-processing
    set(handles.chkboxIgnoreXYBorderCells, 'Value', initialParameters.flagIgnoreXYBorderCells );

    % region-merging
    set(handles.chkboxPerformRegionMerging, 'Value', initialParameters.flagPerformRegionMerging );
    if initialParameters.flagPerformRegionMerging && ~isempty(initialParameters.regionMergingModelFile)
        set(handles.editRegionMergingModelFilePath, 'String', initialParameters.regionMergingModelFile );
        set(handles.editRegionMergingModelFilePath, 'TooltipString', initialParameters.regionMergingModelFile );
    end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellSegmentationQualityAnnotatorParametersGUI wait for user response (see UIRESUME)
uiwait(handles.figCellPatternAnnotatorParametersGUI);


% --- Outputs from this function are returned to the command line.
function varargout = CellSegmentationQualityAnnotatorParametersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);

function editMinCellDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to editMinCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinCellDiameter as text
%        str2double(get(hObject,'String')) returns contents of editMinCellDiameter as a double

    minCellDiameter = str2double( get(hObject,'String') );
    
    if isnan( minCellDiameter ) || minCellDiameter > handles.data.curParameters.cellDiameterRange(2)
        set( handles.editMinCellDiameter, 'String',  num2str(handles.data.curParameters.cellDiameterRange(1)) );
        errordlg( 'minimum cell diameter should be a real number and its value should be less than the maximum cell diameter' );        
        return;
    end

    handles.data.curParameters.cellDiameterRange(1) = minCellDiameter;
    
    minCellVolumeSuggested = 0.75 * (4/3) * pi * (0.5 * handles.data.curParameters.cellDiameterRange(1))^3;    
    set( handles.labelSuggestedMinVolume, 'String', sprintf( 'suggested value - %f', minCellVolumeSuggested ) );
    
    % Update handles structure
    guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function editMinCellDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxCellDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxCellDiameter as text
%        str2double(get(hObject,'String')) returns contents of editMaxCellDiameter as a double

    maxCellDiameter = str2double( get(hObject,'String') );
    
    if isnan( maxCellDiameter ) || maxCellDiameter < handles.data.curParameters.cellDiameterRange(1)
        set( handles.editMaxCellDiameter, 'String',  num2str(handles.data.curParameters.cellDiameterRange(2)) );
        errordlg( 'maximum cell diameter should be a real number and its value should be greater than the minimum cell diameter' );        
        return;
    end

    handles.data.curParameters.cellDiameterRange(2) = maxCellDiameter;
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMaxCellDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupSeedPointDetectionMethod.
function popupSeedPointDetectionMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popupSeedPointDetectionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupSeedPointDetectionMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupSeedPointDetectionMethod

    selectedAlgorithmId = get(hObject,'Value');
    handles.data.curParameters.seedPointDetectionAlgorithm = handles.data.seedPointDetectionAlgorithmList{ selectedAlgorithmId };
    
    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupSeedPointDetectionMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupSeedPointDetectionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMinCellVolume_Callback(hObject, eventdata, handles)
% hObject    handle to editMinCellVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinCellVolume as text
%        str2double(get(hObject,'String')) returns contents of editMinCellVolume as a double

    minCellVolume = str2double( get(hObject,'String') );
    
    if isnan( minCellVolume ) 
        errordlg( 'minimum cell volume should be a real number.' );        
        return;
    end

    minCellDiameter = handles.data.curParameters.cellDiameterRange(1);
    minDiameterSphereVolume = (4.0 * pi * (0.5 * minCellDiameter)^3 / 3.0);
    
    if minCellVolume > minDiameterSphereVolume
        warndlg( 'specified minimum cell volume is greater than the volume of sphere with specified minimum cell diameter' );        
    end
    
    handles.data.curParameters.minCellVolume = minCellVolume;
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMinCellVolume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinCellVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figCellPatternAnnotatorParametersGUI.
function figCellPatternAnnotatorParametersGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figCellPatternAnnotatorParametersGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

% --- Executes on button press in btnApplyParameters.
function btnApplyParameters_Callback(hObject, eventdata, handles)
% hObject    handle to btnApplyParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.data.curParameters.flagPerformRegionMerging && ...
       isempty(handles.data.curParameters.regionMergingModelFile)     
       errordlg( 'Region merging was enabled but model file was not specified' ); 
       return;
    end

    handles.output = handles.data.curParameters;
    
    % Update handles structure
    guidata(hObject, handles);

    % close set parameters dialog
    close(handles.figCellPatternAnnotatorParametersGUI);
    
% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % close set parameters dialog
    close(handles.figCellPatternAnnotatorParametersGUI);

function editRegionMergingModelFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to editRegionMergingModelFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRegionMergingModelFilePath as text
%        str2double(get(hObject,'String')) returns contents of editRegionMergingModelFilePath as a double


% --- Executes during object creation, after setting all properties.
function editRegionMergingModelFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRegionMergingModelFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSetRegionMergingModelFile.
function btnSetRegionMergingModelFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetRegionMergingModelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isempty( handles.data.curParameters.regionMergingModelFile )
        [initDir, ~, ~] = fileparts(handles.data.curParameters.regionMergingModelFile);
    else
        initDir = pwd;
    end
    
    [fileName, pathName] = uigetfile(fullfile(initDir, '*.model'), 'Select Region Merging Model File'); 
    
    if ~fileName 
        return;
    end
    
    regionMergingModelFile = fullfile(pathName, fileName);
    
    handles.data.curParameters.regionMergingModelFile = regionMergingModelFile;
    set(handles.editRegionMergingModelFilePath, 'String', regionMergingModelFile );
    set(handles.editRegionMergingModelFilePath, 'TooltipString', regionMergingModelFile );
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in chkboxIgnoreXYBorderCells.
function chkboxIgnoreXYBorderCells_Callback(hObject, eventdata, handles)
% hObject    handle to chkboxIgnoreXYBorderCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkboxIgnoreXYBorderCells

    handles.data.curParameters.flagIgnoreXYBorderCells = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in chkboxPerformRegionMerging.
function chkboxPerformRegionMerging_Callback(hObject, eventdata, handles)
% hObject    handle to chkboxPerformRegionMerging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkboxPerformRegionMerging

    handles.data.curParameters.flagPerformRegionMerging = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
