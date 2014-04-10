function [varargout] = NucleiSegmentationParametersGUI( varargin )
% NUCLEISEGMENTATIONPARAMETERSGUI MATLAB code for NucleiSegmentationParametersGUI.fig
%      NUCLEISEGMENTATIONPARAMETERSGUI, by itself, creates a new NUCLEISEGMENTATIONPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = NUCLEISEGMENTATIONPARAMETERSGUI returns the handle to a new NUCLEISEGMENTATIONPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      NUCLEISEGMENTATIONPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUCLEISEGMENTATIONPARAMETERSGUI.M with the given input arguments.
%
%      NUCLEISEGMENTATIONPARAMETERSGUI('Property','Value',...) creates a new NUCLEISEGMENTATIONPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NucleiSegmentationParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NucleiSegmentationParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NucleiSegmentationParametersGUI

% Last Modified by GUIDE v2.5 01-Apr-2014 17:24:15

% default parameters
defaultParameters.cellDiameterRange = [12, 20];
defaultParameters.minCellVolume = 400;
defaultParameters.seedPointDetectionAlgorithm = 'AdaptiveMultiscaleLoG';
defaultParameters.flagIgnoreXYBorderCells = true;
defaultParameters.flagPerformRegionMerging = false;
defaultParameters.regionMergingModelFile = [];

if nargin == 1 && ischar(varargin{1})
    
    switch( varargin{1} )
    
        case 'defaultWithRegionMerging'
            
            defaultParameters.cellDiameterRange = [8, 20];
            
        case 'defaultWithoutRegionMerging'
            
            defaultParameters.cellDiameterRange = [12, 20];

        case 'defaultRegionMergingTraining'
            
            defaultParameters.cellDiameterRange = [8, 20];
            defaultParameters.minCellVolume = 201.06;
            
        otherwise
            
            error( 'Invalid parameter code' );
    end
    
    varargout{1} = defaultParameters;
    return;
end

if nargin == 0
    varargin{1} = defaultParameters;
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NucleiSegmentationParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NucleiSegmentationParametersGUI_OutputFcn, ...
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


% --- Executes just before NucleiSegmentationParametersGUI is made visible.
function NucleiSegmentationParametersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NucleiSegmentationParametersGUI (see VARARGIN)

% load parameter presets if it exists
[pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );

if exist( fullfile(pathstr, 'nucleiSegmentationParameterPresets.mat'), 'file' )
    load( fullfile(pathstr, 'nucleiSegmentationParameterPresets.mat') );
    handles.data.presetParameterDict = presetParameterDict;
else
    handles.data.presetParameterDict = containers.Map;
end

% initialize uicontrols of all parameters
initialParameters = varargin{1};
handles.data.curParameters = initialParameters;
handles.output = initialParameters;

updateUIControls( handles );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NucleiSegmentationParametersGUI wait for user response (see UIRESUME)
uiwait(handles.figCellPatternAnnotatorParametersGUI);

function updateUIControls( handles )

    % cell diameter range
    set( handles.editMinCellDiameter, 'String', num2str(handles.data.curParameters.cellDiameterRange(1)) );
    set( handles.editMaxCellDiameter, 'String', num2str(handles.data.curParameters.cellDiameterRange(2)) );
    
    % minimum cell volume
    set( handles.editMinCellVolume, 'String', num2str(handles.data.curParameters.minCellVolume) );
    
    minCellVolumeSuggested = 0.75 * (4/3) * pi * (0.5 * handles.data.curParameters.cellDiameterRange(1))^3;    
    set( handles.labelSuggestedMinVolume, 'String', sprintf( 'suggested value - %f', minCellVolumeSuggested ) );
    
    % cell seed point detection method
    handles.data.seedPointDetectionAlgorithmList = { 'IntensityMaxima', ...
                                                     'MultiscaleLoG', ...
                                                     'AdaptiveMultiscaleLoG', ...
                                                     'MultiscaleLoBG' };

    initMethodId = find( strcmpi( handles.data.seedPointDetectionAlgorithmList, ...
                                  handles.data.curParameters.seedPointDetectionAlgorithm ) );
                            
    if isempty( initMethodId )
        error( 'ERROR: unknown seed point detection algorithm' );
    end
    
    set( handles.popupSeedPointDetectionMethod, 'String', handles.data.seedPointDetectionAlgorithmList );    
    set( handles.popupSeedPointDetectionMethod, 'Value', initMethodId );    
    
    % post-processing
    set(handles.chkboxIgnoreXYBorderCells, 'Value', handles.data.curParameters.flagIgnoreXYBorderCells );

    % region-merging
    set(handles.chkboxPerformRegionMerging, 'Value', handles.data.curParameters.flagPerformRegionMerging );
    if handles.data.curParameters.flagPerformRegionMerging && ~isempty(handles.data.curParameters.regionMergingModelFile)
        set(handles.editRegionMergingModelFilePath, 'String', handles.data.curParameters.regionMergingModelFile );
        set(handles.editRegionMergingModelFilePath, 'TooltipString', handles.data.curParameters.regionMergingModelFile );
    end


% --- Outputs from this function are returned to the command line.
function varargout = NucleiSegmentationParametersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.data.presetParameterDict.isempty() 
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    presetParameterDict = handles.data.presetParameterDict;
    save(fullfile(pathstr, 'nucleiSegmentationParameterPresets.mat'), 'presetParameterDict' );
end

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


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Load_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fileName,pathName] = uigetfile( 'nucleiSegmentationParameters*.mat', 'Select the parameter file' );   
    if ~fileName 
        return;
    end
    
    params = load( fullfile(pathName, fileName) );
    if ~isequal( fieldnames(params), fieldnames(handles.data.curParameters) )
        errordlg( 'Invalid Parameter File' );
        return;
    end
    
    handles.data.curParameters = params;
    
    updateUIControls( handles );

   % Update handles structure
   guidata(hObject, handles);
    
% --------------------------------------------------------------------
function File_Save_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_Save_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fileName,pathName] = uiputfile( '*.mat', 'Select the parameter file' );   

    if ~ischar(fileName)
        return;
    end
    
    outFile = fullfile(pathName, fileName);
    
    if exist( outFile, 'file' )
        answer = questdlg( 'A file with that name already exists. Do you want to overwrite it?', ...
                           '', 'Yes', 'No', 'No' );
                       
        if strcmp( answer, 'No' )
            return;
        end
    end
    
    nucleiSegmentationParams = handles.data.curParameters;
    save( outFile, '-struct', 'nucleiSegmentationParams' );
        
% --------------------------------------------------------------------
function File_Parameter_Presets_Callback(hObject, eventdata, handles)
% hObject    handle to File_Parameter_Presets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Load_Preset_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Preset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   if handles.data.presetParameterDict.isempty()
        warndlg( 'No presets were found' );
        return;
   end
   
   presetNames = handles.data.presetParameterDict.keys();
   [selid, flagOk] = listdlg( 'ListString', presetNames, ...
                              'SelectionMode', 'single', ...
                              'ListSize', [400 300], ...
                              'Name', 'Select the preset you want to load');   
                
   if isempty(selid)
       return;
   end

   handles.data.curParameters = handles.data.presetParameterDict( presetNames{selid} );

   updateUIControls( handles );

   % Update handles structure
   guidata(hObject, handles);

% --------------------------------------------------------------------
function File_Save_Preset_Callback(hObject, eventdata, handles)
% hObject    handle to File_Save_Preset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   presetName = inputdlg( 'Enter preset name', '' ); 
   if ~isempty(presetName)
       handles.data.presetParameterDict(presetName{1}) = handles.data.curParameters;
   end
   
   % Update handles structure
   guidata(hObject, handles);

% --------------------------------------------------------------------
function File_Delete_Preset_Callback(hObject, eventdata, handles)
% hObject    handle to File_Delete_Preset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   if handles.data.presetParameterDict.isempty()
        warndlg( 'No presets were found' );
        return;
   end
   
   presetNames = handles.data.presetParameterDict.keys();
   [selid, flagOk] = listdlg( 'ListString', presetNames, ...
                              'SelectionMode', 'single', ...
                              'Name', 'Select the preset you want to remove');   
                
   if isempty(selid)
       return;
   end

   handles.data.presetParameterDict.remove( presetNames{selid} );

   % Update handles structure
   guidata(hObject, handles);
