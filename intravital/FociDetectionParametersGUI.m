function varargout = FociDetectionParametersGUI(varargin)
% FOCIDETECTIONPARAMETERSGUI MATLAB code for FociDetectionParametersGUI.fig
%      FOCIDETECTIONPARAMETERSGUI, by itself, creates a new FOCIDETECTIONPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = FOCIDETECTIONPARAMETERSGUI returns the handle to a new FOCIDETECTIONPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      FOCIDETECTIONPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOCIDETECTIONPARAMETERSGUI.M with the given input arguments.
%
%      FOCIDETECTIONPARAMETERSGUI('Property','Value',...) creates a new FOCIDETECTIONPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FociDetectionParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FociDetectionParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FociDetectionParametersGUI

% Last Modified by GUIDE v2.5 11-Jun-2014 15:08:13

% default parameters
defaultParameters.version = '1.0.0';
defaultParameters.fociDiameterRange = [1.86, 4.35];
defaultParameters.minSignalToBackgroundRatio = 1.0;
defaultParameters.minDoGResponse = 0.0;
defaultParameters.minHessianEigenRatio21 = 0.0;
defaultParameters.minHessianEigenRatio31 = 0.0;
defaultParameters.minDistanceToROIBoundary = 1.0;
defaultParameters.flagApplyFociDetectionModel = false;
defaultParameters.fociDetectionModelFile = [];

if nargin == 1 && ischar(varargin{1})
    
    switch( varargin{1} )
    
        case 'default'
    
            varargout{1} = defaultParameters;
            
        otherwise
            
            error( 'Invalid parameter code' );
    end
    
    return;
    
end

if nargin == 0
    varargin{1} = defaultParameters;
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FociDetectionParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FociDetectionParametersGUI_OutputFcn, ...
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


% --- Executes just before FociDetectionParametersGUI is made visible.
function FociDetectionParametersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FociDetectionParametersGUI (see VARARGIN)

% Choose default command line output for FociDetectionParametersGUI
handles.output = hObject;

% load parameter presets if it exists
[pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );

if exist( fullfile(pathstr, 'fociDetectionParameterPresets.mat'), 'file' )
    load( fullfile(pathstr, 'fociDetectionParameterPresets.mat') );
    handles.data.presetParameterDict = presetParameterDict;
else
    handles.data.presetParameterDict = containers.Map;
end

% initialize uicontrols of all parameters
initialParameters = varargin{1};
handles.data.curParameters = initialParameters;
handles.output = initialParameters;

% update UI controls
updateUIControls( handles );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FociDetectionParametersGUI wait for user response (see UIRESUME)
uiwait(handles.figFociDetectionParametersGUI);

function updateUIControls( handles )

    % foci diameter range
    set( handles.editMinFociDiameter, 'String', num2str(handles.data.curParameters.fociDiameterRange(1)) );
    set( handles.editMaxFociDiameter, 'String', num2str(handles.data.curParameters.fociDiameterRange(2)) );
    
    % minimum signal to background ratio
    set( handles.editMinSignalToBackgroundRatio, 'String', num2str(handles.data.curParameters.minSignalToBackgroundRatio) );

    % minimum DoG response
    set( handles.editMinDoGResponse, 'String', num2str(handles.data.curParameters.minDoGResponse) );

    % minimum hessian eigen value ratios
    set( handles.editMinHessianEigenRatio21, 'String', num2str(handles.data.curParameters.minHessianEigenRatio21) );
    set( handles.editMinHessianEigenRatio31, 'String', num2str(handles.data.curParameters.minHessianEigenRatio31) );
    
    % minimum distance to ROI boundary as a multiple of foci radius
    set( handles.editMinDistanceToROIBoundary, 'String', num2str(handles.data.curParameters.minDistanceToROIBoundary) );
    
    % region-merging
    set(handles.checkboxApplyFociDetectionModel, 'Value', handles.data.curParameters.flagApplyFociDetectionModel );
    if ~isempty(handles.data.curParameters.fociDetectionModelFile)
        set(handles.editFociDetectionModelFilePath, 'String', handles.data.curParameters.fociDetectionModelFile );
        set(handles.editFociDetectionModelFilePath, 'TooltipString', handles.data.curParameters.fociDetectionModelFile );
    end

% --- Outputs from this function are returned to the command line.
function varargout = FociDetectionParametersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);


function editMinFociDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to editMinFociDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinFociDiameter as text
%        str2double(get(hObject,'String')) returns contents of editMinFociDiameter as a double

    minFociDiameter = str2double( get(hObject,'String') );
    
    if isnan( minFociDiameter ) || minFociDiameter > handles.data.curParameters.fociDiameterRange(2)
        set( handles.editMinfociDiameter, 'String',  num2str(handles.data.curParameters.fociDiameterRange(1)) );
        errordlg( 'minimum foci diameter should be a real number and its value should be less than the maximum foci diameter' );        
        return;
    end

    handles.data.curParameters.fociDiameterRange(1) = minFociDiameter;
    
    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinFociDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinFociDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxFociDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxFociDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxFociDiameter as text
%        str2double(get(hObject,'String')) returns contents of editMaxFociDiameter as a double

    maxFociDiameter = str2double( get(hObject,'String') );
    
    if isnan( maxFociDiameter ) || maxFociDiameter < handles.data.curParameters.fociDiameterRange(1)
        set( handles.editMaxCellDiameter, 'String',  num2str(handles.data.curParameters.fociDiameterRange(2)) );
        errordlg( 'maximum foci diamter should be a real number and its value should be greater than the minimum foci diameter' );        
        return;
    end

    handles.data.curParameters.fociDiameterRange(2) = maxFociDiameter;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMaxFociDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxFociDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMinSignalToBackgroundRatio_Callback(hObject, eventdata, handles)
% hObject    handle to editMinSignalToBackgroundRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinSignalToBackgroundRatio as text
%        str2double(get(hObject,'String')) returns contents of editMinSignalToBackgroundRatio as a double

    minSignalToBackgroundRatio = str2double( get(hObject,'String') );
    
    if isnan( minSignalToBackgroundRatio ) || minSignalToBackgroundRatio < 1.0
        errordlg( 'Should be a real number > 1.' );        
        return;
    end

    handles.data.curParameters.minSignalToBackgroundRatio = minSignalToBackgroundRatio;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinSignalToBackgroundRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinSignalToBackgroundRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinDistanceToROIBoundary_Callback(hObject, eventdata, handles)
% hObject    handle to editMinDistanceToROIBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinDistanceToROIBoundary as text
%        str2double(get(hObject,'String')) returns contents of editMinDistanceToROIBoundary as a double

    minDistanceToROIBoundary = str2double( get(hObject,'String') );
    
    if isnan( minDistanceToROIBoundary ) || minDistanceToROIBoundary <= 0.0
        errordlg( 'Should be a real number > 0.' );        
        return;
    end

    handles.data.curParameters.minDistanceToROIBoundary = minDistanceToROIBoundary;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinDistanceToROIBoundary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinDistanceToROIBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinDoGResponse_Callback(hObject, eventdata, handles)
% hObject    handle to editMinDoGResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinDoGResponse as text
%        str2double(get(hObject,'String')) returns contents of editMinDoGResponse as a double

    minDoGResponse = str2double( get(hObject,'String') );
    
    if isnan( minDoGResponse ) || minDoGResponse < 0
        errordlg( 'Should be a real number >= 0.' );        
        return;
    end

    handles.data.curParameters.minDoGResponse = minDoGResponse;

    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMinDoGResponse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinDoGResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMinHessianEigenRatio21_Callback(hObject, eventdata, handles)
% hObject    handle to editMinHessianEigenRatio21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinHessianEigenRatio21 as text
%        str2double(get(hObject,'String')) returns contents of editMinHessianEigenRatio21 as a double

    minHessianEigenRatio21 = str2double( get(hObject,'String') );
    
    if isnan( minHessianEigenRatio21 ) || minHessianEigenRatio21 < 0
        errordlg( 'Should be a value between [0, 1].' );        
        return;
    end

    handles.data.curParameters.minHessianEigenRatio21 = minHessianEigenRatio21;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinHessianEigenRatio21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinHessianEigenRatio21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMinHessianEigenRatio31_Callback(hObject, eventdata, handles)
% hObject    handle to editMinHessianEigenRatio31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinHessianEigenRatio31 as text
%        str2double(get(hObject,'String')) returns contents of editMinHessianEigenRatio31 as a double

    minHessianEigenRatio31 = str2double( get(hObject,'String') );
    
    if isnan( minHessianEigenRatio31 ) || minHessianEigenRatio31 < 0
        errordlg( 'Should be a value between [0, 1].' );        
        return;
    end

    handles.data.curParameters.minHessianEigenRatio31 = minHessianEigenRatio31;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinHessianEigenRatio31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinHessianEigenRatio31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btnSetFociDetectionModelFile.
function btnSetFociDetectionModelFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetFociDetectionModelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isempty( handles.data.curParameters.fociDetectionModelFile )
        [initDir, ~, ~] = fileparts(handles.data.curParameters.fociDetectionModelFile);
    else
        initDir = pwd;
    end
    
    [fileName, pathName] = uigetfile(fullfile(initDir, '*.model'), 'Select Foci Detection Model File'); 
    
    if ~fileName 
        return;
    end
    
    fociDetectionModelFile = fullfile(pathName, fileName);
    
    handles.data.curParameters.fociDetectionModelFile = fociDetectionModelFile;
    set(handles.editFociDetectionModelFilePath, 'String', fociDetectionModelFile );
    set(handles.editFociDetectionModelFilePath, 'TooltipString', fociDetectionModelFile );
    
    % Update handles structure
    guidata(hObject, handles);

function editFociDetectionModelFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to editFociDetectionModelFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFociDetectionModelFilePath as text
%        str2double(get(hObject,'String')) returns contents of editFociDetectionModelFilePath as a double


% --- Executes during object creation, after setting all properties.
function editFociDetectionModelFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFociDetectionModelFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxApplyFociDetectionModel.
function checkboxApplyFociDetectionModel_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxApplyFociDetectionModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxApplyFociDetectionModel

    handles.data.curParameters.flagApplyFociDetectionModel = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in btnApplyParameters.
function btnApplyParameters_Callback(hObject, eventdata, handles)
% hObject    handle to btnApplyParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.data.curParameters.flagApplyFociDetectionModel && ...
       isempty(handles.data.curParameters.fociDetectionModelFile)     
       errordlg( 'Foci Detection was enabled but model file was not specified' ); 
       return;
    end

    handles.output = handles.data.curParameters;
    
    % Update handles structure
    guidata(hObject, handles);

    % close set parameters dialog
    close(handles.figFociDetectionParametersGUI);

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % close set parameters dialog
    close(handles.figFociDetectionParametersGUI);


% --- Executes when user attempts to close figFociDetectionParametersGUI.
function figFociDetectionParametersGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figFociDetectionParametersGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


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

    [fileName,pathName] = uigetfile( 'fociDetectionParameters*.mat', 'Select the parameter file' );   
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
function File_Parameter_Presets_Callback(hObject, eventdata, handles)
% hObject    handle to File_Parameter_Presets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
    
    fociDetectionParams = handles.data.curParameters;
    save( outFile, '-struct', 'fociDetectionParams' );


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
