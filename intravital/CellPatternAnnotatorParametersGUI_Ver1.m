function [varargout] = CellPatternAnnotatorParametersGUI_Ver1( varargin )
% CELLPATTERNANNOTATORPARAMETERSGUI_VER1 MATLAB code for CellPatternAnnotatorParametersGUI_Ver1.fig
%      CELLPATTERNANNOTATORPARAMETERSGUI_VER1, by itself, creates a new CELLPATTERNANNOTATORPARAMETERSGUI_VER1 or raises the existing
%      singleton*.
%
%      H = CELLPATTERNANNOTATORPARAMETERSGUI_VER1 returns the handle to a new CELLPATTERNANNOTATORPARAMETERSGUI_VER1 or the handle to
%      the existing singleton*.
%
%      CELLPATTERNANNOTATORPARAMETERSGUI_VER1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLPATTERNANNOTATORPARAMETERSGUI_VER1.M with the given input arguments.
%
%      CELLPATTERNANNOTATORPARAMETERSGUI_VER1('Property','Value',...) creates a new CELLPATTERNANNOTATORPARAMETERSGUI_VER1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellPatternAnnotatorParametersGUI_Ver1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellPatternAnnotatorParametersGUI_Ver1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellPatternAnnotatorParametersGUI_Ver1

% Last Modified by GUIDE v2.5 27-Mar-2013 17:20:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellPatternAnnotatorParametersGUI_Ver1_OpeningFcn, ...
                   'gui_OutputFcn',  @CellPatternAnnotatorParametersGUI_Ver1_OutputFcn, ...
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


% --- Executes just before CellPatternAnnotatorParametersGUI_Ver1 is made visible.
function CellPatternAnnotatorParametersGUI_Ver1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellPatternAnnotatorParametersGUI_Ver1 (see VARARGIN)

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
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellPatternAnnotatorParametersGUI_Ver1 wait for user response (see UIRESUME)
uiwait(handles.figCellPatternAnnotatorParametersGUI);


% --- Outputs from this function are returned to the command line.
function varargout = CellPatternAnnotatorParametersGUI_Ver1_OutputFcn(hObject, eventdata, handles) 
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

    minCellDiameter = handles.data.curParameters.cellDiamterRange(1);
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
