function [varargout] = CellCycleStateIdentificationParametersGUI( varargin )
% CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI MATLAB code for CellCycleStateIdentificationParametersGUI.fig
%      CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI, by itself, creates a new CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI returns the handle to a new CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI.M with the given input arguments.
%
%      CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI('Property','Value',...) creates a new CELLCYCLESTATEIDENTIFICATIONPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellCycleStateIdentificationParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellCycleStateIdentificationParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellCycleStateIdentificationParametersGUI

% Last Modified by GUIDE v2.5 03-Apr-2014 13:18:23

% default parameters
defaultParameters.minCellROIOverlap = 0.5;
defaultParameters.flagPerformCellCycleStateIdentification = false;
defaultParameters.cellCycleStateIdentificationModelDir = [];
defaultParameters.cellCycleStateIdentificationModelType = 'G1_S_G2_M';

if nargin == 1 && ischar(varargin{1}) && strcmpi( varargin{1}, 'default' )
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
                   'gui_OpeningFcn', @CellCycleStateIdentificationParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellCycleStateIdentificationParametersGUI_OutputFcn, ...
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


% --- Executes just before CellCycleStateIdentificationParametersGUI is made visible.
function CellCycleStateIdentificationParametersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellCycleStateIdentificationParametersGUI (see VARARGIN)

% get initial parameters
initialParameters = varargin{1};

% set default output
handles.output = initialParameters;

% initialize uicontrols of all parameters
handles.data.initialParameters = initialParameters;
handles.data.curParameters = initialParameters;

    set(handles.editMinCellROIOverlap, 'String', num2str(initialParameters.minCellROIOverlap) );

    set(handles.chkboxPerformCellCycleStateClassification, 'Value', initialParameters.flagPerformCellCycleStateIdentification );

    if initialParameters.flagPerformCellCycleStateIdentification && ...
       ~isempty(initialParameters.cellCycleStateIdentificationModelDir)
        set(handles.editCellCycleStateIdenficationModelDir, 'String', initialParameters.cellCycleStateIdentificationModelDir );
        set(handles.editCellCycleStateIdenficationModelDir, 'TooltipString', initialParameters.cellCycleStateIdentificationModelDir );
    end
    
    handles.cellCycleStateIdentificationModelTypeList = get( handles.popupCellCycleStateIdentificationModelType, 'String');    

    initModelTypeId = find( strcmpi( handles.cellCycleStateIdentificationModelTypeList, ...
                                     initialParameters.cellCycleStateIdentificationModelType ) );
    
    if isempty( initModelTypeId )
        error( 'ERROR: invalid model type' );
    end
                                 
    set(handles.popupCellCycleStateIdentificationModelType, 'Value', initModelTypeId);
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellCycleStateIdentificationParametersGUI wait for user response (see UIRESUME)
uiwait(handles.figCellCycleStateIdentificationParametersGUI);

% --- Outputs from this function are returned to the command line.
function varargout = CellCycleStateIdentificationParametersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);

% --- Executes when user attempts to close figCellCycleStateIdentificationParametersGUI.
function figCellCycleStateIdentificationParametersGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figCellCycleStateIdentificationParametersGUI (see GCBO)
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

    if handles.data.curParameters.flagPerformCellCycleStateIdentification && ...
       isempty(handles.data.curParameters.cellCycleStateIdentificationModelDir)     
       errordlg( 'Cell cycle state classification was enabled but model diretory was not specified' ); 
       return;
    end
    
    handles.output = handles.data.curParameters;
    
    % Update handles structure
    guidata(hObject, handles);

    % close set parameters dialog
    close(handles.figCellCycleStateIdentificationParametersGUI);
    
% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % close set parameters dialog
    close(handles.figCellCycleStateIdentificationParametersGUI);

% --- Executes on button press in chkboxPerformCellCycleStateClassification.
function chkboxPerformCellCycleStateClassification_Callback(hObject, eventdata, handles)
% hObject    handle to chkboxPerformCellCycleStateClassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chkboxPerformCellCycleStateClassification

    handles.data.curParameters.flagPerformCellCycleStateIdentification = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in btnSetCellCycleStateIdentificationModelDir.
function btnSetCellCycleStateIdentificationModelDir_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetCellCycleStateIdentificationModelDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isempty( handles.data.curParameters.cellCycleStateIdentificationModelDir  )
        initDir = handles.data.curParameters.cellCycleStateIdentificationModelDir;
    else
        initDir = pwd;
    end
    
    modelDir = uigetdir(initDir, 'Select Directory Containing the Cell Cycle State Identification Model'); 
    
    if ~modelDir
        return;
    end

    if ~exist( fullfile(modelDir, 'G1_S_G2_M.model'), 'file')
        errordlg( 'cell cycle state classification model dir must contain a file named G1_S_G2_M.model' );
        return;
    end
    
    handles.data.curParameters.cellCycleStateIdentificationModelDir = modelDir;
    set(handles.editCellCycleStateIdenficationModelDir, 'String', modelDir );
    set(handles.editCellCycleStateIdenficationModelDir, 'TooltipString', modelDir);
    
    % Update handles structure
    guidata(hObject, handles);

function editCellCycleStateIdenficationModelDir_Callback(hObject, eventdata, handles)
% hObject    handle to editCellCycleStateIdenficationModelDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCellCycleStateIdenficationModelDir as text
%        str2double(get(hObject,'String')) returns contents of editCellCycleStateIdenficationModelDir as a double


% --- Executes during object creation, after setting all properties.
function editCellCycleStateIdenficationModelDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCellCycleStateIdenficationModelDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMinCellROIOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to editMinCellROIOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinCellROIOverlap as text
%        str2double(get(hObject,'String')) returns contents of editMinCellROIOverlap as a double

    minCellROIOverlap = str2double(get(hObject,'String'));
    
    if isnan(minCellROIOverlap) || minCellROIOverlap < 0.0 || minCellROIOverlap >= 1.0
        errordlg( 'minimum cell ROI overlap must be a value between 0 and 1.' );        
        return;
    end
    
    handles.data.curParameters.minCellROIOverlap = minCellROIOverlap;

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function editMinCellROIOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinCellROIOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupCellCycleStateIdentificationModelType.
function popupCellCycleStateIdentificationModelType_Callback(hObject, eventdata, handles)
% hObject    handle to popupCellCycleStateIdentificationModelType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupCellCycleStateIdentificationModelType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupCellCycleStateIdentificationModelType

    selectedModelTypeId = get(hObject,'Value');
    handles.data.curParameters.cellCycleStateIdentificationModelType = handles.cellCycleStateIdentificationModelTypeList{ selectedModelTypeId };
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupCellCycleStateIdentificationModelType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupCellCycleStateIdentificationModelType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
