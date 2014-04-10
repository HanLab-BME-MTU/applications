function varargout = InvivoCytometer(varargin)
% INVIVOCYTOMETER MATLAB code for InvivoCytometer.fig
%      INVIVOCYTOMETER, by itself, creates a new INVIVOCYTOMETER or raises the existing
%      singleton*.
%
%      H = INVIVOCYTOMETER returns the handle to a new INVIVOCYTOMETER or the handle to
%      the existing singleton*.
%
%      INVIVOCYTOMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INVIVOCYTOMETER.M with the given input arguments.
%
%      INVIVOCYTOMETER('Property','Value',...) creates a new INVIVOCYTOMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InvivoCytometer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InvivoCytometer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InvivoCytometer

% Last Modified by GUIDE v2.5 07-Apr-2014 12:40:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InvivoCytometer_OpeningFcn, ...
                   'gui_OutputFcn',  @InvivoCytometer_OutputFcn, ...
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


% --- Executes just before InvivoCytometer is made visible.
function InvivoCytometer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InvivoCytometer (see VARARGIN)

    % Choose default command line output for InvivoCytometer
    handles.output = hObject;

    handles.parameters.flagParallelize = logical( get(handles.checkboxParallelize, 'Value') );
    handles.parameters.flagDebugMode = logical( get(handles.checkboxParallelize, 'Value') );
    
    % open matlab pool for parallel processing    
    handles.parameters.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.parameters.flagParallelize && ~handles.parameters.flagPoolOpenedAlready 
        matlabpool open;
    end

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes InvivoCytometer wait for user response (see UIRESUME)
    % uiwait(handles.InvivoCytometer);

% --- Outputs from this function are returned to the command line.
function varargout = InvivoCytometer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnNucleiSegmentation.
function btnNucleiSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to btnNucleiSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    CellSegmentationQualityAnnotator( 'flagParallelize', handles.parameters.flagParallelize, ...
                                      'flagDebugMode',  handles.parameters.flagDebugMode ); 

% --- Executes on button press in btnCellCycleAnalysis.
function btnCellCycleAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to btnCellCycleAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    CellStateAnalyzer( 'flagParallelize', handles.parameters.flagParallelize, ...
                       'flagDebugMode',  handles.parameters.flagDebugMode ); 

% --- Executes on button press in btnAnnotateDataForRegionMerging.
function btnAnnotateDataForRegionMerging_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnnotateDataForRegionMerging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    CellSegmentationQualityAnnotator( 'mode', 'training', ...
                                      'flagParallelize', handles.parameters.flagParallelize, ...
                                      'flagDebugMode',  handles.parameters.flagDebugMode ); 

% --- Executes on button press in btnBuildRegionMergingModel.
function btnBuildRegionMergingModel_Callback(hObject, eventdata, handles)
% hObject    handle to btnBuildRegionMergingModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    trainRegionMergingModel_3D
    
% --- Executes on button press in btnAnnotateDataForCellCycleStateIdentification.
function btnAnnotateDataForCellCycleStateIdentification_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnnotateDataForCellCycleStateIdentification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    CellPatternAnnotator( 'flagParallelize', handles.parameters.flagParallelize, ...
                          'flagDebugMode',  handles.parameters.flagDebugMode ); 
    
% --- Executes on button press in btnBuildCellStateIdentificationModel.
function btnBuildCellStateIdentificationModel_Callback(hObject, eventdata, handles)
% hObject    handle to btnBuildCellStateIdentificationModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    trainCellStateClassificationModel_3D

% --- Executes during object deletion, before destroying properties.
function InvivoCytometer_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to InvivoCytometer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % close matlab pool
    if ~handles.parameters.flagPoolOpenedAlready && matlabpool( 'size' ) > 0
        matlabpool close;
    end

% --- Executes on button press in checkboxParallelize.
function checkboxParallelize_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxParallelize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxParallelize
    
    handles.parameters.flagParallelize = logical( get(hObject,'Value') );
    
    if handles.parameters.flagParallelize
        if matlabpool('size') == 0 
            matlabpool open;
        end
    end
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in checkboxDebugMode.
function checkboxDebugMode_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDebugMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDebugMode
    handles.parameters.flagDebugMode = logical( get(hObject,'Value') );

    % Update handles structure
    guidata(hObject, handles);
    
