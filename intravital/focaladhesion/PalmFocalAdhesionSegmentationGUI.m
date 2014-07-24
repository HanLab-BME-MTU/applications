function varargout = PalmFocalAdhesionSegmentationGUI(varargin)
% PALMFOCALADHESIONSEGMENTATIONGUI MATLAB code for PalmFocalAdhesionSegmentationGUI.fig
%      PALMFOCALADHESIONSEGMENTATIONGUI, by itself, creates a new PALMFOCALADHESIONSEGMENTATIONGUI or raises the existing
%      singleton*.
%
%      H = PALMFOCALADHESIONSEGMENTATIONGUI returns the handle to a new PALMFOCALADHESIONSEGMENTATIONGUI or the handle to
%      the existing singleton*.
%
%      PALMFOCALADHESIONSEGMENTATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PALMFOCALADHESIONSEGMENTATIONGUI.M with the given input arguments.
%
%      PALMFOCALADHESIONSEGMENTATIONGUI('Property','Value',...) creates a new PALMFOCALADHESIONSEGMENTATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PalmFocalAdhesionSegmentationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PalmFocalAdhesionSegmentationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PalmFocalAdhesionSegmentationGUI

% Last Modified by GUIDE v2.5 13-Dec-2012 11:36:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PalmFocalAdhesionSegmentationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PalmFocalAdhesionSegmentationGUI_OutputFcn, ...
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


% --- Executes just before PalmFocalAdhesionSegmentationGUI is made visible.
function PalmFocalAdhesionSegmentationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PalmFocalAdhesionSegmentationGUI (see VARARGIN)

    % Choose default command line output for PalmFocalAdhesionSegmentationGUI
    handles.output = hObject;

    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'PalmFASegmentationGUIHistory.mat' );
    if exist( historyFile, 'file' )
       handles.history = load( historyFile );
    else
       handles.history.lastAnalyzedDir = pathstr;
    end

    % initialize global variables
    handles.flagDataLoaded = false;
    handles.flagUseLog = true;
    handles.flagShowFASegMask = true;
    handles.flagShowFABoundingBox = true;
    handles.flagSegmentationDone = false;
    
    set( handles.checkboxUseLog, 'Value', handles.flagUseLog );
    set( handles.checkboxShowFASegMask, 'Value', handles.flagShowFASegMask );
    set( handles.checkboxShowFABoundingBox, 'Value', handles.flagShowFABoundingBox );
    
    % initialize values of uicontrols
    
        % local molecule density
        handles.data.localMoleculeDensityThresholdRange = [0, 100];
        handles.data.localMoleculeDensityThreshold = 2.0;
        
        set( handles.editMinLocalMoleculeDensityThreshold, 'String', num2str(handles.data.localMoleculeDensityThresholdRange(1)) );
        set( handles.editMaxLocalMoleculeDensityThreshold, 'String', num2str(handles.data.localMoleculeDensityThresholdRange(2)) );
        set( handles.editCurLocalMoleculeDensityThreshold, 'String', num2str(handles.data.localMoleculeDensityThreshold) );
        
        set( handles.sliderLocalMoleculeDensity, 'Min',  handles.data.localMoleculeDensityThresholdRange(1) );
        set( handles.sliderLocalMoleculeDensity, 'Max',  handles.data.localMoleculeDensityThresholdRange(2) );
        set( handles.sliderLocalMoleculeDensity, 'Value', handles.data.localMoleculeDensityThreshold );   

        % adhesion area
        handles.data.adhesionAreaThreshold = 0.0;
        
        set( handles.editMinAdhesionAreaThreshold, 'String', num2str(0.0) );
        set( handles.editMaxAdhesionAreaThreshold, 'String', num2str(400.0) );
        set( handles.editCurAdhesionAreaThreshold, 'String', handles.data.adhesionAreaThreshold );
        
        set( handles.sliderAdhesionAreaThreshold, 'Min',  0.0 );
        set( handles.sliderAdhesionAreaThreshold, 'Max',  400.0 );
        set( handles.sliderAdhesionAreaThreshold, 'Val',  0.0 );
    
    % Update handles structure
    guidata(hObject, handles);

% UIWAIT makes PalmFocalAdhesionSegmentationGUI wait for user response (see UIRESUME)
% uiwait(handles.figPalmFocalAdhesionSegmentationGUI);

% --- Outputs from this function are returned to the command line.
function varargout = PalmFocalAdhesionSegmentationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderLocalMoleculeDensity_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLocalMoleculeDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    if ~handles.flagDataLoaded
        return;
    end
    
    handles.data.localMoleculeDensityThreshold = get(hObject,'Value');
    set( handles.editCurLocalMoleculeDensityThreshold, 'String', num2str( handles.data.localMoleculeDensityThreshold ) );    
    
    % do preprocessing
    handles.data.imNMapPreprocessed = PrepocessNMap( handles.data.imNMap,  handles.data.localMoleculeDensityThreshold );
        
    handles.flagSegmentationDone = false;
    
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sliderLocalMoleculeDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLocalMoleculeDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderAdhesionAreaThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    if ~handles.flagDataLoaded
        return;
    end
    
    if ~handles.flagSegmentationDone
        return;
    end

    handles.data.adhesionAreaThreshold = get(hObject,'Value');
    set( handles.editCurAdhesionAreaThreshold, 'String', num2str( handles.data.adhesionAreaThreshold ) );    

    % re-do post processing
    [ handles ] = PostprocessFASegmentationResult( handles );
    
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderAdhesionAreaThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in buttonExecuteSegmentationPipeline.
function buttonExecuteSegmentationPipeline_Callback(hObject, eventdata, handles)
% hObject    handle to buttonExecuteSegmentationPipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Segment focal adhesions using mean-shift clustering
    [ handles ] = SegmentFAUsingMeanShift( handles );
    
    % Post-process segmentation result
    [ handles ] = PostprocessFASegmentationResult( handles );
    
    % Update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);

% --------------------------------------------------------------------
function File_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Load_Palm_Data_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Palm_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the use to select the mat file containing palm images
    [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedDir, '*.mat' ), ...
                                     'Select iPalm data file' );   

    if ~fileName 
        return;
    end
    
    % load the data
    try 
        dataFilePath = fullfile( pathName, fileName );
        dataFileContents = load( dataFilePath );
        imNMap = double( dataFileContents.imNMap );
        
        if size(imNMap,1) < size(imNMap,2)
            imNMap = imrotate(imNMap, 90);
        end
        
    catch err
        return;
    end
    
    handles.flagDataLoaded = true;
    handles.data.dataFilePath = dataFilePath;
    handles.history.lastAnalyzedDir = pathName;
    
    handles.data.imNMap = imNMap;    
    handles.data.displayRangeNMap = ComputeImageDynamicRange(handles.data.imNMap, 98.0);   
    
    handles.data.imNMapLog = ComputeImageLogTransformForDisplay( handles.data.imNMap );
    handles.data.displayRangeNMapLog = [ min(handles.data.imNMapLog(:)), max(handles.data.imNMapLog(:)) ];   

    handles.flagSegmentationDone = false;
    
    % do preprocessing
    handles.data.imNMapPreprocessed = PrepocessNMap( handles.data.imNMap,  handles.data.localMoleculeDensityThreshold );
        
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);

    
% --------------------------------------------------------------------    
function [ imNMapPreprocessed ] = PrepocessNMap( imNMap, localMoleculeDensityThreshold )
    
    imNMapAverage = imfilter( imNMap, fspecial('average', [5, 5]) );
    imNMapPreprocessed = double( imNMapAverage > localMoleculeDensityThreshold & imNMap > 0 );         
    
% --------------------------------------------------------------------    
function [ handles ] = SegmentFAUsingMeanShift( handles )

    hStatusDialog = waitbar(0, 'Segmenting Focal Adhesions');
    
    % Segment focal adhesions using mean-shift clustering
    handles.data.minClusterDistance = 30.0;
    handles.data.bandwidth = handles.data.minClusterDistance / 2.0;
    handles.data.ptIndFA = find( handles.data.imNMapPreprocessed );
    [yind, xind] = ind2sub( size( handles.data.imNMapPreprocessed ), handles.data.ptIndFA );    
    handles.data.ptAdhesion = [ xind, yind ];
    [handles.data.clusterInfoUnpruned, ...
     handles.data.pointToClusterMapUnpruned] = MeanShiftClustering(handles.data.ptAdhesion, handles.data.bandwidth, ... 
                                                                   'flagDebug', false, ...
                                                                   'kernel', 'gaussian', ...
                                                                   'method', 'optimized', ...
                                                                   'minClusterDistance', handles.data.minClusterDistance, ...
                                                                   'flagUseKDTree', true );        
                                                       

    % prepare segmentation mask
    imAdhesionSegLabel = zeros( size( handles.data.imNMap ) );
    imAdhesionSegLabel(handles.data.ptIndFA) = handles.data.pointToClusterMapUnpruned;
    
    % compute adhesion statistics
    handles.data.adhesionStatsUnpruned = regionprops( imAdhesionSegLabel, ...
                                                      {'Area', 'Centroid', ...
                                                       'MajorAxisLength', 'MinorAxisLength', 'Orientation', ...
                                                       'ConvexHull', 'ConvexArea'} );
                                                       
    close( hStatusDialog );

% --------------------------------------------------------------------    
function [ handles ] = PostprocessFASegmentationResult( handles )
    
    handles.data.pointToClusterMap = handles.data.pointToClusterMapUnpruned;
    handles.data.clusterInfo = handles.data.clusterInfoUnpruned;
    handles.data.adhesionStats = handles.data.adhesionStatsUnpruned;
    
    % perform post-processing -- delete small clusters
    if handles.data.adhesionAreaThreshold > 0.0 
        
        numSignificantClusters = 0;
        flagSmallCluster = false( numel(handles.data.clusterInfo), 1 );
        for i = 1:numel( handles.data.clusterInfo )          
            curAdhesionDensity = handles.data.adhesionStats(i).Area / handles.data.adhesionStats(i).ConvexArea;
            if handles.data.adhesionStats(i).ConvexArea < handles.data.adhesionAreaThreshold || curAdhesionDensity < 0.2            
                flagSmallCluster(i) = true;
                handles.data.pointToClusterMap( handles.data.clusterInfo(i).ptIdData ) = 0;
            else
                numSignificantClusters = numSignificantClusters + 1;
                handles.data.pointToClusterMap( handles.data.clusterInfo(i).ptIdData ) = numSignificantClusters;
            end
        end
        handles.data.clusterInfo( flagSmallCluster ) = []; % delete small clusters
    end
    
    % update segmentation mask
    handles.data.imAdhesionSegLabel = zeros( size( handles.data.imNMap ) );
    handles.data.imAdhesionSegLabel(handles.data.ptIndFA) = handles.data.pointToClusterMap;
    [handles.data.imAdhesionSegMaskRGB, ...
     handles.data.adhesionLabelToColorMap] = label2rgbND( handles.data.imAdhesionSegLabel ); 
    
    % compute adhesion statistics
    handles.data.adhesionStats = regionprops( handles.data.imAdhesionSegLabel, ...
                                              {'Area', 'Centroid', ...
                                               'MajorAxisLength', 'MinorAxisLength', 'Orientation', ...
                                               'ConvexHull', 'ConvexArea'} );
        
    handles.flagSegmentationDone = true;
    
% --------------------------------------------------------------------
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 99.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
% --------------------------------------------------------------------
function UpdateDisplay(handles)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    % display NMap
    if handles.flagUseLog
        imCurNMap = mat2gray( handles.data.imNMapLog, handles.data.displayRangeNMapLog );        
    else        
        imCurNMap = mat2gray( handles.data.imNMap, handles.data.displayRangeNMap );        
    end    
    
    cla( handles.Axes_PalmImage );
    image( repmat(imCurNMap, [1,1,3]), 'Parent', handles.Axes_PalmImage );
    set( handles.Axes_PalmImage, 'XTickLabel', [], 'YTickLabel', [] );
    
    % display pre-processing result
    cla( handles.Axes_Preprocessing );
    image( repmat(handles.data.imNMapPreprocessed, [1, 1, 3]), 'Parent', handles.Axes_Preprocessing );
    set( handles.Axes_Preprocessing, 'XTickLabel', [], 'YTickLabel', [] );
    
    % display segmentation result
    cla( handles.Axes_SegmentationResult );
    
    if ~handles.flagSegmentationDone        
        return;
    end
    
    if handles.flagShowFASegMask
        image( genImageRGBMaskOverlay( imCurNMap, handles.data.imAdhesionSegMaskRGB, 0.4 ), ...
               'Parent', handles.Axes_SegmentationResult );
    else
        image( repmat(imCurNMap, [1,1,3]), ...
               'Parent', handles.Axes_SegmentationResult );
    end
    
    if handles.flagShowFABoundingBox
    
        hold on;
        %hold( handles.Axes_SegmentationResult, 'on' );
        
            for i = 1:max(handles.data.imAdhesionSegLabel(:))
               ptCenter = handles.data.adhesionStats(i).Centroid;
               a = 0.5 * handles.data.adhesionStats(i).MajorAxisLength;
               b = 0.5 * handles.data.adhesionStats(i).MinorAxisLength;

               theta = handles.data.adhesionStats(i).Orientation;
               ptRect = ([ -a, -b; ...
                            a, -b; ...
                            a, b; ...
                           -a, b; ...
                           -a, -b ])';
               ptRect(end+1,:) = 1;
               transMat = [ cosd(theta), sind(theta), ptCenter(1); ...
                           -sind(theta) cosd(theta), ptCenter(2); ...
                            0, 0, 1 ];
               ptRect = round( transMat * ptRect );
               
               plot( ptRect(1,:), ptRect(2,:), '-', ...
                     'Color', handles.data.adhesionLabelToColorMap(i,:), ...
                     'LineWidth', 2.0 ); 
                 
               plot( ptCenter(1), ptCenter(2), 'o', ...
                     'Color', [0,0,0], 'MarkerFaceColor', ...
                     handles.data.adhesionLabelToColorMap(i,:), 'LineWidth', 2.0 ); 
               
            end
        
       hold off;
       %hold( handles.Axes_SegmentationResult, 'off' );
        
    end
    
    set( handles.Axes_SegmentationResult, 'XTickLabel', [], 'YTickLabel', [] );
        
% --------------------------------------------------------------------
function File_Close_Callback(hObject, eventdata, handles)
% hObject    handle to File_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    close(gcf);

function editMinAdhesionAreaThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editMinAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinAdhesionAreaThreshold as text
%        str2double(get(hObject,'String')) returns contents of editMinAdhesionAreaThreshold as a double


% --- Executes during object creation, after setting all properties.
function editMinAdhesionAreaThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxAdhesionAreaThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxAdhesionAreaThreshold as text
%        str2double(get(hObject,'String')) returns contents of editMaxAdhesionAreaThreshold as a double


% --- Executes during object creation, after setting all properties.
function editMaxAdhesionAreaThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCurAdhesionAreaThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editCurAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurAdhesionAreaThreshold as text
%        str2double(get(hObject,'String')) returns contents of editCurAdhesionAreaThreshold as a double


% --- Executes during object creation, after setting all properties.
function editCurAdhesionAreaThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCurAdhesionAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinLocalMoleculeDensityThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editMinLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinLocalMoleculeDensityThreshold as text
%        str2double(get(hObject,'String')) returns contents of editMinLocalMoleculeDensityThreshold as a double


% --- Executes during object creation, after setting all properties.
function editMinLocalMoleculeDensityThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxLocalMoleculeDensityThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxLocalMoleculeDensityThreshold as text
%        str2double(get(hObject,'String')) returns contents of editMaxLocalMoleculeDensityThreshold as a double


% --- Executes during object creation, after setting all properties.
function editMaxLocalMoleculeDensityThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCurLocalMoleculeDensityThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editCurLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurLocalMoleculeDensityThreshold as text
%        str2double(get(hObject,'String')) returns contents of editCurLocalMoleculeDensityThreshold as a double


% --- Executes during object creation, after setting all properties.
function editCurLocalMoleculeDensityThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCurLocalMoleculeDensityThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxUseLog.
function checkboxUseLog_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUseLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUseLog

    handles.flagUseLog = get(hObject,'Value');
    
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figPalmFocalAdhesionSegmentationGUI_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figPalmFocalAdhesionSegmentationGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    history = handles.history;
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'PalmFASegmentationGUIHistory.mat' );
    save( historyFile, '-struct', 'history' );  


% --- Executes on button press in checkboxShowFASegMask.
function checkboxShowFASegMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowFASegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowFASegMask

    handles.flagShowFASegMask = get(hObject,'Value');
    
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in checkboxShowFABoundingBox.
function checkboxShowFABoundingBox_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowFABoundingBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowFABoundingBox

    handles.flagShowFABoundingBox = get(hObject,'Value');
    
    % update display
    UpdateDisplay(handles);

    % Update handles structure
    guidata(hObject, handles);
