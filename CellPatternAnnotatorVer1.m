function varargout = CellPatternAnnotatorVer1(varargin)
%CELLPATTERNANNOTATORVER1 M-file for CellPatternAnnotatorVer1.fig
%      CELLPATTERNANNOTATORVER1, by itself, creates a new CELLPATTERNANNOTATORVER1 or raises the existing
%      singleton*.
%
%      H = CELLPATTERNANNOTATORVER1 returns the handle to a new CELLPATTERNANNOTATORVER1 or the handle to
%      the existing singleton*.
%
%      CELLPATTERNANNOTATORVER1('Property','Value',...) creates a new CELLPATTERNANNOTATORVER1 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CellPatternAnnotatorVer1_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CELLPATTERNANNOTATORVER1('CALLBACK') and CELLPATTERNANNOTATORVER1('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CELLPATTERNANNOTATORVER1.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellPatternAnnotatorVer1

% Last Modified by GUIDE v2.5 28-Sep-2012 18:12:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellPatternAnnotatorVer1_OpeningFcn, ...
                   'gui_OutputFcn',  @CellPatternAnnotatorVer1_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before CellPatternAnnotatorVer1 is made visible.
function CellPatternAnnotatorVer1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

    % Choose default command line output for CellPatternAnnotatorVer1
    handles.output = hObject;

    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
    if exist( historyFile, 'file' )
       handles.data.history = load( historyFile );
    else
       handles.data.history.lastAnalyzedDir = pathstr;
    end

    % Initialize global variables
    handles.data.cellDisplaySize = [80, 80]; % [width, height]
    
    handles.data.flagDataLoaded = false;
    handles.data.flagUseLOG = true;
    handles.data.flagShowGlobalSegMask = true;
    handles.data.flagShowLocalSegMask = true;
    handles.data.flagShowCellBBox = true;

    % open matlab pool for parallel processing
    handles.data.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if ~handles.data.flagPoolOpenedAlready 
        matlabpool open;
    end            
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes CellPatternAnnotatorVer1 wait for user response (see UIRESUME)
    % uiwait(handles.CellPatternAnnotatorVer1);


% --- Outputs from this function are returned to the command line.
function varargout = CellPatternAnnotatorVer1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


% --- Executes on selection change in ListboxCellPatternSelector.
function ListboxCellPatternSelector_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxCellPatternSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns ListboxCellPatternSelector contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from ListboxCellPatternSelector


% --- Executes during object creation, after setting all properties.
function ListboxCellPatternSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListboxCellPatternSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in ShowPreviousCell.
function ShowPreviousCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowPreviousCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end
    
    % decrement cell id
    handles.data.curCellId = max(1, handles.data.curCellId - 1);

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
% --- Executes on button press in ShowNextCell.
function ShowNextCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNextCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end

    % increment cell id
    handles.data.curCellId = min( numel(handles.data.cellStats), handles.data.curCellId + 1);

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


function CellCountDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to CellCountDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellCountDisplay as text
%        str2double(get(hObject,'String')) returns contents of CellCountDisplay as a double


% --- Executes during object creation, after setting all properties.
function CellCountDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellCountDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --------------------------------------------------------------------
function File_Open_Callback(hObject, eventdata, handles)
% hObject    handle to File_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the use to select the oif file for the nucleus channel
    [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastAnalyzedDir, '*.oif' ), 'Select the data file for the nucleus channel' );   
    
    if ~fileName 
        return;
    end
    
    handles.data.dataFilePath = fullfile( pathName, fileName );
    handles.data.history.lastAnalyzedDir = pathName;

    % ask the use to select the oif file for the red/green channel
    [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastAnalyzedDir, '*.oif' ), 'Select the data file for the red/green channel' );   
    
    if ~fileName 
        handles.data.flagRedGreenChannel = false;
    else
        handles.data.flagRedGreenChannel = true;
    end
    
    % load nucleus channel data
    PrettryPrintStepDescription( 'Loading Nuclear Channel Data' );
    imageSeries = loadIntravitalDataset( handles.data.dataFilePath  );    
    
    if (imageSeries.metadata.voxelSpacing(1)/ imageSeries.metadata.voxelSpacing(3)) >= 1
        imageSeries.metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
    end
    
    handles.data.flagDataLoaded = true;
    handles.data.imageData{1} = imageSeries(1).imageData{1,1};
    handles.data.imDisplayRange = [ min(handles.data.imageData{1}(:)), max(handles.data.imageData{1}(:))];    
    
    handles.data.imageDataLOG{1} = ComputeImageLogTransform( handles.data.imageData{1} );
    handles.data.imLogDisplayRange = [ min(handles.data.imageDataLOG{1}(:)), max(handles.data.imageDataLOG{1}(:))];    

    handles.data.metadata = imageSeries(1).metadata;

    % load red/green channel data
    if handles.data.flagRedGreenChannel
        
        PrettryPrintStepDescription( 'Loading Red Green Channel Data' );
        imageSeries = loadIntravitalDataset( handles.data.dataFilePath  );    

        if (imageSeries.metadata.voxelSpacing(1)/ imageSeries.metadata.voxelSpacing(3)) >= 1
            imageSeries.metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
        end

        for i = 1:2
            handles.data.imageData{i+1} = imageSeries(1).imageData{1,i};                
            handles.data.imDisplayRange(i+1,:) = [ min(handles.data.imageData{i+1}(:)), max(handles.data.imageData{i+1}(:))];   
            
            handles.data.imageDataLOG{i+1} = ComputeImageLogTransform( handles.data.imageData{i+1} );
            handles.data.imLogDisplayRange(i+1,:) = [ min(handles.data.imageDataLOG{i+1}(:)), max(handles.data.imageDataLOG{i+1}(:))];   
        end
        
    end
    
    % segment cells
    PrettryPrintStepDescription( 'Running Cell Sementation Algorithm' );
    handles.data.imLabelCellSeg = segmentCellsInIntravitalData( handles.data.imageData{1}, ...
                                                                handles.data.metadata.voxelSpacing, true );

    [handles.data.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    
    % compute properties of each cell
    handles.data.cellStats = regionprops( handles.data.imLabelCellSeg, ...
                                          'Centroid', 'BoundingBox' );
    
    % set current cell id to first cell
    handles.data.curCellId = 1;
    
    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --------------------------------------------------------------------
function UpdateCellDisplay(handles)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end

    curCellStats = handles.data.cellStats( handles.data.curCellId );
    curCellCentroid = round( curCellStats.Centroid );    
    
    % display global MIPs    
    if handles.data.flagUseLOG               
        imGlobalXY = handles.data.imageDataLOG{1}( :, :, curCellCentroid(3) );       
        imGlobalXZ = squeeze( handles.data.imageDataLOG{1}( curCellCentroid(2), :, : ) );
        displayrange = handles.data.imLogDisplayRange;
    else        
        imGlobalXY = handles.data.imageData{1}( :, :, curCellCentroid(3) );
        imGlobalXZ = squeeze( handles.data.imageData{1}( curCellCentroid(2), :, : ) );
        displayrange = handles.data.imDisplayRange;
    end
    imGlobalXY = mat2gray(imGlobalXY, displayrange(1,:) );
    imGlobalXZ = mat2gray(imGlobalXZ', displayrange(1,:) );

    if handles.data.flagShowGlobalSegMask 
        
        imGlobalXYSegMaskRGB = squeeze( handles.data.imCellSegRGBMask( :, :, curCellCentroid(3), : ) );
        imGlobalXZSegMaskRGB = squeeze( handles.data.imCellSegRGBMask(curCellCentroid(2), :, :, :) );
        imGlobalXZSegMaskRGB = permute(imGlobalXZSegMaskRGB, [2,1,3] );
        
        image( genImageRGBMaskOverlay( imGlobalXY, imGlobalXYSegMaskRGB, 0.2 ), ...
               'Parent', handles.Axes_Global_XY );
        image( genImageRGBMaskOverlay( imGlobalXZ, imGlobalXZSegMaskRGB, 0.2 ), ....
               'Parent', handles.Axes_Global_XZ );
           
    else
        
        image( repmat( imGlobalXY, [1,1,3] ), ...
               'Parent', handles.Axes_Global_XY );
        image( repmat( imGlobalXZ, [1,1,3] ), ....
               'Parent', handles.Axes_Global_XZ );
        
    end

    % set axes aspect ratio 
    set( handles.Axes_Global_XY, ...
         'XTickLabel', [], 'YTickLabel', [], ...
         'DataAspectRatio', [ handles.data.metadata.voxelSpacing([2,1]), 1 ], ...
         'PlotBoxAspectRatio', [ handles.data.metadata.voxelSpacing([2,1]), 1 ] );
     
    set( handles.Axes_Global_XZ, ...
         'XTickLabel', [], 'YTickLabel', [], ...
         'DataAspectRatio', [ handles.data.metadata.voxelSpacing([3,1]), 1], ...
         'PlotBoxAspectRatio', [ handles.data.metadata.voxelSpacing([3,1]), 1] );
    
    % draw bounding box around each cell
    if handles.data.flagShowCellBBox
        
        imsize = size( handles.data.imageData{1} );
        hold( handles.Axes_Global_XY, 'on' );
            w = handles.data.cellDisplaySize(1:2);
            ptCorner = curCellCentroid(1:2) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];        
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(1), 2 ) = imsize(1);
            plot( handles.Axes_Global_XY, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
        hold( handles.Axes_Global_XY, 'off' );    

        hold( handles.Axes_Global_XZ, 'on' );
            w = [handles.data.cellDisplaySize(1), curCellStats.BoundingBox(6)];
            ptCorner = curCellCentroid([1,3]) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(3), 2 ) = imsize(3);
            plot( handles.Axes_Global_XZ, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
        hold( handles.Axes_Global_XZ, 'off' );    
        
    end
        
    % extract image within a bounding box around the cell
    subinds = cell(1,3);
    for i = 1:2
        sxi = handles.data.cellDisplaySize(i);
        xi = curCellCentroid(i) - ceil(0.5 * sxi);
        subinds{i} = xi:(xi+sxi-1);
    end
    subinds{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));
    
    [X, Y, Z] = meshgrid( subinds{:} );
    
    imCellCropped = zeros( [handles.data.cellDisplaySize, curCellStats.BoundingBox(6)] );    
    imCellCropped(:) = interpn( handles.data.imageData{1}, Y(:), X(:), Z(:), 'nearest', 0 );

    imCellSegCropped = zeros( [handles.data.cellDisplaySize, curCellStats.BoundingBox(6)] );
    imCellSegCropped(:) = interpn( handles.data.imLabelCellSeg, Y(:), X(:), Z(:), 'nearest', 0 );
    imCellSegCropped = double( imCellSegCropped == handles.data.curCellId );
    
    imCellCropped( ~imCellSegCropped ) = 0;
    
    % display Local XY MIPs
    imCurCellMipXY = mat2gray( max( imCellCropped, [], 3 ) );

    if handles.data.flagShowLocalSegMask 
        
        imCurCellSegMaskMipXY = max( imCellSegCropped, [], 3 );

        curCellColor = handles.data.CellSegColorMap( handles.data.curCellId, : );
        image( genImageMaskOverlay( imCurCellMipXY, imCurCellSegMaskMipXY, curCellColor, 0.2 ), ...
               'Parent', handles.Axes_Nucleus_XY );

        image( genImageMaskOverlay( imCurCellMipXY, imCurCellSegMaskMipXY, curCellColor, 0.2 ), ... 
               'Parent', handles.Axes_RedGreen_XY );
           
    else
        
        image( repmat( imCurCellMipXY, [1,1,3] ), ...
               'Parent', handles.Axes_Nucleus_XY );

        image( repmat( imCurCellMipXY, [1,1,3] ), ... 
               'Parent', handles.Axes_RedGreen_XY );
        
    end

    clear X Y Z;    
       
    set(  handles.Axes_Nucleus_XY, 'XTickLabel', [], 'YTickLabel', [] );
    set( handles.Axes_RedGreen_XY, 'XTickLabel', [], 'YTickLabel', [] );

    % display cell id    
    strtmp = sprintf('%d / %d', handles.data.curCellId, numel(handles.data.cellStats) );
    set(handles.CellCountDisplay, 'String', strtmp);    

% --------------------------------------------------------------------    
function [ imMaskOverlay ] = genImageMaskOverlay( im, mask, maskColor, maskAlpha )

    imr = mat2gray( im );
    img = imr;
    imb = imr;
    mask = logical(mask);
    
    imr(mask) = double( (1 - maskAlpha) * imr(mask) + maskAlpha * maskColor(1) );
    img(mask) = double( (1 - maskAlpha) * img(mask) + maskAlpha * maskColor(2) );
    imb(mask) = double( (1 - maskAlpha) * imb(mask) + maskAlpha * maskColor(3) );
    
    imMaskOverlay = cat(3, imr, img, imb );
    imMaskOverlay( imMaskOverlay > 1 ) = 1;

function [ imMaskOverlay ] = genImageRGBMaskOverlay( im, rgbMask, maskAlpha )

    imMaskOverlay = repmat( mat2gray(im), [1,1,3] );
    blnMask = repmat( max( rgbMask, [], 3 ) > 0, [1, 1, 3] );
    imMaskOverlay(blnMask) = (1 - maskAlpha) * imMaskOverlay(blnMask) + maskAlpha * rgbMask(blnMask);
    imMaskOverlay( imMaskOverlay > 1 ) = 1;
    
% --------------------------------------------------------------------
function File_SaveAnnotation_Callback(~, eventdata, handles)
% hObject    handle to File_SaveAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Close_Callback(~, eventdata, handles)
% hObject    handle to File_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    close(gcf);
    
% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to file_close CellPatternAnnotatorVer1.
function CellPatternAnnotator_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotatorVer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
    
    % file_close
    delete(hObject);
    
function PrettryPrintStepDescription( strStepDescription )

    strStar = strStepDescription;
    strStar(:) = '*';
    strStar = [ '****', strStar, '****' ];
    fprintf( '\n\n%s', strStar );
    fprintf( '\n    %s    ', strStepDescription );
    fprintf( '\n%s\n\n', strStar );


% --- Executes on button press in ShowLastCell.
function ShowLastCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowLastCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end

    % set cell id as last one
    handles.data.curCellId = numel( handles.data.cellStats );

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes on button press in ShowFirstCell.
function ShowFirstCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFirstCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end

    % decrement cell id
    handles.data.curCellId = 1;

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);



function CellDescription_Callback(hObject, eventdata, handles)
% hObject    handle to CellDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellDescription as text
%        str2double(get(hObject,'String')) returns contents of CellDescription as a double


% --- Executes during object creation, after setting all properties.
function CellDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

% --- Executes on button press in CheckboxGlobalSegMask.
function CheckboxGlobalSegMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxGlobalSegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxGlobalSegMask

    % change mask use mode
    handles.data.flagShowGlobalSegMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
% --- Executes on button press in CheckboxGlobalLog.
function CheckboxGlobalLog_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxGlobalLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxGlobalLog

    % change log use mode
    handles.data.flagUseLOG = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    

% --- Executes on button press in CheckboxLocalSegMask.
function CheckboxLocalSegMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxLocalSegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxLocalSegMask

    % change mask use mode
    handles.data.flagShowLocalSegMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


% --------------------------------------------------------------------
function File_Set_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_Set_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function CellPatternAnnotator_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotatorVer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    history = handles.data.history;
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
    save( historyFile, '-struct', 'history' );  

    % close matlab pool
    if handles.data.flagPoolOpenedAlready
        matlabpool close;
    end    

% --- Executes on button press in CheckboxCellBBox.
function CheckboxCellBBox_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxCellBBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxCellBBox

    % change mask use mode
    handles.data.flagShowCellBBox = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
