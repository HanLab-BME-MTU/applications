function varargout = CellSegmentationQualityAnnotator(varargin)
%CELLSEGMENTATIONQUALITYANNOTATOR M-file for CellSegmentationQualityAnnotator.fig
%      CELLSEGMENTATIONQUALITYANNOTATOR, by itself, creates a new CELLSEGMENTATIONQUALITYANNOTATOR or raises the existing
%      singleton*.
%
%      H = CELLSEGMENTATIONQUALITYANNOTATOR returns the handle to a new CELLSEGMENTATIONQUALITYANNOTATOR or the handle to
%      the existing singleton*.
%
%      CELLSEGMENTATIONQUALITYANNOTATOR('Property','Value',...) creates a new CELLSEGMENTATIONQUALITYANNOTATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CellSegmentationQualityAnnotator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CELLSEGMENTATIONQUALITYANNOTATOR('CALLBACK') and CELLSEGMENTATIONQUALITYANNOTATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CELLSEGMENTATIONQUALITYANNOTATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellSegmentationQualityAnnotator

% Last Modified by GUIDE v2.5 06-Sep-2014 01:35:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellSegmentationQualityAnnotator_OpeningFcn, ...
                   'gui_OutputFcn',  @CellSegmentationQualityAnnotator_OutputFcn, ...
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

% --- Executes just before CellSegmentationQualityAnnotator is made visible.
function CellSegmentationQualityAnnotator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

    p = inputParser;
    p.addParamValue('mode', 'analysis', @(x) (ischar(x) && ismember(x, {'training', 'analysis'})) );
    p.addParamValue('flagParallelize', false, @isscalar);
    p.addParamValue('flagDebugMode', true, @isscalar);
    p.parse( varargin{:} );    
    PARAMETERS = p.Results;

    % Choose default command line output for CellSegmentationQualityAnnotator
    handles.output = hObject;

    % make sure weka is in the path    
    AddWekaClassesToPath()
    
    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellSegmentationQualityAnnotatorHistory.mat' );
    if exist( historyFile, 'file' )
       handles.history = load( historyFile );
    else
       handles.history.lastAnalyzedNucleusDir = pathstr;
    end

    % Initialize global variables
    handles.cellDisplaySize = 70; 
    
    handles.flagDataLoaded = false;
    handles.flagUseLOG = get(handles.CheckboxGlobalLog, 'Value');
    handles.flagShowGlobalSegMask = get(handles.CheckboxGlobalSegMask, 'Value');
    handles.flagShowCellBBox = get(handles.CheckboxCellBBox, 'Value');

    handles.defaultCellPatternTypes = cellstr( get(handles.ListboxSegmentationQualitySelector, 'String') );
    
    set( handles.poplistCellClassSelector, 'String', handles.defaultCellPatternTypes );
    handles.data.cellPatternTypes = handles.defaultCellPatternTypes;
    
    handles.imarisAppCellSeg = [];
    handles.imarisAppCellPattern = [];
    handles.imarisAppCellSegCropped = [];
    
    % default parameters
    handles.defaultParameters.flagParallelize = PARAMETERS.flagParallelize;
    handles.defaultParameters.flagDebugMode = PARAMETERS.flagDebugMode;
    
    if strcmp(PARAMETERS.mode, 'analysis')
        f = rdir( fullfile(fileparts(mfilename('fullpath')), '**', 'regionMerging.model') );
        if ~isempty(f)
            handles.defaultParameters.nucleiSegmentation = NucleiSegmentationParametersGUI( 'defaultWithRegionMerging' );
            handles.defaultParameters.nucleiSegmentation.flagPerformRegionMerging = true;
            handles.defaultParameters.nucleiSegmentation.regionMergingModelFile = f(1).name;
        else
            handles.defaultParameters.nucleiSegmentation = NucleiSegmentationParametersGUI( 'defaultWithoutRegionMerging' );
        end
    else
        handles.defaultParameters.nucleiSegmentation = NucleiSegmentationParametersGUI( 'defaultRegionMergingTraining' );
    end
    
    handles.parameters = handles.defaultParameters;
    
    % open matlab pool for parallel processing    
    handles.parameters.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.parameters.flagParallelize && ~handles.parameters.flagPoolOpenedAlready 
        matlabpool open;
    end
    
    % Set callbacks
    set(handles.CellSegmentationQualityAnnotator, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes CellSegmentationQualityAnnotator wait for user response (see UIRESUME)
    % uiwait(handles.CellSegmentationQualityAnnotator);

% --- Outputs from this function are returned to the command line.
function varargout = CellSegmentationQualityAnnotator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

% --- Executes on selection change in ListboxSegmentationQualitySelector.
function ListboxSegmentationQualitySelector_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxSegmentationQualitySelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns ListboxSegmentationQualitySelector contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from ListboxSegmentationQualitySelector

    curCellId = handles.dataDisplay.curCellId;
    curCellStats = handles.data.cellStats( curCellId );
    
    % assigned selected pattern type to current cell    
    selPatternId = get(hObject, 'Value');
    curCellStats.cellPatternId = selPatternId;
    curCellStats.cellPatternType = handles.data.cellPatternTypes{selPatternId};
    
    handles.data.cellStats(curCellId) = curCellStats;
    
    % if pattern type is not equal to 'None' remove it from the list of
    % unannotated cells, if not add it to the list
    if ~strcmpi( handles.data.cellPatternTypes{selPatternId}, 'None' )
        handles.data.unannotatedCellList = setdiff( handles.data.unannotatedCellList, curCellId ); 
    else
        handles.data.unannotatedCellList = union( handles.data.unannotatedCellList, curCellId ); 
    end

    if isempty( handles.data.unannotatedCellList )
        set( handles.ShowNextUnannotatedCell, 'Enable', 'off' );
    else
        set( handles.ShowNextUnannotatedCell, 'Enable', 'on' );
    end
    
    % navigate to next cell automatically
    handles.dataDisplay.curCellId = min( numel(handles.data.cellStats), handles.dataDisplay.curCellId + 1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function ListboxSegmentationQualitySelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListboxSegmentationQualitySelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


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
function File_Load_Image_Data_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Image_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % ask the user to select the oif file for the nucleus channel
    status = bfCheckJavaPath(1);
    [fileName,pathName] = uigetfile(bfGetFileExtensions, ...
                                    'Select the file contaning the nuclear marker channel - histone-2B', ...
                                    handles.history.lastAnalyzedNucleusDir);
    
    if ~fileName 
        return;
    end
    
    dataFilePath = fullfile( pathName, fileName );
    handles.history.lastAnalyzedNucleusDir = pathName;

    % load image data
    PrettyPrintStepDescription( 'Loading Image Data' );

    [imageData, metadata] = uiGetBioImageData(dataFilePath, 'channelDescriptionsAndNamesList', {{'Nuclear Marker Channel', 'channelIdNuclei'}});
    
    if isempty(imageData)
        return;
    end
    
    metadata
        
    % store image data in handles structures    
    clear handles.data;
    
    handles.flagDataLoaded = true;
    handles.data.dataFilePath = dataFilePath;
    handles.data.metadata = metadata;
    handles.data.imageData = imageData;
    
    handles.data.metadata.channelColors = ones(1,3);
    
    set( handles.CellSegmentationQualityAnnotator, 'Name', sprintf( 'Cell Nuclei Segmenter - %s', handles.data.dataFilePath ) );
    
    % pre-compute data needed for display
    handles = ComputeDisplayData( handles );
    
    % Update handles structure
    guidata(hObject, handles);

    % Run analysis
    RunAnalysis(hObject, handles);

% --------------------------------------------------------------------    
function [ handles ] = ComputeDisplayData(handles)

    for i = 1:numel( handles.data.imageData )
        handles.data.imageData{i} = double( handles.data.imageData{i} );
        handles.dataDisplay.imDisplayRange(i,:) = ComputeImageDynamicRange( handles.data.imageData{i}, 98.0 );   
        handles.dataDisplay.imageDataLOG{i} = ComputeImageLogTransformForDisplay( handles.data.imageData{i} );
        handles.dataDisplay.imLogDisplayRange(i,:) = [ min(handles.dataDisplay.imageDataLOG{i}(:)), max(handles.dataDisplay.imageDataLOG{i}(:))];   
    end
    
% --------------------------------------------------------------------    
function RunAnalysis(hObject, handles)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    % segment cells
    PrettyPrintStepDescription( 'Running Cell Sementation Algorithm' );    
    hStatusDialog = waitbar(0, 'Segmentating Cells');
   
    regionMergingModelFile = [];
    if handles.parameters.nucleiSegmentation.flagPerformRegionMerging
        regionMergingModelFile = handles.parameters.nucleiSegmentation.regionMergingModelFile;
    end
    
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints] = segmentCellsInIntravitalData( handles.data.imageData{handles.data.metadata.channelIdNuclei}, ...
                                                                    handles.data.metadata.pixelSize, ...                                                                      
                                                                    'flagParallelize', handles.parameters.flagParallelize, ...
                                                                    'flagDebugMode', handles.parameters.flagDebugMode, ...
                                                                    'cellDiameterRange', handles.parameters.nucleiSegmentation.cellDiameterRange, ...
                                                                    'thresholdingAlgorithm', 'MinErrorPoissonSliceBySliceLocal', ...
                                                                    'minSignalToBackgroundRatio', 3.0, ...
                                                                    'seedPointDetectionAlgorithm', handles.parameters.nucleiSegmentation.seedPointDetectionAlgorithm, ...
                                                                    'minCellVolume', handles.parameters.nucleiSegmentation.minCellVolume, ...
                                                                    'flagIgnoreCellsOnXYBorder', handles.parameters.nucleiSegmentation.flagIgnoreXYBorderCells, ...
                                                                    'regionMergingModelFile', regionMergingModelFile);

%                                                                     'minCellBBoxSizePhysp', [5, 5, 2] .* handles.data.metadata.pixelSize, ...
                                                                
    [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.dataDisplay.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints, ones(3,3,3) );
    closeStatusDialog(hStatusDialog);
    
    % compute properties of each cell
    handles.data.cellStats = ComputeCellProperties( handles );
    
    % initialize annotation to none
    for i = 1:numel(handles.data.cellStats)
        
        handles.data.cellStats(i).cellPatternId = 1; %Not Annotated
        handles.data.cellStats(i).cellPatternType = handles.data.cellPatternTypes{1};
    
    end
        
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % add all cells to unannotated list
    handles.data.unannotatedCellList = 1:numel(handles.data.cellStats);
    set( handles.ShowNextUnannotatedCell, 'Enable', 'on' );
    
    % close progress bar
    closeStatusDialog(hStatusDialog);
    
    % Update handles structure
    guidata(hObject, handles);
   
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

% --------------------------------------------------------------------
function [cellStats] = ComputeCellProperties( handles )

    PrettyPrintStepDescription( 'Computing Properties of Segmentated Cells' );    
    
    hStatusDialog = waitbar(0, 'Computing properties of segmented cells');
    
    cellStats = regionprops( handles.data.imLabelCellSeg, ...
                             'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
    
    for i = 1:numel(cellStats)

        % intensity descriptors
        cellPixelIntensities = handles.data.imageData{handles.data.metadata.channelIdNuclei}( cellStats(i).PixelIdxList );
        cellStats(i).meanIntensity = mean( cellPixelIntensities );
        cellStats(i).stdIntensity = std( double(cellPixelIntensities) );
        cellStats(i).minIntensity = min( cellPixelIntensities );
        cellStats(i).maxIntensity = max( cellPixelIntensities );
        
        % volume
        cellStats(i).AreaPhysp = cellStats(i).Area * prod( handles.data.metadata.pixelSize );
        
        % update progress
        hStatusDialog = waitbar(i/numel(cellStats), hStatusDialog, 'Computing properties of segmented cells');    

    end
    
    closeStatusDialog(hStatusDialog);    

% --------------------------------------------------------------------
function File_Run_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Run_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    RunAnalysis(hObject, handles);

% --------------------------------------------------------------------
function UpdateCellDisplay(handles)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellSliceId = handles.dataDisplay.curCellSliceId;
    curCellCentroid = round( curCellStats.Centroid );  
    curChannelId = handles.data.metadata.channelIdNuclei;
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );
    
    % display global cross section images  
    if handles.flagUseLOG               
        imGlobalXY = handles.dataDisplay.imageDataLOG{curChannelId}( :, :, curCellSliceId(3) );       
        imGlobalXZ = squeeze( handles.dataDisplay.imageDataLOG{curChannelId}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.dataDisplay.imageDataLOG{curChannelId}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imLogDisplayRange;
    else        
        imGlobalXY = handles.data.imageData{curChannelId}( :, :, curCellSliceId(3) );
        imGlobalXZ = squeeze( handles.data.imageData{curChannelId}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.data.imageData{curChannelId}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imDisplayRange;
    end
    imGlobalXY = mat2gray(imGlobalXY, displayrange(1,:) );
    imGlobalXZ = mat2gray(imGlobalXZ', displayrange(1,:) );
    imGlobalYZ = mat2gray(imGlobalYZ, displayrange(1,:) );

    if handles.flagShowGlobalSegMask 
        
        imGlobalXYSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask( :, :, curCellSliceId(3), : ) );
        
        imGlobalXZSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask(curCellSliceId(1), :, :, :) );
        imGlobalXZSegMaskRGB = permute(imGlobalXZSegMaskRGB, [2,1,3] );
        
        imGlobalYZSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask(:, curCellSliceId(2), :, :) );

        imGlobalXYDisplay = genImageRGBMaskOverlay( imGlobalXY, imGlobalXYSegMaskRGB, 0.2 );
        imGlobalXZDisplay = genImageRGBMaskOverlay( imGlobalXZ, imGlobalXZSegMaskRGB, 0.2 );
        imGlobalYZDisplay = genImageRGBMaskOverlay( imGlobalYZ, imGlobalYZSegMaskRGB, 0.2 );
        
    else
        
        imGlobalXYDisplay = repmat( imGlobalXY, [1,1,3] );
        imGlobalXZDisplay = repmat( imGlobalXZ, [1,1,3] );
        imGlobalYZDisplay = repmat( imGlobalYZ, [1,1,3] );
        
    end
        
    cla( handles.Axes_Global_XY, 'reset' );
    cla( handles.Axes_Global_XZ, 'reset' );
    cla( handles.Axes_Global_YZ, 'reset' );
    
    image( imGlobalXYDisplay, 'Parent', handles.Axes_Global_XY );           
    image( imGlobalXZDisplay, 'Parent', handles.Axes_Global_XZ );
    image( imGlobalYZDisplay, 'Parent', handles.Axes_Global_YZ );           

    set( handles.Axes_Global_XY, 'XTickLabel', [], 'YTickLabel', [] );
    set( handles.Axes_Global_XZ, 'XTickLabel', [], 'YTickLabel', [] );
    set( handles.Axes_Global_YZ, 'XTickLabel', [], 'YTickLabel', [] );        

    % draw bounding box around each cell
    if handles.flagShowCellBBox
        
        imsize = size( handles.data.imageData{curChannelId} );
        hold( handles.Axes_Global_XY, 'on' );
            
            w = curCellDisplaySize([1,1]);
            ptCorner = round(curCellCentroid(1:2) - 0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1)-1, 0; w-1 ; 0, w(2)-1; 0, 0 ];        
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(1), 2 ) = imsize(1);            
            plot( handles.Axes_Global_XY, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
            
            plot( handles.Axes_Global_XY, curCellSliceId(2) * [1, 1], [1, size(imGlobalXY, 1)], 'g-' );
            plot( handles.Axes_Global_XY, [1, size(imGlobalXY, 2)], curCellSliceId(1) * [1, 1], 'g-' );
            
        hold( handles.Axes_Global_XY, 'off' );    

        hold( handles.Axes_Global_XZ, 'on' );
        
            w = [curCellDisplaySize(1), curCellStats.BoundingBox(6)];
            ptCorner = curCellCentroid([1,3]) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(3), 2 ) = imsize(3);
            plot( handles.Axes_Global_XZ, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 

            plot( handles.Axes_Global_XZ, [1, size(imGlobalXZ,2)], curCellSliceId(3) * [1, 1], 'g-' );
            plot( handles.Axes_Global_XZ, curCellSliceId(2) * [1, 1], [1, size(imGlobalXZ,1)], 'g-' );
            
        hold( handles.Axes_Global_XZ, 'off' );    

        hold( handles.Axes_Global_YZ, 'on' );
        
            w = [curCellStats.BoundingBox(6), curCellDisplaySize(1)];
            ptCorner = curCellCentroid([3,2]) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(3), 1 ) = imsize(3);
            ptBBox( ptBBox(:,2) > imsize(1), 2 ) = imsize(1);
            plot( handles.Axes_Global_YZ, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
            
            plot( handles.Axes_Global_YZ, curCellSliceId(3) * [1, 1], [1, size(imGlobalYZ,1)], 'g-' );
            plot( handles.Axes_Global_YZ, [1, size(imGlobalYZ,2)], curCellSliceId(1) * [1, 1], 'g-' );
            
        hold( handles.Axes_Global_XZ, 'off' );    
        
    end
        
    % extract image within a bounding box around the cell
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{curChannelId});
    for i = 1:2
        
        xi = round(curCellCentroid(3-i) - 0.5 * curCellDisplaySize);
        
        xi_low = xi;
        if xi_low < 1 
            xi_low = 1;
        end
        
        xi_high = xi + curCellDisplaySize - 1;
        if xi_high > imsize(i)
            xi_high = imsize(i);
        end
        
        subinds{i} = xi_low:xi_high;
        
    end    
    subinds{3} = curCellSliceId(3);    
    
    if handles.flagUseLOG
        imCellCropped = mat2gray( handles.dataDisplay.imageDataLOG{1}(subinds{:}), handles.dataDisplay.imLogDisplayRange(1,:) );        
    else
        imCellCropped = mat2gray( handles.data.imageData{curChannelId}(subinds{:}), handles.dataDisplay.imDisplayRange(1,:) );        
    end
    imCellSegCropped = handles.data.imLabelCellSeg( subinds{:} );
    imCellSegCropped = double( imCellSegCropped == handles.dataDisplay.curCellId );    
    
    imCellSeedCropped = handles.dataDisplay.imCellSeedPoints( subinds{:} );
    imCellSeedCropped(~imCellSegCropped) = 0;
    imCellSeedCropped(:, :, :, 2:3) = 0;
    
    % display the nuclei channel of the extracted cell image
    if handles.flagShowGlobalSegMask 
        
        curCellColor = handles.data.CellSegColorMap( handles.dataDisplay.curCellId, : );
        imNucleusXYDisplay = genImageMaskOverlay( imCellCropped, {imCellSegCropped, imCellSeedCropped}, [curCellColor; 1 0 0], [0.2, 0.6] );
        
    else
        
        imNucleusXYDisplay = repmat( imCellCropped, [1,1,3] );
        
    end

    cla( handles.Axes_Nucleus_XY, 'reset' );
    image( imNucleusXYDisplay, 'Parent', handles.Axes_Nucleus_XY );
    set(  handles.Axes_Nucleus_XY, 'XTickLabel', [], 'YTickLabel', [] ); 
    
    % display Nucleus MIP
    subindsMIP = subinds;
    subindsMIP{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));

    imCurCellSegMIP = max( double( handles.data.imLabelCellSeg( subindsMIP{:} ) == handles.dataDisplay.curCellId ), [], 3);
    imCurCellMIP = mat2gray( max( double( handles.data.imageData{curChannelId}( subindsMIP{:} ) ), [], 3 ) .* imCurCellSegMIP );
    imHistoneMIPDisplay = repmat( imCurCellMIP, [1,1,3] );
    
    cla( handles.Axes_Histone_MIP, 'reset' );
    image( imHistoneMIPDisplay, 'Parent', handles.Axes_Histone_MIP );
    set( handles.Axes_Histone_MIP, 'XTickLabel', [], 'YTickLabel', [] );
    
    % display cell id    
    strtmp = sprintf('%d / %d', handles.dataDisplay.curCellId, numel(handles.data.cellStats) );
    set(handles.CellCountDisplay, 'String', strtmp);    
    
    % display cell slice id
    strtmp = sprintf('x-slice: %.3d / %.3d\ny-slice: %.3d / %.3d\nz-slice: %.3d / %.3d', ...
                     curCellSliceId(2), imsize(2), ...
                     curCellSliceId(1), imsize(1), ...
                     curCellSliceId(3), imsize(3) );
    set(handles.EditboxSliceId, 'String', strtmp);    
    
% --------------------------------------------------------------------
function UpdateCellDescriptors(handles)

    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );
    curCellCentroid = round(curCellStats.Centroid);
    
    % display cell information
    strCellDescription = sprintf( '' );
    strCellDescription = sprintf( '%s\nVolume (cu um): %.2f', strCellDescription, curCellStats.AreaPhysp );
    strCellDescription = sprintf( '%s\n\nMean-std Intensity: [%.2f, %.2f]', strCellDescription, curCellStats.meanIntensity, curCellStats.stdIntensity );
    strCellDescription = sprintf( '%s\n\nIntensity Range: [%d, %d]', strCellDescription, curCellStats.minIntensity, curCellStats.maxIntensity );
    strCellDescription = sprintf( '%s\n\nCentroid: [%d, %d, %d]', strCellDescription, curCellCentroid(1), curCellCentroid(2), curCellCentroid(3) );

    set(handles.LabelCellDescription, 'String', strCellDescription); 

% --------------------------------------------------------------------    
function [ imMultichannelOverlay ] = genMultiChannelOverlay( im, displaycolors )

    imsize = size( im(:,:,1) );
    displaycolors_hsv = rgb2hsv( displaycolors );
    
    imMultichannelOverlay = zeros( [imsize, 3] );
    
    cur_im_hsv = zeros( [imsize, 3] );
    for i = 1:size(im,3)
        curGrayImage = im(:,:,i);
        cur_im_hsv(:,:,1) = displaycolors_hsv(i,1);
        cur_im_hsv(:,:,2) = displaycolors_hsv(i,2);
        cur_im_hsv(:,:,3) = curGrayImage;  
        cur_im_rgb = hsv2rgb( cur_im_hsv );
        imMultichannelOverlay = imMultichannelOverlay + cur_im_rgb;
        imMultichannelOverlay( imMultichannelOverlay > 1 ) = 1;
    end
    
% --------------------------------------------------------------------    
function [ imMaskOverlay ] = genImageMaskOverlay( im, masks, maskColors, maskAlphas )

    imr = mat2gray( im );
    img = imr;
    imb = imr;
    
    if ~iscell( masks )
        masks = { masks };
    end
       
    for i = 1:numel(masks)
        
        curMask = logical(masks{i});
        curMaskColor = maskColors(i,:);
        curMaskAlpha = maskAlphas(i);
                
        imr(curMask) = double( (1 - curMaskAlpha) * imr(curMask) + curMaskAlpha * curMaskColor(1) );
        img(curMask) = double( (1 - curMaskAlpha) * img(curMask) + curMaskAlpha * curMaskColor(2) );
        imb(curMask) = double( (1 - curMaskAlpha) * imb(curMask) + curMaskAlpha * curMaskColor(3) );
        
    end
    
    imMaskOverlay = cat(3, imr, img, imb );
    imMaskOverlay( imMaskOverlay > 1 ) = 1;

function [ imMaskOverlay ] = genImageRGBMaskOverlay( im, rgbMask, maskAlpha )

    imMaskOverlay = repmat( mat2gray(im), [1,1,3] );
    blnMask = repmat( max( rgbMask, [], 3 ) > 0, [1, 1, 3] );
    imMaskOverlay(blnMask) = (1 - maskAlpha) * imMaskOverlay(blnMask) + maskAlpha * rgbMask(blnMask);
    imMaskOverlay( imMaskOverlay > 1 ) = 1;
    
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 99.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
% --------------------------------------------------------------------
function File_Load_Annotation_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % get annotation file from user
    if isfield( handles.history, 'lastOutputDir' )
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastOutputDir, 'CellSegmentationQualityAnnotation.mat' ), 'Select annotation file' );   
    else
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedNucleusDir, 'CellSegmentationQualityAnnotation.mat' ), 'Select annotation file' );   
    end

    if ~fileName 
        return;
    end
    
    hStatusDialog = waitbar(0, 'Loading selected annotation file ...');
    
    % load annotation data from file
    annotationFile = fullfile(pathName, fileName);
    annotationData = load( annotationFile );
    
    % retrieve data needed for this tool
    handles.data = [];
    
        % basic data
        if iscell(annotationData.dataFilePath)
            handles.data.dataFilePath = annotationData.dataFilePath{1};
        else
            handles.data.dataFilePath = annotationData.dataFilePath;
        end
        
        % metadata
        if ~isfield(annotationData.metadata, 'pixelSize') % check if old metadata struct
        
            metadata.version = '2.0.0';
            metadata.dataFilePath = annotationData.dataFilePath([1, 2, 2]);
            metadata.format = repmat({'Olympus'}, 1, 3);
            
            metadata.numSeries = ones(1,3);
            metadata.seriesId = ones(1,3);
            
            metadata.numTimePoints = 1;
            metadata.timePointId = 1;
            
            metadata.numChannels = 3;
            metadata.channelId = [1, 1, 2];
            metadata.channelNames = {'Nuclei', 'FUCCI Geminin', 'FUCCI Cdt1'};
            
            metadata.imageSize = annotationData.metadata.volSize;
            metadata.pixelSize = annotationData.metadata.voxelSpacing;
            
            metadata.pixelType = 'uint16';
            metadata.bitsPerPixel = 12;
            
            handles.data.metadata = metadata;
            
        else
            handles.data.metadata = annotationData.metadata;
        end
        
        % image data
        handles.data.imageData = annotationData.imageData;

        % basic display data
        handles = ComputeDisplayData(handles);
    
        % segmentation stuff
        handles.data.imLabelCellSeg = annotationData.imLabelCellSeg;        
        handles.data.imCellSeedPoints = annotationData.imCellSeedPoints;        

        if ~isfield( annotationData, 'CellSegColorMap' )
            [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND( handles.data.imLabelCellSeg );
        else
            handles.data.CellSegColorMap = annotationData.CellSegColorMap;        
            handles.dataDisplay.imCellSegRGBMask = label2rgbND( handles.data.imLabelCellSeg, handles.data.CellSegColorMap );
        end
        
        handles.dataDisplay.imCellSeedPoints = imdilate(annotationData.imCellSeedPoints, ones(3,3,3));
        
        % annotated stuff
        handles.data.unannotatedCellList = annotationData.unannotatedCellList;
        handles.data.cellStats = ComputeCellProperties( handles );
        handles.data.cellPatternTypes = handles.defaultCellPatternTypes;

        for i = 1:numel(handles.data.cellStats)

            curCellPatternType = annotationData.cellStats(i).cellPatternType;
            if ~ismember(curCellPatternType, handles.data.cellPatternTypes)                
                curCellPatternType = 'Good_Segmentation';
            end
            curCellPatternId = find( strcmpi(handles.data.cellPatternTypes, curCellPatternType) );
            
            handles.data.cellStats(i).cellPatternType = curCellPatternType;
            handles.data.cellStats(i).cellPatternId = curCellPatternId; 

        end
        
    % data ready for display
    handles.flagDataLoaded = true;
    
    % set cell pattern type list and cell class selector list
    %set( handles.ListboxSegmentationQualitySelector, 'String', handles.data.cellPatternTypes );

    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % change window name
    set( handles.CellSegmentationQualityAnnotator, 'Name', sprintf( 'Cell Nuclei Segmenter - %s', annotationFile ) );
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % close status dialog
    closeStatusDialog(hStatusDialog);    
    
% --------------------------------------------------------------------
function File_SaveAnnotation_Callback(~, eventdata, handles)
% hObject    handle to File_SaveAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % check if there are any cells that have not been annotated - if so
    % prompt the user that he hasnt finished annotation
    flagCellPatternAnnotation = ([handles.data.cellStats.cellPatternId] > 1);
    
    if any(~flagCellPatternAnnotation)
        numCellsUnannotated = numel(find(~flagCellPatternAnnotation));
        strQuestion = sprintf( '%d/%d cells have not been annotated yet.\nDo you still want to save to resume annotation later?', numCellsUnannotated, numel(flagCellPatternAnnotation) );
        button = questdlg( strQuestion, 'Annotation Incomplete', 'Yes', 'No', 'No' );
        
        if strcmp(button, 'No')
           return; 
        end
    end
    
    % ask user to select the directory in which to save the data    
    if isfield( handles.history, 'lastOutputDir' )
        outputDir = uigetdir( handles.history.lastOutputDir, 'Select annotation output directory');
    else
        outputDir = uigetdir( handles.history.lastAnalyzedNucleusDir, 'Select annotation output directory');
    end

    if ~outputDir
        return;
    end

    handles.history.lastOutputDir = outputDir;
    
    % save data
    [pathstr, name, ext] = fileparts( handles.data.dataFilePath );
    outputDir = fullfile(outputDir, name);

    if ~isdir( outputDir )
        mkdir( outputDir );
    end
    
    data = handles.data;
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save( fullfile(outputDir, 'CellSegmentationQualityAnnotation.mat'), '-struct', 'data' );
    close(h);
    
    clear data;

    % output summary file
    summary_fid = fopen( fullfile(outputDir, 'annotationSummary.txt'), 'w' );    
        
    PrettyPrintStepDescription( 'Segmentation Quality Annotation Summary', summary_fid );
    
        % print some information about the dataset
        fprintf( summary_fid, '\n>> Dataset Description:\n' );
        
        fprintf( summary_fid, '\n\tHistone data file -- %s\n', handles.data.dataFilePath );
        fprintf( summary_fid, '\n\tImage Size - [ %s ]\n', sprintf( ' %d ', handles.data.metadata.imageSize) );
        fprintf( summary_fid, '\n\tImage Spacing - [ %s ]\n', sprintf( ' %.2f ', handles.data.metadata.pixelSize) );
        
        % print information about the annotation
        fprintf( summary_fid, '\n>> Annotation Summary:\n' );
        
        numTotalCells = numel( handles.data.cellStats );
        numCellsAnnotated = numel(find(flagCellPatternAnnotation));
        
        fprintf( summary_fid, '\n\t%d cells were found by the segmentation algorithm\n', numTotalCells );
        
        fprintf( summary_fid, '\n\t%d/%d (%.2f%%) cells have been annotated\n', numCellsAnnotated, numTotalCells, 100*numCellsAnnotated/numTotalCells );
        
        for i = 2:numel(handles.data.cellPatternTypes)          

            numCellsInCurClass = numel(find([handles.data.cellStats.cellPatternId] == i));
            fprintf( summary_fid, '\n\t-- %s -- %d (%.2f%%) cells\n', ...
                     handles.data.cellPatternTypes{i}, ...
                     numCellsInCurClass, 100*numCellsInCurClass/numTotalCells );

            if numCellsInCurClass > 2
                
                curClassCellStats = handles.data.cellStats([handles.data.cellStats.cellPatternId] == i);

                fprintf( summary_fid, '\n\t\tIntensity Statistics:\n');
                fprintf( summary_fid, '\n\t\tMin-Max: [%.2f, %.2f]\n', ...
                                      min( [curClassCellStats.meanIntensity] ), ...
                                      max( [curClassCellStats.meanIntensity] ) );
                fprintf( summary_fid, '\n\t\tMean-std: [%.2f, %.2f]\n', ...
                                      mean( [curClassCellStats.meanIntensity] ), ...
                                      std( [curClassCellStats.meanIntensity] ) );
                                      
                fprintf( summary_fid, '\n\t\tCell Volume (cu um):\n' );
                fprintf( summary_fid, '\n\t\t\tMin-Max: [%.2f, %.2f]\n', min( [curClassCellStats.AreaPhysp] ), max( [curClassCellStats.AreaPhysp] ) );
                fprintf( summary_fid, '\n\t\t\tMean-std: [%.2f, %.2f]\n', mean( [curClassCellStats.AreaPhysp] ), std( [curClassCellStats.AreaPhysp] ) );

            end
            
        end             

    fclose( summary_fid );
    
    % ask if the user wants to save annotated cell images
    strQuestion = 'Do you want save images of the annotated cell patterns?';
    button = questdlg( strQuestion, 'Save annotated cell images', 'Yes', 'No', 'Yes' );
    flagSaveImages = strcmp( button, 'Yes' ); 
    szOutputImage = [100, 100];
    
    h = waitbar(0, 'Saving Images of Annotated Cell Patterns ... Please Wait' );
    if flagSaveImages
        
        % create sub-directories for each cell class
        for i = 1:numel(handles.data.cellPatternTypes)
            if isdir( fullfile(outputDir, handles.data.cellPatternTypes{i}) )
                rmdir( fullfile(outputDir, handles.data.cellPatternTypes{i}), 's' );                   
            end
            mkdir( fullfile(outputDir, handles.data.cellPatternTypes{i}) );
        end

        for cellId = 1:numel(handles.data.cellStats)

            curCellStats = handles.data.cellStats(cellId);
            curCellCentroid = curCellStats.Centroid;
            curCellBoundingBox = curCellStats.BoundingBox;
            curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

            curCellPatternType = handles.data.cellStats(cellId).cellPatternType;
            curCellOutputDir = fullfile(outputDir, curCellPatternType); 

            subinds = cell(1,3);
            imsize = size(handles.data.imageData{handles.data.metadata.channelIdNuclei});
            for i = 1:2

                xi = round(curCellCentroid(3-i) - 0.5 * curCellDisplaySize);

                xi_low = xi;
                if xi_low < 1 
                    xi_low = 1;
                end

                xi_high = xi + curCellDisplaySize - 1;
                if xi_high > imsize(i)
                    xi_high = imsize(i);
                end

                subinds{i} = xi_low:xi_high;

            end    
            subinds{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));
            
            % MIP
            imCurCellCropped = handles.data.imageData{handles.data.metadata.channelIdNuclei}(subinds{1:2}, :);
            imCurCellSegCropped = (handles.data.imLabelCellSeg(subinds{1:2}, :) == cellId);
            
            imCurCellMIP = imresize( mat2gray(max(imCurCellCropped .* imCurCellSegCropped, [], 3)), szOutputImage);

            imwrite( imCurCellMIP, fullfile(curCellOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   

            % Mid slices
            imCurCellSegMidSliceBndCropped = imresize( bwperim( imCurCellSegCropped(:, :, round(curCellCentroid(3))) ), szOutputImage, 'nearest'); 
            
            imCurCellCroppedMidSlice = mat2gray( handles.data.imageData{handles.data.metadata.channelIdNuclei}( subinds{1:2}, round(curCellCentroid(3)) ), handles.dataDisplay.imDisplayRange(1,:) );
            imCurCellCroppedMidSlice = imresize(imCurCellCroppedMidSlice, szOutputImage);

            imwrite( genImageMaskOverlay(imCurCellCroppedMidSlice, imCurCellSegMidSliceBndCropped, [1, 0, 0], 0.5), ...
                     fullfile(curCellOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   

            waitbar( cellId/numel(handles.data.cellStats), h);
            
        end
        
    end
    
    closeStatusDialog(h);
    
    % Update handles structure
    guidata(gcbo, handles);

% --------------------------------------------------------------------
function [imMIP] = genMIPImageForDisplay( im, mask )

    if exist( 'mask', 'var' )

        imMIP_XY = max( im, [], 3);
        imMIP_XZ = (squeeze( max( im, [], 1) ))';
        imMIP_YZ = (squeeze( max( im, [], 2) ));

        imMaskMIP_XY = max( mask, [], 3);
        imMaskMIP_XZ = (squeeze( max( mask, [], 1) ))';
        imMaskMIP_YZ = (squeeze( max( mask, [], 2) ));
        
        gap = 5;
        imsize = size(imMIP_XY);
        zsize = size(im,3);
        
        imMIP = ones( imsize + gap + zsize );
        imMIP(1:imsize(1),1:imsize(2)) = imMIP_XY;
        imMIP( 1:imsize(1), imsize(2) + gap + (1:zsize) ) = imMIP_YZ;
        imMIP( imsize(1) + gap + (1:zsize), 1:imsize(2) ) = imMIP_XZ;  

        imMaskMIP = zeros( imsize + gap + zsize );
        imMaskMIP(1:imsize(1),1:imsize(2)) = imMaskMIP_XY;
        imMaskMIP( 1:imsize(1), imsize(2) + gap + (1:zsize) ) = imMaskMIP_YZ;
        imMaskMIP( imsize(1) + gap + (1:zsize), 1:imsize(2) ) = imMaskMIP_XZ;  
        
        imMIP = genImageMaskOverlay( imMIP, imMaskMIP, [0, 1, 0], 0.2 );
        
    else

        imMIP_XY = max( im, [], 3);
        imMIP_XZ = (squeeze( max( im, [], 1) ))';
        imMIP_YZ = (squeeze( max( im, [], 2) ));

        gap = 5;
        imsize = size(imMIP_XY);
        zsize = size(im,3);
        imMIP = ones( imsize + gap + zsize );
        imMIP(1:imsize(1),1:imsize(2)) = imMIP_XY;
        imMIP( 1:imsize(1), imsize(2) + gap + (1:zsize) ) = imMIP_YZ;
        imMIP( imsize(1) + gap + (1:zsize), 1:imsize(2) ) = imMIP_XZ;  
        
    end    
    
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


% --- Executes when user attempts to file_close CellSegmentationQualityAnnotator.
function CellSegmentationQualityAnnotator_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CellSegmentationQualityAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
    
    button = questdlg('Are you sure you want to close?', 'Closing Annotation Tool', 'Yes', 'No', 'No');
    if strcmp( button, 'No' )
        return;
    end
    
    % file_close
    delete(hObject);
    
function PrettyPrintStepDescription( strStepDescription, fid )

    if ~exist( 'fid', 'var' )        
        fid = 1; % print to command-line
    end
    
    strStar = strStepDescription;
    strStar(:) = '*';
    strStar = [ '****', strStar, '****' ];
    fprintf( fid, '\n\n%s', strStar );
    fprintf( fid, '\n    %s    ', strStepDescription );
    fprintf( fid, '\n%s\n\n', strStar );

% --- Executes on button press in ShowPreviousCell.
function ShowPreviousCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowPreviousCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    % decrement cell id
    handles.dataDisplay.curCellId = max(1, handles.dataDisplay.curCellId - 1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in ShowNextCell.
function ShowNextCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNextCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % increment cell id
    handles.dataDisplay.curCellId = min( numel(handles.data.cellStats), handles.dataDisplay.curCellId + 1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in ShowNextUnannotatedCell.
function ShowNextUnannotatedCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNextUnannotatedCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded || isempty(handles.data.unannotatedCellList) 
        return;
    end

    % set cell id to next unannotated cell
    handles.dataDisplay.curCellId = handles.data.unannotatedCellList(1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in ShowLastCell.
function ShowLastCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowLastCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % set cell id as last one
    handles.dataDisplay.curCellId = numel( handles.data.cellStats );
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in ShowFirstCell.
function ShowFirstCell_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFirstCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % decrement cell id
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in CheckboxGlobalSegMask.
function CheckboxGlobalSegMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxGlobalSegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxGlobalSegMask

    % change mask use mode
    handles.flagShowGlobalSegMask = get(hObject,'Value');
    
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
    handles.flagUseLOG = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    

% --------------------------------------------------------------------
function File_Set_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_Set_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters.nucleiSegmentation = NucleiSegmentationParametersGUI( handles.parameters.nucleiSegmentation );

    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object deletion, before destroying properties.
function CellSegmentationQualityAnnotator_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellSegmentationQualityAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    if isfield(handles, 'history')
        history = handles.history;
        [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
        historyFile = fullfile( pathstr, 'CellSegmentationQualityAnnotatorHistory.mat' );
        save( historyFile, '-struct', 'history' );  
    end

    % close matlab pool
    if ~handles.parameters.flagPoolOpenedAlready && matlabpool( 'size' ) > 0
        matlabpool close;
    end

% --- Executes on button press in CheckboxCellBBox.
function CheckboxCellBBox_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxCellBBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxCellBBox

    % change mask use mode
    handles.flagShowCellBBox = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --------------------------------------------------------------------
function FnSliceScroll_Callback(hSrc, eventdata)

    handles = guidata(hSrc);
    
    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{handles.data.metadata.channelIdNuclei});

    if IsMouseInsideAxes(handles.Axes_Global_XZ)
        viewId = 1; % y-slice
    elseif IsMouseInsideAxes(handles.Axes_Global_YZ)
        viewId = 2; % x-slize
    else
        viewId = 3; % z-slice
    end
    
    if eventdata.VerticalScrollCount > 0
        if handles.dataDisplay.curCellSliceId(viewId) < imsize(viewId); 
            handles.dataDisplay.curCellSliceId(viewId) = handles.dataDisplay.curCellSliceId(viewId) + 1;
        end
    elseif eventdata.VerticalScrollCount < 0
        if handles.dataDisplay.curCellSliceId > 1
            handles.dataDisplay.curCellSliceId(viewId) = handles.dataDisplay.curCellSliceId(viewId) - 1;
        end
    end
    
    guidata(hSrc, handles);    
    UpdateCellDisplay(handles);

% --- Executes on button press in CheckboxShowNucleusChannel.
function CheckboxShowNucleusChannel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxShowNucleusChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxShowNucleusChannel

    % change log use mode
    handles.data.flagShowNucleusChannel = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes during object creation, after setting all properties.
function EditboxSliceId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditboxSliceId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function View_Cell_Segmentation_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to View_Cell_Segmentation_In_Imaris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % make sure data is loaded
    if ~handles.flagDataLoaded
        return;
    end
    
    % make sure we are in windows
    if ~ispc 
        errordlg( 'Sorry. This feature is only available in windows' );
        return;
    end

    % make sure imaris is installed
    try
        t = actxserver('Imaris.Application');
        clear t;
    catch err
        errordlg( 'Unable to connect to imaris. To use this feature, make sure imaris is installed' );
        return;
    end
    
    curCellId = handles.dataDisplay.curCellId;
    curCellStats = handles.data.cellStats( curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

    % create crop indices
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{handles.data.metadata.channelIdNuclei});
    for i = 1:2
        
        xi = round(curCellCentroid(3-i) - 0.5 * curCellDisplaySize);
        
        xi_low = xi;
        if xi_low < 1 
            xi_low = 1;
        end
        
        xi_high = xi + curCellDisplaySize - 1;
        if xi_high > imsize(i)
            xi_high = imsize(i);
        end
        
        subinds{i} = xi_low:xi_high;
        
    end    
    subinds{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));    
    
    % generate visualization
    imCellCropped = cell( size(handles.data.imageData) );
    for i = 1:numel(handles.data.imageData)
        imCellCropped{i} = handles.data.imageData{i}(subinds{:});
    end

    imvis = ImarisDataVisualizer(cat(4, imCellCropped{:}), ...
                                 'spacing', handles.data.metadata.pixelSize);
    
    handles.imarisAppCellSegCropped = imvis;
    hSegmentation = imvis.AddDataContainer();

        % compute isosurface geometry for cells in each pattern
        imCurCellMask = handles.data.imLabelCellSeg(subinds{:}) == handles.dataDisplay.curCellId;
        surfaceQuality = 1.0;
        
        curCellGeometry = ImarisDataVisualizer.generateSurfaceFromMask(imCurCellMask, 'surfaceQuality', surfaceQuality);

        imvis.AddSurfaces(curCellGeometry, hSegmentation, ...
                          'name', sprintf( 'CellSeg_%d', handles.dataDisplay.curCellId), ...
                          'color', handles.data.CellSegColorMap(handles.dataDisplay.curCellId, :) );

        % display cell seed points
        imCellSeedCropped = handles.dataDisplay.imCellSeedPoints(subinds{:});
        imCellSeedCropped(~imCurCellMask) = 0;
        stats = regionprops( bwlabeln(imCellSeedCropped), 'Centroid' );
        cellSeedPointLocations = cat( 1, stats.Centroid );
    
        imvis.AddSpots(cellSeedPointLocations, 0, ...
                       'hContainer', hSegmentation, ...
                       'name', 'Nuclei Seed Points', 'color', [1, 0, 0]);
    
                      
    % Update handles structure
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function View_Full_Segmentation_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to View_Full_Segmentation_In_Imaris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % make sure data is loaded
    if ~handles.flagDataLoaded
        return;
    end
    
    % make sure we are in windows
    if ~ispc 
        errordlg( 'Sorry. This feature is only available in windows' );
        return;
    end

    % make sure imaris is installed
    try
        t = actxserver('Imaris.Application');
        clear t;
    catch err
        errordlg( 'Unable to connect to imaris. To use this feature, make sure imaris is installed' );
        return;
    end

    % generate visualization
    imvis = ImarisDataVisualizer( cat(4, handles.data.imageData{:}), 'spacing', handles.data.metadata.pixelSize );
    handles.imarisAppCellSeg = imvis;
    
    hSegmentation = imvis.AddDataContainer();

        % display cell seed point locations
        stats = regionprops( bwlabeln( handles.dataDisplay.imCellSeedPoints ), 'Centroid' );
        cellSeedPointLocations = cat( 1, stats.Centroid );
    
        imvis.AddSpots(cellSeedPointLocations, 0, ...
                       'hContainer', hSegmentation, ...
                       'name', 'Nuclei Seed Points', 'color', [1, 0, 0]);

    
        % compute isosurface geometry for cells in each pattern
        hStatusDlg = waitbar( 0, 'computing surface geometry of cells' );
        surfaceQuality = 1.0;
        for cid = 1:numel(handles.data.cellStats)

            imCurCellMask = handles.data.imLabelCellSeg == cid;
            curCellGeometry = ImarisDataVisualizer.generateSurfaceFromMask(imCurCellMask, 'surfaceQuality', surfaceQuality);
            
            imvis.AddSurfaces(curCellGeometry, hSegmentation, ...
                              'name', sprintf( 'CellSeg_%d', cid), ...
                              'color', handles.data.CellSegColorMap(cid, :) );
            
            waitbar( cid/numel( handles.data.cellStats ), hStatusDlg );

        end
        closeStatusDialog( hStatusDlg );
                   
    % Update handles structure
    guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function CellDescriptorsTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellDescriptorsTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function LabelCellDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelCellDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function CellSegmentationQualityAnnotator_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to CellSegmentationQualityAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{handles.data.metadata.channelIdNuclei});
    newCellSliceId = handles.dataDisplay.curCellSliceId;
    
    if IsMouseInsideAxes(handles.Axes_Global_XZ)
        
        viewId = 1; % y-slice
        relMousePos = getRelativeMousePos(handles.Axes_Global_XZ);
        pixelPos = floor(1 + imsize([2,3]) .* relMousePos);
        newCellSliceId([2, 3]) = pixelPos;
        
    elseif IsMouseInsideAxes(handles.Axes_Global_YZ)
        
        viewId = 2; % x-slize
        relMousePos = getRelativeMousePos(handles.Axes_Global_YZ);
        pixelPos = fliplr(floor( 1 + imsize([3,1]) .* relMousePos ));
        newCellSliceId([1, 3]) = pixelPos;
        
    elseif IsMouseInsideAxes(handles.Axes_Global_XY)
        
        viewId = 3; % z-slice
        relMousePos = getRelativeMousePos(handles.Axes_Global_XY);
        pixelPos = fliplr(floor( 1 + imsize([2,1]) .* relMousePos ));
        newCellSliceId([1, 2]) = pixelPos;
        
    else
        return;
    end

    % get cell id under mouse click
    indViewSlicer = { ':', ':', ':' };
    indViewSlicer{viewId} = handles.dataDisplay.curCellSliceId(viewId);
    imCurLabelCellSeg = squeeze( handles.data.imLabelCellSeg( indViewSlicer{:} ) );

    pixelPos = num2cell( pixelPos );
    curCellId = imCurLabelCellSeg( pixelPos{:} );
    handles.dataDisplay.curCellSliceId = newCellSliceId;
    
    if curCellId >= 1
        
        % set cell id
        handles.dataDisplay.curCellId = imCurLabelCellSeg( pixelPos{:} );

        % update cell pattern listbox
        set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    end
        
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);    


% --- Executes on button press in btnShowPreviousCellInSelectedClass.
function btnShowPreviousCellInSelectedClass_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowPreviousCellInSelectedClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded
        return;
    end

    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indPrevCellInCurClass = find( cellPatternVec(handles.dataDisplay.curCellId-1:-1:1) == get(handles.poplistCellClassSelector,'Value') );
    
    if isempty( indPrevCellInCurClass )
        return;
    end
    
    indPrevCellInCurClass = handles.dataDisplay.curCellId - indPrevCellInCurClass;
    
    % decrement cell id
    handles.dataDisplay.curCellId = indPrevCellInCurClass(1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in btnShowNextCellInSelectedClass.
function btnShowNextCellInSelectedClass_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowNextCellInSelectedClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded
        return;
    end

    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indNextCellInCurClass = find( cellPatternVec(handles.dataDisplay.curCellId+1:end) == get(handles.poplistCellClassSelector,'Value') );
        
    if isempty( indNextCellInCurClass )
        return;
    end
    
    indNextCellInCurClass = handles.dataDisplay.curCellId + indNextCellInCurClass;
    
    % decrement cell id
    handles.dataDisplay.curCellId = indNextCellInCurClass(1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on selection change in poplistCellClassSelector.
function poplistCellClassSelector_Callback(hObject, eventdata, handles)
% hObject    handle to poplistCellClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistCellClassSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistCellClassSelector

    % first check if data has been loaded
    if ~handles.flagDataLoaded
        return;
    end
    
    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indFirstCellInCurClass = find( cellPatternVec == get(hObject,'Value') );
    
    if isempty( indFirstCellInCurClass )
        return;
    end
    
    % decrement cell id
    handles.dataDisplay.curCellId = indFirstCellInCurClass(1);
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % update cell pattern listbox
    set(handles.ListboxSegmentationQualitySelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function poplistCellClassSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistCellClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
