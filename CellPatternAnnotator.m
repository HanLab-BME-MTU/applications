function varargout = CellPatternAnnotator(varargin)
%CELLPATTERNANNOTATOR M-file for CellPatternAnnotator.fig
%      CELLPATTERNANNOTATOR, by itself, creates a new CELLPATTERNANNOTATOR or raises the existing
%      singleton*.
%
%      H = CELLPATTERNANNOTATOR returns the handle to a new CELLPATTERNANNOTATOR or the handle to
%      the existing singleton*.
%
%      CELLPATTERNANNOTATOR('Property','Value',...) creates a new CELLPATTERNANNOTATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CellPatternAnnotator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CELLPATTERNANNOTATOR('CALLBACK') and CELLPATTERNANNOTATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CELLPATTERNANNOTATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellPatternAnnotator

% Last Modified by GUIDE v2.5 15-Apr-2014 12:14:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellPatternAnnotator_OpeningFcn, ...
                   'gui_OutputFcn',  @CellPatternAnnotator_OutputFcn, ...
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

% --- Executes just before CellPatternAnnotator is made visible.
function CellPatternAnnotator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

    p = inputParser;
    p.addParamValue('flagParallelize', false, @isscalar);
    p.addParamValue('flagDebugMode', false, @isscalar);
    p.parse( varargin{:} );    
    PARAMETERS = p.Results;

    % Choose default command line output for CellPatternAnnotator
    handles.output = hObject;

    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
    if exist( historyFile, 'file' )
       handles.history = load( historyFile );
    else
       handles.history.lastAnalyzedNucleusDir = pathstr;
       handles.history.lastAnalyzedRedGreenDir = pathstr;
    end

    % Initialize global variables
    handles.cellDisplaySize = 70;     
    
    handles.flagDataLoaded = false;
    handles.flagReviewMode = false;
    
    handles.flagUseLOG = get(handles.CheckboxGlobalLog, 'Value');
    handles.flagShowGlobalSegMask = get(handles.CheckboxGlobalSegMask, 'Value');
    handles.flagShowCellBBox = get(handles.CheckboxCellBBox, 'Value');
    handles.flagShowInvalidROIMask = get(handles.CheckboxValidROIMask, 'Value');
    
    handles.defaultCellPatternTypes = cellstr( get(handles.ListboxCellPatternSelector, 'String') );
    
    handles.imarisApp = [];
    handles.imarisAppCellCropped = [];
    handles.imarisAppCellSegCropped = [];
    
    % default parameters
    f = rdir( fullfile(fileparts(mfilename('fullpath')), '**', 'regionMerging.model') );
    if ~isempty(f)
        handles.defaultParameters.segmentation = NucleiSegmentationParametersGUI( 'defaultWithRegionMerging' );
        handles.defaultParameters.segmentation.flagPerformRegionMerging = true;
        handles.defaultParameters.segmentation.regionMergingModelFile = f(1).name;
    else
        handles.defaultParameters.segmentation = NucleiSegmentationParametersGUI( 'defaultWithoutRegionMerging' );
    end
    
    handles.defaultParameters.flagParallelize = PARAMETERS.flagParallelize;
    handles.defaultParameters.flagDebugMode = PARAMETERS.flagDebugMode;
    
    handles.parameters = handles.defaultParameters;
    
    % open matlab pool for parallel processing    
    handles.parameters.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.parameters.flagParallelize && ~handles.parameters.flagPoolOpenedAlready 
        matlabpool open;
    end
    
    % make sure weka is in the path
    AddWekaClassesToPath();
    
    % Set callbacks
    set(handles.CellPatternAnnotator, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Initialize the state of some UI controls
    set( handles.LabelAllChannels, 'String', 'All Channel Overlay' );
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes CellPatternAnnotator wait for user response (see UIRESUME)
    % uiwait(handles.CellPatternAnnotator);

% --- Outputs from this function are returned to the command line.
function varargout = CellPatternAnnotator_OutputFcn(hObject, eventdata, handles)
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

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

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

    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    
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
    [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedNucleusDir, '*.oif; *.oib' ), ...
                                               'Select the data file containing Histone-2B CFP' );   
    
    if ~fileName 
        return;
    end
    
    dataFilePath{1} = fullfile( pathName, fileName );
    handles.history.lastAnalyzedNucleusDir = pathName;
    handles.history.lastAnalyzedRedGreenDir = pathName;
    
    % ask the use to select the oif file for the red/green channel
    [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedRedGreenDir, '*.oif; *.oib' ), ...
                                               'Select the data file containing the FUCCI Cell Cycle Reporter' );   
        
    if ~fileName 
        return;
    end
    
    dataFilePath{2} = fullfile( pathName, fileName );        
    handles.history.lastAnalyzedRedGreenDir = pathName;
    
    % load nucleus channel data
    PrettyPrintStepDescription( 'Loading Histone Channel Data' );
    
    hStatusDialog = waitbar(0, 'Loading Histone Channel Data');
    try        
        imageSeriesNucleus = loadIntravitalDataset( dataFilePath{1}  );    
        if (imageSeriesNucleus(1).metadata.voxelSpacing(1)/ imageSeriesNucleus(1).metadata.voxelSpacing(3)) >= 1
            imageSeriesNucleus(1).metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
        end        
        imageSeriesNucleus(1).metadata.channelColors = [ 0 , 0 , 1 ];
    catch err
        err
        err.message
        fprintf( 'ERROR: could not load nucleus channel data from file %s', dataFilePath{1} );
        errordlg( sprintf( 'Could not load nucleus channel data from file %s\n\n%s', dataFilePath{1}, err.message) ); 
        closeStatusDialog( hStatusDialog );
        return;
    end
    
    metadata_NucleusChannel = imageSeriesNucleus(1).metadata
    
    % load red/green channel data
    PrettyPrintStepDescription( 'Loading Red/Green Channel Data' );

    waitbar(0, hStatusDialog, 'Loading Red/Green Channel Data');
    try

        imageSeriesFUCCI = loadIntravitalDataset( dataFilePath{2}  );    

        if imageSeriesFUCCI(1).metadata.channelExcitationWavelength(1) > imageSeriesFUCCI(1).metadata.channelExcitationWavelength(2)

            % first channel is not green so swap it 
            warning( 'GFP was not the first channel in confocal data. It will be swapped to the first place');
            imageSeriesFUCCI(1).imageData = imageSeriesFUCCI(1).imageData(:,[2,1]);
            imageSeriesFUCCI(1).metadata.channelExcitationWavelength = imageSeriesFUCCI(1).metadata.channelExcitationWavelength([2,1]);

        end

        imageSeriesFUCCI(1).metadata.channelColors = [ 0 1 0; 1 0 0 ];

        metadata_FUCCI = imageSeriesFUCCI(1).metadata

    catch err

        err
        err.message
        fprintf( 'ERROR: could not load red-green channel data from file %s', dataFilePath{2} );
        errordlg( sprintf( 'Could not load red-green channel data from file %s\n\n%s', dataFilePath{2}, err.message ) ); 
        closeStatusDialog( hStatusDialog );
        return;
    end

    if any( imageSeriesFUCCI(1).metadata.numChannels ~= 2  )
        errordlg( sprintf( 'Red-Green channel data is expected to contain 2 channels. The file selected by you contains %d channels', imageSeriesFUCCI(1).metadata.numChannels ) );
        return;
    end

    if any( imageSeriesFUCCI(1).metadata.volSize ~= imageSeriesNucleus(1).metadata.volSize )
        errordlg( sprintf( 'Volume Size of Nucleus data doesnt match with red-green channel data' ) );
        closeStatusDialog( hStatusDialog );
        return;
    end

    if (imageSeriesFUCCI(1).metadata.voxelSpacing(1)/ imageSeriesFUCCI(1).metadata.voxelSpacing(3)) >= 1
        imageSeriesFUCCI(1).metadata.voxelSpacing = imageSeriesNucleus(1).metadata.voxelSpacing; % incorrect spacing in metadata - use something meaningful
    end

    if any( imageSeriesFUCCI(1).metadata.voxelSpacing ~= imageSeriesNucleus(1).metadata.voxelSpacing )
        
        strError = sprintf( '\nPixel spacing of nuclear marker data doesnt match with FUCCI data.\n\n\tSpacing of Histone data - [ %s ]\n\tSpacing of FUCCI data  - [ %s ]', ...
                               sprintf( ' %f ', imageSeriesNucleus(1).metadata.voxelSpacing ), ...
                               sprintf( ' %f ', imageSeriesFUCCI(1).metadata.voxelSpacing ) );
                           
        errordlg( strError, 'Incompatible voxel size accross channels');
        closeStatusDialog( hStatusDialog );
        return;
    end
    
    % check if the user wants to process only a substack    
    PrettyPrintStepDescription( 'Selected Range of Slices:' );
    
    imseriesshow( imageSeriesNucleus(1).imageData{1,1} );
    set(gcf, 'Name', 'Histone data before slice cropping');
    
    numSlices = metadata_NucleusChannel.volSize(3);
    strSliceRange = sprintf( '(1 - %d)', numSlices );
    dlgOptions.Resize = 'on';
    dlgOptions.WindowStyle = 'normal';
    selectedSliceRange = inputdlg( { sprintf('Start Slice %s:', strSliceRange), sprintf('End Slice %s:', strSliceRange)}, ...
                              'Select slice range of interest', 1, ...
                              { '1', num2str(numSlices) }, dlgOptions );
                          
    try
        
      selectedSliceRange = [ str2num(selectedSliceRange{1}), str2num(selectedSliceRange{2}) ];    
      assert( ~any( (selectedSliceRange - floor(selectedSliceRange)) > 0 ) ); % must be integers
      assert( all( selectedSliceRange >= 1 ) && all( selectedSliceRange <= numSlices ) );
      assert( selectedSliceRange(1) < selectedSliceRange(2) ); % start slice should be less than end slice
      
    catch        
        
        errordlg( sprintf( 'Invalid start/end slice. Must be integers between 1 - %d', numSlices ) );
        closeStatusDialog( hStatusDialog );
        return;
        
    end
    
    selectedSliceRange
    
    imageSeriesNucleus(1).imageData{1} = imageSeriesNucleus(1).imageData{1}(:, :, selectedSliceRange(1):selectedSliceRange(2));
    for i = 1:2
        imageSeriesFUCCI(1).imageData{1,i} = imageSeriesFUCCI(1).imageData{1,i}(:, :, selectedSliceRange(1):selectedSliceRange(2));
    end
    
    imageSeriesNucleus(1).metadata.volSize = size( imageSeriesNucleus(1).imageData{1} );
    
    % store image data in handles structures    
    handles.flagDataLoaded = true;
    handles.data = [];
    handles.data.cellPatternTypes = handles.defaultCellPatternTypes;
    handles.data.dataFilePath = dataFilePath;
    handles.data.metadata = imageSeriesNucleus(1).metadata;
    
    set( handles.CellPatternAnnotator, 'Name', sprintf( 'Cell Pattern Annotator - %s', handles.data.dataFilePath{1} ) );
    
    handles.data.metadata.channelExcitationWavelength = [ imageSeriesNucleus(1).metadata.channelExcitationWavelength, ...
                                                          imageSeriesFUCCI(1).metadata.channelExcitationWavelength ];

    handles.data.metadata.channelColors = cat( 1 , imageSeriesNucleus(1).metadata.channelColors, ...
                                                   imageSeriesFUCCI(1).metadata.channelColors );             

    metadata = handles.data.metadata
    
    % correct stage-shift        
    PrettyPrintStepDescription( 'Correcting shift between confocal and two-photon data' );

    waitbar(0, hStatusDialog, 'Correcting shift between confocal and two-photon data');
    [handles.data.imageData, ...
     handles.data.imRegValidMask] = CorrectTwoPhotonConfocalStageShift( imageSeriesNucleus(1).imageData, imageSeriesNucleus(1).metadata.voxelSpacing, ...
                                                                        imageSeriesFUCCI(1).imageData, imageSeriesFUCCI(1).metadata.voxelSpacing, ...
                                                                        'flagParallelize', handles.parameters.flagParallelize, ...
                                                                        'flagDebugMode', false);
                                                                    
    % pre-compute data needed for display
    handles = ComputeDisplayData( handles );
                                                                    
    % close progress bar
    closeStatusDialog( hStatusDialog );

    % Update handles structure
    guidata(hObject, handles);

    % Run analysis
    RunAnalysis(hObject, handles);
    
% --------------------------------------------------------------------    
function [ handles ] = ComputeDisplayData(handles)

    for i = 1:3
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
    if handles.parameters.flagPerformRegionMerging
        regionMergingModelFile = handles.parameters.segmentation.regionMergingModelFile;
    end
    
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints, ...
     handles.data.segAlgoParameters ] = segmentCellsInIntravitalData( handles.data.imageData{1}, ...
                                                                      handles.data.metadata.voxelSpacing, ...                                                                      
                                                                      'flagParallelize', handles.parameters.flagParallelize, ...
                                                                      'flagDebugMode', handles.parameters.flagDebugMode, ...
                                                                      'cellDiameterRange', handles.parameters.segmentation.cellDiameterRange, ...
                                                                      'thresholdingAlgorithm', 'MinErrorPoissonSliceBySliceLocal', ...
                                                                      'seedPointDetectionAlgorithm', handles.parameters.segmentation.seedPointDetectionAlgorithm, ...
                                                                      'minCellVolume', handles.parameters.segmentation.minCellVolume, ...
                                                                      'flagIgnoreCellsOnXYBorder', handles.parameters.segmentation.flagIgnoreXYBorderCells, ...
                                                                      'roiMask', handles.data.imRegValidMask, ...
                                                                      'regionMergingModelFile', regionMergingModelFile);

    [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.dataDisplay.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints > 0, ones(3,3,3) );
    closeStatusDialog(hStatusDialog);
    
    % compute cell properties
    handles.data.cellStats = ComputeCellProperties( handles );
    
    % initialize annotation
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
    
    % switch off review mode if segmentation is done again
    handles.flagReviewMode = false;
    
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
        cellPixelIntensities = handles.data.imageData{1}( cellStats(i).PixelIdxList );
        cellStats(i).meanIntensity = mean( cellPixelIntensities );
        cellStats(i).stdIntensity = std( double(cellPixelIntensities) );
        cellStats(i).minIntensity = min( cellPixelIntensities );
        cellStats(i).maxIntensity = max( cellPixelIntensities );
        
        % volume
        cellStats(i).AreaPhysp = cellStats(i).Area * prod( handles.data.metadata.voxelSpacing );
        
        % fit ellipsoid and store its radii
        [yind, xind, zind] = ind2sub( size(handles.data.imageData{1}), cellStats(i).PixelIdxList ); 
        ptCell = [xind, yind, zind] .* repmat( handles.data.metadata.voxelSpacing, [numel(xind), 1] );
        ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
        [U, S, V] = svd( (ptCell' * ptCell) / size(ptCell,1) );
        
        cellStats(i).ellipsoidRadiusPhysp = zeros(1,3); 
        for j = 1:3
            cellStats(i).ellipsoidRadiusPhysp(j) = 2 * sqrt(S(j,j)); % eigen-values are a measure of variance
        end
        
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
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );
    
    % display global cross section images  
    if handles.flagUseLOG               
        imGlobalXY = handles.dataDisplay.imageDataLOG{1}( :, :, curCellSliceId(3) );       
        imGlobalXZ = squeeze( handles.dataDisplay.imageDataLOG{1}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.dataDisplay.imageDataLOG{1}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imLogDisplayRange;
    else        
        imGlobalXY = handles.data.imageData{1}( :, :, curCellSliceId(3) );
        imGlobalXZ = squeeze( handles.data.imageData{1}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.data.imageData{1}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imDisplayRange;
    end
    imGlobalXY = mat2gray(imGlobalXY, displayrange(1,:) );
    imGlobalXZ = mat2gray(imGlobalXZ', displayrange(1,:) );
    imGlobalYZ = mat2gray(imGlobalYZ, displayrange(1,:) );
    
    overlayMasksGlobalXY = {};
    overlayMaskAlphaGlobalXY = [];
    
    overlayMasksGlobalXZ = {};
    overlayMaskAlphaGlobalXZ = [];
    
    overlayMasksGlobalYZ = {};
    overlayMaskAlphaGlobalYZ = [];
    
    if handles.flagShowGlobalSegMask 
        
        segMaskAlpha = 0.2;
        
        imGlobalXYSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask( :, :, curCellSliceId(3), : ) );
        overlayMasksGlobalXY{end+1} = imGlobalXYSegMaskRGB;
        overlayMaskAlphaGlobalXY(end+1) = segMaskAlpha;
        
        imGlobalXZSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask(curCellSliceId(1), :, :, :) );
        imGlobalXZSegMaskRGB = permute(imGlobalXZSegMaskRGB, [2,1,3] );
        overlayMasksGlobalXZ{end+1} = imGlobalXZSegMaskRGB;
        overlayMaskAlphaGlobalXZ(end+1) = segMaskAlpha;
        
        imGlobalYZSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask(:, curCellSliceId(2), :, :) );
        overlayMasksGlobalYZ{end+1} = imGlobalYZSegMaskRGB;
        overlayMaskAlphaGlobalYZ(end+1) = segMaskAlpha;
        
    end
    
    if handles.flagShowInvalidROIMask

        roiMaskAlpha = 0.5;
        
        imGlobalXYInvalidROIMask = 1 - repmat( squeeze( handles.data.imRegValidMask( :, :, curCellSliceId(3), : ) ), [1, 1, 3] );
        imGlobalXYInvalidROIMask(:,:,2) = 0;
        overlayMasksGlobalXY{end+1} = imGlobalXYInvalidROIMask;
        overlayMaskAlphaGlobalXY(end+1) = roiMaskAlpha;
        
        imGlobalXZInvalidROIMask = 1 - squeeze( handles.data.imRegValidMask(curCellSliceId(1), :, :, :) );
        imGlobalXZInvalidROIMask = repmat( permute(imGlobalXZInvalidROIMask, [2,1,3] ), [1, 1, 3] );
        imGlobalXZInvalidROIMask(:,:,2) = 0;
        overlayMasksGlobalXZ{end+1} = imGlobalXZInvalidROIMask;
        overlayMaskAlphaGlobalXZ(end+1) = roiMaskAlpha;
        
        imGlobalYZInvalidROIMask = 1 - repmat( squeeze( handles.data.imRegValidMask(:, curCellSliceId(2), :, :) ), [1, 1, 3] );
        imGlobalYZInvalidROIMask(:,:,2) = 0;
        overlayMasksGlobalYZ{end+1} = imGlobalYZInvalidROIMask;
        overlayMaskAlphaGlobalYZ(end+1) = roiMaskAlpha;
        
    end
    
    if ~isempty(overlayMasksGlobalXY)
        
        imGlobalXYDisplay = genImageRGBMaskOverlay( imGlobalXY, overlayMasksGlobalXY, overlayMaskAlphaGlobalXY );
        imGlobalXZDisplay = genImageRGBMaskOverlay( imGlobalXZ, overlayMasksGlobalXZ, overlayMaskAlphaGlobalXZ );
        imGlobalYZDisplay = genImageRGBMaskOverlay( imGlobalYZ, overlayMasksGlobalYZ, overlayMaskAlphaGlobalYZ );
        
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
    
    set( handles.Axes_Global_XZ, 'XTickLabel', [], 'YTickLabel', []);
    set( handles.Axes_Global_XZ, 'XLim', [1 512], 'YLim', [1 size(imGlobalXZDisplay,1)]);
    
    set( handles.Axes_Global_YZ, 'XTickLabel', [], 'YTickLabel', []);        
    set( handles.Axes_Global_YZ, 'YLim', [1 512], 'XLim', [1 size(imGlobalYZDisplay,2)]);
    
    % draw bounding box around each cell
    if handles.flagShowCellBBox
        
        imsize = size( handles.data.imageData{1} );
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
    imsize = size(handles.data.imageData{1});
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
        imCellCropped = mat2gray( handles.data.imageData{1}(subinds{:}), handles.dataDisplay.imDisplayRange(1,:) );        
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
    imCurCellMIP = mat2gray( max( double( handles.data.imageData{1}( subindsMIP{:} ) ), [], 3 ) .* imCurCellSegMIP );
    imHistoneMIPDisplay = repmat( imCurCellMIP, [1,1,3] );
    
    cla( handles.Axes_Histone_MIP, 'reset' );
    image( imHistoneMIPDisplay, 'Parent', handles.Axes_Histone_MIP );
    set( handles.Axes_Histone_MIP, 'XTickLabel', [], 'YTickLabel', [] );
    
    % display FUCCI stuff
    imCellChannelData = [];
    displaycolors = [];
    for i = 1:3
       imCurChannelData = mat2gray(handles.data.imageData{i}(subinds{:}), handles.dataDisplay.imDisplayRange(i,:));
       imCellChannelData = cat(3, imCellChannelData, imCurChannelData );
       displaycolors = [ displaycolors; handles.data.metadata.channelColors(i,:) ];
    end
    
        % display FUCCI RFP
        imFUCCI_RFP = genMultiChannelOverlay(imCellChannelData(:,:,3), displaycolors(3,:));            
        
        cla( handles.Axes_FUCCI_RFP, 'reset' );
        image( imFUCCI_RFP, 'Parent', handles.Axes_FUCCI_RFP );
        set( handles.Axes_FUCCI_RFP, 'XTickLabel', [], 'YTickLabel', [] );
        
        % display FUCCI GFP
        imFUCCI_GFP = genMultiChannelOverlay(imCellChannelData(:,:,2), displaycolors(2,:));            
        
        cla( handles.Axes_FUCCI_GFP, 'reset' );
        image( imFUCCI_GFP, 'Parent', handles.Axes_FUCCI_GFP );
        set( handles.Axes_FUCCI_GFP, 'XTickLabel', [], 'YTickLabel', [] );        
        
        % display FUCCI overlay
        imFUCCIOverlay = genMultiChannelOverlay(imCellChannelData(:,:,2:3), displaycolors(2:3,:));            
        
        cla( handles.Axes_FUCCI_Overlay, 'reset' );
        image( imFUCCIOverlay, 'Parent', handles.Axes_FUCCI_Overlay );
        set( handles.Axes_FUCCI_Overlay, 'XTickLabel', [], 'YTickLabel', [] );                
    
        % display All channel overlay
        imAllChannelOverlay = genMultiChannelOverlay(imCellChannelData, displaycolors);            
        
        cla( handles.Axes_AllChannel_XY, 'reset' );
        image( imAllChannelOverlay, 'Parent', handles.Axes_AllChannel_XY );
        set( handles.Axes_AllChannel_XY, 'XTickLabel', [], 'YTickLabel', [] );                
        
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
function UpdateRedGreenDifferenceDisplay(handles)
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );   
    
    prev_axis = gca;
    
    % show histograms of individual channels    
    axes( handles.Axes_Channel_Histogram );
    cla reset;    
    hold on;
    
    maxval = 4096;
    for chid = 1:3
        
        channelPixelIntensties{chid} = maxval * mat2gray(handles.data.imageData{chid}( curCellStats.PixelIdxList), ...
                                                       handles.dataDisplay.imDisplayRange(chid,:));
                                                   
    end        
    
    warning( 'off', 'all' );
    nhist( channelPixelIntensties, ...
           'nolegend', 'noerror', 'pdf', 'smooth', 'proportion', ...
           'minx', 0, 'maxx', maxval, ...
           'maxbins', 16, ...
           'color', handles.data.metadata.channelColors );        
    warning( 'on', 'all' );
    
    hold off;
    
    legend( {'CFP', 'GFP', 'RFP'}, 'Location', 'best', 'FontSize', 8.0 );
    title( 'Channel Histograms of Cell Pixels', 'FontSize', 11.0 ); 
    
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
    strCellDescription = sprintf( '%s\n\nFitted Ellipsoid Radii (um): [%.2f, %.2f, %.2f]', ...
                                  strCellDescription, ...
                                  curCellStats.ellipsoidRadiusPhysp(1), ...
                                  curCellStats.ellipsoidRadiusPhysp(2), ...
                                  curCellStats.ellipsoidRadiusPhysp(3) );
    strCellDescription = sprintf( '%s\n\nFitted Ellipsoid Volume (cu um): %.2f', ...
                                  strCellDescription, ...
                                  (4 * pi * prod(curCellStats.ellipsoidRadiusPhysp)/ 3) );
                  
    % display difference between pixel intensities of red ang green channels 
    curGreenPixelIntensities = mat2gray(handles.data.imageData{2}( curCellStats.PixelIdxList), handles.dataDisplay.imDisplayRange(2,:));
    curRedPixelIntensities = mat2gray(handles.data.imageData{3}( curCellStats.PixelIdxList), handles.dataDisplay.imDisplayRange(3,:));

    gr_ratio = mat2gray(curGreenPixelIntensities ./ (eps + curRedPixelIntensities));        
    strCellDescription = sprintf( '%s\n\nGreed-Red Median Intensity Ratio: %.2f', ...
                                  strCellDescription, median( gr_ratio ));                          

    ratioOfMedianGreenRedIntensity = min(4096, median(curGreenPixelIntensities) / (eps + median(curRedPixelIntensities)));  
    strCellDescription = sprintf( '%s\n\nMedian-Green Median-Red Ratio: %.2f', ...
                                  strCellDescription, ratioOfMedianGreenRedIntensity);                          

    strCellDescription = sprintf( '%s\n\nMax Seed Strength: %.2f', strCellDescription, max(handles.data.imCellSeedPoints(curCellStats.PixelIdxList)) );
                              
    set(handles.LabelCellDescription, 'String', strCellDescription); 

    UpdateRedGreenDifferenceDisplay( handles );
    
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
    end
    imMultichannelOverlay( imMultichannelOverlay > 1 ) = 1;
    
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
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastOutputDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    else
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedNucleusDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    end

    if ~fileName 
        return;
    end
    
    annotationFile = fullfile(pathName, fileName);
    
    % load annotation data from file
    annotationData = load( annotationFile );
    
    if (isfield(annotationData, 'flagDataLoaded') && ~annotationData.flagDataLoaded) || ...
       (isfield(annotationData, 'flagRedGreenChannel') && ~annotationData.flagRedGreenChannel)
        errordlg( 'Invalid Annotation File' );
        return;
    end
    
    % extract all the data needed
    handles.data = [];

        % basic data
        handles.data.dataFilePath = annotationData.dataFilePath;
        handles.data.metadata = annotationData.metadata;
        handles.data.imageData = annotationData.imageData;
        handles.data.imRegValidMask = annotationData.imRegValidMask;

        if ~isfield( handles.data.metadata, 'channelColors' )
            handles.data.metadata.channelColors = [ 0 0 1; 0 1 0; 1 0 0 ];
        end
        handles = ComputeDisplayData(handles);
        
        % segmentation stuff
        handles.data.imLabelCellSeg = annotationData.imLabelCellSeg;        
        handles.data.CellSegColorMap = annotationData.CellSegColorMap;        
        handles.data.imCellSeedPoints = annotationData.imCellSeedPoints;        
        
        % display data
        handles.dataDisplay.imCellSeedPoints = imdilate(annotationData.imCellSeedPoints, ones(3,3,3));
        handles.dataDisplay.imCellSegRGBMask = label2rgbND( handles.data.imLabelCellSeg, handles.data.CellSegColorMap );
        
        % annotated stuff
        handles.data.cellStats = annotationData.cellStats;
        handles.data.cellPatternTypes = annotationData.cellPatternTypes;
        handles.data.unannotatedCellList = annotationData.unannotatedCellList;
        
        % parameters
        if isfield( annotationData, 'segAlgoParameters' )
            handles.data.segAlgoParameters = annotationData.segAlgoParameters;
        end        

        if isfield( annotationData, 'parameters' )
            handles.dataDisplay.oldParameters = annotationData.parameters;
        end        
        
    % data is now ready for display
    handles.flagDataLoaded = true;
    handles.flagReviewMode = true;    
    
    % change window name
    set( handles.CellPatternAnnotator, 'Name', sprintf( 'Cell Pattern Annotator - %s', annotationFile) );
        
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % initialize the state of UI controls
    set( handles.LabelAllChannels, 'String', 'All Channel Overlay' );
    
    % set cell pattern type list and cell class selector list
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );
    set( handles.poplistCellClassSelector, 'String', handles.data.cellPatternTypes );
    
    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % check if there are any annotated cells -- if not disable the ShowNextUnannotatedCell button
    if isempty( handles.data.unannotatedCellList )
        set( handles.ShowNextUnannotatedCell, 'Enable', 'off' );
    else
        set( handles.ShowNextUnannotatedCell, 'Enable', 'on' );
    end
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

% --------------------------------------------------------------------
function imageDataCorrected = CorrectDepthAttentuation( imageData )

    imageDataCorrected = imageData;
    for i = 2:numel(imageData)
        imageDataCorrected{i} = CorrectConfocalDepthSignalVariation_Opening( imageData{i}, 60 );        
    end

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
    [pathstr, name, ext] = fileparts( handles.data.dataFilePath{1} );
    outputDir = strtrim( fullfile(outputDir, name) );

    if ~isdir( outputDir )
        mkdir( outputDir );
    end
    
    data = handles.data;
    
    if handles.flagReviewMode 
        if isfield(handles.dataDisplay, 'oldParameters')
            data.parameters = handles.dataDisplay.oldParameters;
        end
    else
        data.parameters = handles.parameters;
    end
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save( fullfile(outputDir, 'CellPatternAnnotation.mat'), '-struct', 'data' );
    closeStatusDialog(h);
    
    clear data;

    % write summary file
    summary_fid = fopen( fullfile(outputDir, 'cellPatternAnnotationSummary.txt'), 'w' );    
    
        % print some information about the dataset
        fprintf( summary_fid, '\n>> Dataset Description:\n' );
        
        fprintf( summary_fid, '\n\tTwo-photon data file -- %s\n', handles.data.dataFilePath{1} );
    
        fprintf( summary_fid, '\n\tConfocal data file -- %s\n', handles.data.dataFilePath{2} );
        
        fprintf( summary_fid, '\n\tImage Size - [ %s ]\n', sprintf( ' %d ', handles.data.metadata.volSize) );
        fprintf( summary_fid, '\n\tImage Spacing - [ %s ]\n', sprintf( ' %.2f ', handles.data.metadata.voxelSpacing) );
    
        % print information about the annotation
        fprintf( summary_fid, '\n>> Annotation Summary:\n' );
        
        numTotalCells = numel( handles.data.cellStats );
        numCellsAnnotated = numel(find(flagCellPatternAnnotation));
        numWellSegmentedCells = numel(find([handles.data.cellStats.cellPatternId] > 4));
        
        fprintf( summary_fid, '\n\t%d cells were found by the segmentation algorithm\n', numTotalCells );
        
        fprintf( summary_fid, '\n\t%d/%d (%.2f%%) cells have been annotated\n', numCellsAnnotated, numTotalCells, 100*numCellsAnnotated/numTotalCells );
        
        cellPatternStats = cell( numel(handles.data.cellPatternTypes)+1, 4 );        
        cellPatternStats(1,1:end) = { 'CellPattern', 'Count', 'TotalPercentage', 'GoodCellPercentage' };
        
        for i = 1:numel(handles.data.cellPatternTypes)          

            numCellsInCurClass = numel(find([handles.data.cellStats.cellPatternId] == i));
            fprintf( summary_fid, '\n\t-- %s -- %d (%.2f%%) cells\n', ...
                     handles.data.cellPatternTypes{i}, ...
                     numCellsInCurClass, 100*numCellsInCurClass/numTotalCells );

            % note counts and percentages -- will be written to an excel file
            if i > 1 && numCellsAnnotated == numTotalCells
                
                % pattern name
                cellPatternStats{i, 1} = handles.data.cellPatternTypes{i};
                
                % raw count
                cellPatternStats{i, 2} = numCellsInCurClass;
                
                % percentage over all objects detected
                cellPatternStats{i, 3} = numCellsInCurClass * 100.0 / numTotalCells;
                
                % percentage over well-segmented cells
                if i > 4
                    cellPatternStats{i, 4} = numCellsInCurClass * 100.0 / numWellSegmentedCells;
                end
            
            end
            
            % get summary info of cells in this class
            if numCellsInCurClass > 2
                
                curClassCellStats = handles.data.cellStats([handles.data.cellStats.cellPatternId] == i);
                curClassPixelList = ismember( handles.data.imLabelCellSeg, find( [handles.data.cellStats.cellPatternId] == i ) );  
                for j = 1:3
                    
                    fprintf( summary_fid, '\n\t\tChannel-%d mean intensity:\n', j );
                    
                    fprintf( summary_fid, '\n\t\t\tMin-Max: [%.2f, %.2f]\n', ...
                                          min( handles.data.imageData{j}(curClassPixelList) ), ...
                                          max( handles.data.imageData{j}(curClassPixelList) ) );
                    fprintf( summary_fid, '\n\t\t\tMean-std: [%.2f, %.2f]\n', ...
                                          mean( handles.data.imageData{j}(curClassPixelList) ), ...
                                          std( double( handles.data.imageData{j}(curClassPixelList) ) ) );
                end
                
                fprintf( summary_fid, '\n\t\tCell Volume (cu um):\n' );
                fprintf( summary_fid, '\n\t\t\tMin-Max: [%.2f, %.2f]\n', min( [curClassCellStats.AreaPhysp] ), max( [curClassCellStats.AreaPhysp] ) );
                fprintf( summary_fid, '\n\t\t\tMean-std: [%.2f, %.2f]\n', mean( [curClassCellStats.AreaPhysp] ), std( [curClassCellStats.AreaPhysp] ) );

                fprintf( summary_fid, '\n\t\tFitted Ellipsoid Radii (um):\n' );
                ellipsoidRadius = cat( 1, curClassCellStats.ellipsoidRadiusPhysp );            
                fprintf( summary_fid, '\n\t\t\tMin Radii - [ %s ]\n', sprintf( ' %.2f ', min( ellipsoidRadius, [], 1 ) ) );
                fprintf( summary_fid, '\n\t\t\tMax Radii - [ %s ]\n', sprintf( ' %.2f ', max( ellipsoidRadius, [], 1 ) ) );
                fprintf( summary_fid, '\n\t\t\tMean Radii - [ %s ]\n', sprintf( ' %.2f ', mean( ellipsoidRadius, 1 ) ) );
                fprintf( summary_fid, '\n\t\t\tStddev Radii - [ %s ]\n', sprintf( ' %.2f ', std( ellipsoidRadius, 1 ) ) );

                fprintf( summary_fid, '\n\t\tFitted Ellipsoid Mean Radius (um):\n' );            
                meanEllipsoidRadius = mean( ellipsoidRadius, 2 );
                fprintf( summary_fid, '\n\t\t\tMin-Max: [%.2f, %.2f]\n', min( meanEllipsoidRadius ), max( meanEllipsoidRadius ) );
                fprintf( summary_fid, '\n\t\t\tMean-std: [%.2f, %.2f]\n', mean( meanEllipsoidRadius ), std( meanEllipsoidRadius ) );
                        
            end
            
        end             
        
        cellPatternStats(end,1:3) = { 'GoodSegmentation', numWellSegmentedCells, numWellSegmentedCells * 100.0 / numTotalCells };
        cellPatternStats(end+1,1:2) = { 'TotalCellsDetected', numTotalCells };
        cellPatternStats(end+2,1:2) = { 'Dataset', handles.data.dataFilePath{1} };
        
    fclose( summary_fid );
    
    if numCellsAnnotated == numTotalCells 
        xlswrite( fullfile(outputDir, 'cellPatternStats.xlsx'), cellPatternStats );
    end
    
    % ask if the user wants to save annotated cell images
    strQuestion = 'Do you want save images of the annotated cell patterns?';
    button = questdlg( strQuestion, 'Save annotated cell images', 'Yes', 'No', 'No' );
    flagSaveImages = strcmp( button, 'Yes' ); 
    
    h = waitbar(0, 'Saving Images of Annotated Cell Patterns ... Please Wait' );
    if flagSaveImages
        
        outputDirImages = fullfile(outputDir, 'images');
        if isdir( outputDirImages )
            rmdir( outputDirImages, 's' );                   
        end
        mkdir( outputDirImages );
        
        for i = 1:numel(handles.data.cellPatternTypes)
            mkdir( fullfile(outputDirImages, handles.data.cellPatternTypes{i}) );
        end

        szOutputImage = [100, 100];
        
        for cellId = 1:numel(handles.data.cellStats)

            curCellStats = handles.data.cellStats(cellId);
            curCellCentroid = curCellStats.Centroid;
            curCellBoundingBox = curCellStats.BoundingBox;
            curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

            curCellPatternType = handles.data.cellStats(cellId).cellPatternType;

            subinds = cell(1,3);
            imsize = size(handles.data.imageData{1});
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
            imCurCellCropped = double( handles.data.imageData{1}(subinds{1:2}, :) );
            imCurCellSegCropped = (handles.data.imLabelCellSeg(subinds{1:2}, :) == cellId);
            
            imCurCellMIP = imresize( mat2gray(max(imCurCellCropped .* imCurCellSegCropped, [], 3)), szOutputImage);

            imwrite( imCurCellMIP, fullfile(outputDirImages, curCellPatternType, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   
                 
            % Mid slices
            imCurCellSegMidSliceBndCropped = imresize( bwperim( imCurCellSegCropped(:, :, round(curCellCentroid(3))) ), szOutputImage, 'nearest'); 
            
            imCurCellMidSliceAllChannel = [];

            for chid = 1:3
                imCurChannelCropped = mat2gray( handles.data.imageData{chid}( subinds{1:2}, round(curCellCentroid(3)) ), handles.dataDisplay.imDisplayRange(chid,:) );
                imCurChannelCropped = imresize(imCurChannelCropped, szOutputImage);
                imCurCellMidSliceAllChannel = cat( 3, imCurCellMidSliceAllChannel, imCurChannelCropped );
            end

            channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];

            imwrite( genImageMaskOverlay(imCurCellMidSliceAllChannel(:,:,1), imCurCellSegMidSliceBndCropped, [1, 0, 0], 0.5), ...
                     fullfile(outputDirImages, curCellPatternType, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   

            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
                     fullfile(outputDirImages, curCellPatternType, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   

            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
                     fullfile(outputDirImages, curCellPatternType, sprintf('CellMidSliceAllChannel_%.3d.png', cellId)), 'png' );   

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

    close(handles.CellPatternAnnotator);
        
% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to file_close CellPatternAnnotator.
function CellPatternAnnotator_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
    
    button = questdlg('Are you sure you want to close?', 'Closing Annotation Tool', 'Yes', 'No', 'No');
    if strcmp( button, 'No' )
        return;
    end
    
    % file_close
    delete(hObject);
    
function PrettyPrintStepDescription( strStepDescription )

    strStar = strStepDescription;
    strStar(:) = '*';
    strStar = [ '****', strStar, '****' ];
    fprintf( '\n\n%s', strStar );
    fprintf( '\n    %s    ', strStepDescription );
    fprintf( '\n%s\n\n', strStar );

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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function LabelCellDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelCellDescription (see GCBO)
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

    handles.parameters.segmentation = NucleiSegmentationParametersGUI( handles.parameters.segmentation );

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
% --- Executes during object deletion, before destroying properties.
function CellPatternAnnotator_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    if isfield(handles, 'history')
        history = handles.history;
        [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
        historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
        save( historyFile, '-struct', 'history' );  
    end

%     if ~isempty(handles.imarisApp) && isobject(handles.imarisApp)
%         delete( handles.imarisApp );
%     end       
    
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
function [isInside] = IsMouseInsideAxes(hAxes)

    mousePos = get(hAxes,'CurrentPoint');  % The current point w.r.t the figure.
    xylim = [get(hAxes, 'xlim'); get(hAxes, 'ylim')];
    
    if any( mousePos(1,1:2) < xylim(1:2) ) || any( mousePos(1,1:2) > xylim(3:4) )
        isInside = false;
    else
        isInside = true;
    end

% --------------------------------------------------------------------
function FnSliceScroll_Callback(hSrc, eventdata)

    handles = guidata(hSrc);
    
    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{1});

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

function EditboxSliceId_Callback(hObject, eventdata, handles)
% hObject    handle to EditboxSliceId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditboxSliceId as text
%        str2double(get(hObject,'String')) returns contents of EditboxSliceId as a double


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


% --- Executes on button press in CheckboxShowNucleusChannel.
function CheckboxShowNucleusChannel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxShowNucleusChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxShowNucleusChannel

    % change log use mode
    handles.flagShowNucleusChannel = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes on button press in CheckboxShowGreenChannel.
function CheckboxShowGreenChannel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxShowGreenChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxShowGreenChannel

    % change log use mode
    handles.flagShowGreenChannel = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes on button press in CheckboxShowRedChannel.
function CheckboxShowRedChannel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxShowRedChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxShowRedChannel

    % get value of red-channel flag
    handles.flagShowRedChannel = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


% --- Executes on button press in CheckboxLocalMIP.
function CheckboxLocalMIP_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxLocalMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxLocalMIP

    % get value of red-channel flag
    handles.flagShowLocalMIP = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

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
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

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


% --------------------------------------------------------------------
function Inspection_View_Full_Seg_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to Inspection_View_Full_Seg_In_Imaris (see GCBO)
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

    if ~isempty(handles.imarisApp) && isobject(handles.imarisApp)
        delete( handles.imarisApp );
    end       

    % compute isosurface geometry for cells in each pattern
    hStatusDlg = waitbar( 0, 'computing surface geometry for cells in each pattern' );
    cellSurfaceObjectList = {};
    
    for cid = 1:numel(handles.data.cellStats)
        
        waitbar( cid/numel( handles.data.cellStats ), hStatusDlg );
        
        % name
        curCellIsoSurface.name = sprintf( 'CellSeg_%d', cid);
        
        % color
        curCellIsoSurface.color = handles.data.CellSegColorMap(cid, :);
        
        % compute surface geometry
        curCellStats = handles.data.cellStats( cid );        
        curCellCentroid = round( curCellStats.Centroid );    

        curCellBoundingBox = curCellStats.BoundingBox;
        curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );
       
        % create crop indices
        subinds = cell(1,3);
        imsize = size(handles.data.imageData{1});
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
        
        % crop segmentation mask
        imCurCellSegCropped = padarray( double(handles.data.imLabelCellSeg(subinds{:}) == cid), ones(1,3), 0 );                                        
        imCurCellSegSmoothed = smooth3( imCurCellSegCropped );            
        curCellSurfaceGeometry = isosurface( imCurCellSegSmoothed, 0.5 );
        curCellSurfaceGeometry.normals = isonormals( imCurCellSegSmoothed, curCellSurfaceGeometry.vertices );
        
        % correct vertex positions by adding offset
        curCellSurfaceGeometry.vertices(:,1) = subinds{2}(1) - 1 + curCellSurfaceGeometry.vertices(:,1);
        curCellSurfaceGeometry.vertices(:,2) = subinds{1}(1) - 1 + curCellSurfaceGeometry.vertices(:,2);
        curCellSurfaceGeometry.vertices(:,3) = subinds{3}(1) - 1 + curCellSurfaceGeometry.vertices(:,3);
        
        % add cell surface to cell pattern surface object
        curCellIsoSurface.surfaces = curCellSurfaceGeometry;
        cellSurfaceObjectList{cid} = curCellIsoSurface;
        
    end
    closeStatusDialog( hStatusDlg );
    
    % get cell seed point locations
    stats = regionprops( bwlabeln( handles.data.imCellSeedPoints ), 'Centroid' );
    cellSeedPointLocations = cat( 1, stats.Centroid );

    % Display everything in imaris
    handles.imarisApp = DisplayMultichannel3DDataInImaris( handles.data.imageData, ...
                                                           'spacing', handles.data.metadata.voxelSpacing, ...
                                                           'spotLocations', cellSeedPointLocations, ...
                                                           'spotRadius', 3, ...
                                                           'surfaceObjects', cellSurfaceObjectList, ...
                                                           'displayRanges', handles.dataDisplay.imDisplayRange, ...
                                                           'displayColors', handles.data.metadata.channelColors );
    % Update handles structure
    guidata(hObject, handles);                                                    


% --------------------------------------------------------------------
function Inspection_Callback(hObject, eventdata, handles)
% hObject    handle to Inspection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Inspection_View_Cell_Seg_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to Inspection_View_Cell_Seg_In_Imaris (see GCBO)
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
    
    if ~isempty(handles.imarisAppCellSegCropped) && isobject(handles.imarisAppCellSegCropped)
        delete( handles.imarisAppCellSegCropped );
    end       
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

    % create crop indices
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{1});
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

    % crop cell bounding box from whole volume
    imCellCropped = cell(1,3);
    for i = 1:3
        %imCellCropped{i} = handles.data.imageData{i}(subinds{:}) .* (handles.data.imLabelCellSeg( subinds{:} ) == handles.dataDisplay.curCellId);
        imCellCropped{i} = handles.data.imageData{i}(subinds{:});
    end

    % crop cell segmentation mask
    imCellSegCropped = handles.data.imLabelCellSeg( subinds{:} );
    imCellSegCropped = padarray( imCellSegCropped, ones(1,3) );
    imCellSegCropped = double( imCellSegCropped == handles.dataDisplay.curCellId );    
    
    % create isosurface
    imCellSegSmoothed = smooth3( imCellSegCropped );
    curCellSurfaceGeometry = isosurface( imCellSegSmoothed, 0.5 );
    curCellSurfaceGeometry.normals = isonormals( imCellSegSmoothed, curCellSurfaceGeometry.vertices );
    cellIsoSurface.surfaces = curCellSurfaceGeometry;
    cellIsoSurface.name = curCellStats.cellPatternType;
    cellIsoSurface.color = MapCellPatternToColor( curCellStats.cellPatternType );
    
    % get seed points
    imCellSeedCropped = padarray( handles.data.imCellSeedPoints( subinds{:} ), ones(1,3), 0 );
    imCellSeedCropped( ~imCellSegCropped ) = 0;
    stats = regionprops( bwlabeln( imCellSeedCropped ), 'Centroid' );
    cellSeedPointLocations = cat( 1, stats.Centroid );
    
    % display in imaris
    curCellDisplayRange = handles.dataDisplay.imDisplayRange;
    curCellDisplayColor = handles.data.metadata.channelColors;
    if numel(stats) > 0
        
        handles.imarisAppCellSegCropped = DisplayMultichannel3DDataInImaris( imCellCropped, ...
                                                                             'spacing', handles.data.metadata.voxelSpacing, ...
                                                                             'displaycolors', curCellDisplayColor, ...
                                                                             'displayranges', curCellDisplayRange, ...
                                                                             'spotLocations', cellSeedPointLocations, ...
                                                                             'spotRadius', 3, ...                                                                      
                                                                             'surfaceObjects', cellIsoSurface );
    else
        
        handles.imarisAppCellSegCropped = DisplayMultichannel3DDataInImaris( imCellCropped, ...
                                                                             'spacing', handles.data.metadata.voxelSpacing, ...
                                                                             'displaycolors', curCellDisplayColor, ...
                                                                             'displayranges', curCellDisplayRange, ...
                                                                             'surfaceObjects', cellIsoSurface );
    end
    
    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Inspect_Check_Seed_Detection_Callback(hObject, eventdata, handles)
% hObject    handle to Inspect_Check_Seed_Detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~handles.flagDataLoaded
        return;
    end
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

    % create crop indices
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{1});
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
    
    % crop cell bounding box from whole volume
    imCellCropped = handles.data.imageData{1}( subinds{:} );
    
    % crop the segmentation mask 
    imCellSegCropped = handles.data.imLabelCellSeg( subinds{:} );
    
    % Call Seed Point Detection GUI
    if isfield( handles.data, 'segAlgoParameters' )
        cellDiameterRange = handles.data.segAlgoParameters.cellDiameterRange;
    else
        cellDiameterRange = [12, 20];
    end
    
    BlobDetectionGUI( imCellCropped, cellDiameterRange, ...
                      'foregroundMask', double( imCellSegCropped > 0 ), ...  
                      'spacing', handles.data.metadata.voxelSpacing );
    

% --------------------------------------------------------------------
function [ cellColor ] = MapCellPatternToColor( cellPatternType )

    switch cellPatternType

        case { 'Under_Segmentation', 'Over_Segmentation' }

            cellColor = [1.0, 1.0, 1.0, 0.5]; % gray

        case { 'Mono_Pre_G1', 'Multi_Pre_G1' } 

            cellColor = [0.0, 0.0, 1.0, 0.5]; % blue

        case { 'Mono_G1', 'Multi_G1' } 

            cellColor = [1.0, 0.0, 0.0, 0.5]; % red

        case { 'Mono_S', 'Multi_S' } 

            cellColor = [1.0, 1.0, 0.0, 0.5]; % yellow

        case { 'Mono_G2', 'Multi_G2' } 

           cellColor = [0.0, 1.0, 0.0, 0.5]; % green

        case { 'Mono_Prophase', 'Multi_Prophase' }     

            cellColor = [0.0, 1.0, 1.0, 0.5]; % green

        otherwise

            cellColor = [0.25, 0.20, 0.32, 0.5]; % purple
            
    end        

% --------------------------------------------------------------------
function Inspect_View_Cell_Patterns_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to Inspect_View_Cell_Patterns_In_Imaris (see GCBO)
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

    if ~isempty(handles.imarisApp) && isobject(handles.imarisApp)
        delete( handles.imarisApp );
    end           
    
    imageDataPadded = cell(1,3);
    for i = 1:3
        imageDataPadded{i} = padarray( handles.data.imageData{i}, ones(1,3), min(handles.data.imageData{i}(:)) );
    end
    
    % compute isosurface geometry for cells in each pattern
    hStatusDlg = waitbar( 0, 'computing surface geometry for cells in each pattern' );
    cellPatternSurfaceObjectList = {};
    
    for pid = 1:numel( handles.data.cellPatternTypes )        
        
        waitbar( pid/numel( handles.data.cellPatternTypes ), hStatusDlg );
        
        % ignore unannotated cells and errors
        if ismember( handles.data.cellPatternTypes{pid}, { 'None', 'Bad_Detection' } )
           continue; 
        end
        
        % get all the cells with the current pattern id
        curPatternCellIds = find( [handles.data.cellStats.cellPatternId] == pid );        
        
        % check if there any cells in this pattern
        if isempty( curPatternCellIds )
           continue; 
        end
        
        % name
        curPatternIsoSurface.name = handles.data.cellPatternTypes{pid};
        
        % color
        curPatternIsoSurface.color = MapCellPatternToColor( handles.data.cellPatternTypes{pid} );
        
        % compute surface geometry of each cell of current pattern        
        curPatternIsoSurface.surfaces = [];
        for cid = curPatternCellIds
            
            curCellStats = handles.data.cellStats( cid );        
            curCellCentroid = round( curCellStats.Centroid );    

            curCellBoundingBox = curCellStats.BoundingBox;
            curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );

            % check if pattern is same
            if curCellStats.cellPatternId ~= pid
                error( 'ERROR - pattern of cell and surface object dont match' );
            end
            
            % create crop indices
            subinds = cell(1,3);
            imsize = size(handles.data.imageData{1});
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
            
            % crop segmentation mask
            imCurCellSegCropped = padarray( double(handles.data.imLabelCellSeg(subinds{:}) == cid), ones(1,3), 0 );                                        
            imCurCellSegSmoothed = smooth3( imCurCellSegCropped );            
            curCellSurfaceGeometry = isosurface( imCurCellSegSmoothed, 0.5 );
            curCellSurfaceGeometry.normals = isonormals( imCurCellSegSmoothed, curCellSurfaceGeometry.vertices );
            
            % correct vertex positions by adding offset
            curCellSurfaceGeometry.vertices(:,1) = subinds{2}(1) - 1 + curCellSurfaceGeometry.vertices(:,1);
            curCellSurfaceGeometry.vertices(:,2) = subinds{1}(1) - 1 + curCellSurfaceGeometry.vertices(:,2);
            curCellSurfaceGeometry.vertices(:,3) = subinds{3}(1) - 1 + curCellSurfaceGeometry.vertices(:,3);
            
            % add cell surface to cell pattern surface object
            curPatternIsoSurface.surfaces = [ curPatternIsoSurface.surfaces; curCellSurfaceGeometry ];
            
        end
        
        % add cell pattern surface to the list
       cellPatternSurfaceObjectList{end+1} = curPatternIsoSurface;
       
    end    
    closeStatusDialog( hStatusDlg );
    
    % display cell pattern distribution in imaris
    handles.imarisAppCellSegCropped = DisplayMultichannel3DDataInImaris( imageDataPadded, ...
                                                                         'spacing', handles.data.metadata.voxelSpacing, ...
                                                                         'surfaceObjects', cellPatternSurfaceObjectList, ...
                                                                         'displayRanges', handles.dataDisplay.imDisplayRange, ...
                                                                         'displayColors', handles.data.metadata.channelColors );
    
    % Update handles structure
    guidata(hObject, handles);
    


% --- Executes on button press in CheckboxInvalidROIMask.
function CheckboxValidROIMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxInvalidROIMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxInvalidROIMask

    handles.flagShowInvalidROIMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function CellPatternAnnotator_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{1});
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
        set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    end
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    