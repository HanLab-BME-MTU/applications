function varargout = CellStateAnalyzer(varargin)
%CELLSTATEANALYZER M-file for CellStateAnalyzer.fig
%      CELLSTATEANALYZER, by itself, creates a new CELLSTATEANALYZER or raises the existing
%      singleton*.
%
%      H = CELLSTATEANALYZER returns the handle to a new CELLSTATEANALYZER or the handle to
%      the existing singleton*.
%
%      CELLSTATEANALYZER('Property','Value',...) creates a new CELLSTATEANALYZER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CellStateAnalyzer_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CELLSTATEANALYZER('CALLBACK') and CELLSTATEANALYZER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CELLSTATEANALYZER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellStateAnalyzer

% Last Modified by GUIDE v2.5 15-Apr-2014 12:12:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellStateAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @CellStateAnalyzer_OutputFcn, ...
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

% --- Executes just before CellStateAnalyzer is made visible.
function CellStateAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
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
    
    % Choose default command line output for CellStateAnalyzer
    handles.output = hObject;

    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellStateAnalyzerHistory.mat' );
    if exist( historyFile, 'file' )
       handles.history = load( historyFile );
    else
       handles.history.lastAnalyzedNucleusDir = pathstr;
       handles.history.lastAnalyzedRedGreenDir = pathstr;
    end

    % Initialize global variables
    handles.cellDisplaySize = 70;     
    
    handles.flagDataLoaded = false;
    handles.flagUseLOG = get(handles.CheckboxGlobalLog, 'Value');
    handles.flagShowGlobalSegMask = get(handles.CheckboxGlobalSegMask, 'Value');
    handles.flagShowCellBBox = get(handles.CheckboxCellBBox, 'Value');
    
    handles.defaultCellPatternTypes = cellstr( get(handles.ListboxCellPatternSelector, 'String') );    
    set( handles.poplistPredictedCellClassSelector, 'String', {handles.defaultCellPatternTypes{:}, 'Any'} );
    
    handles.imarisAppCellSeg = [];
    handles.imarisAppCellPattern = [];
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

    handles.defaultParameters.classification = CellCycleStateIdentificationParametersGUI( 'default' );
    f = rdir( fullfile(fileparts(mfilename('fullpath')), '**', 'G1_S_G2_M.model') );
    if ~isempty(f)
        handles.defaultParameters.classification.flagPerformCellCycleStateIdentification = true;
        handles.defaultParameters.classification.cellCycleStateIdentificationModelDir = fileparts(f(1).name);
    end
    
    handles.defaultParameters.flagParallelize = PARAMETERS.flagParallelize;
    handles.defaultParameters.flagDebugMode = PARAMETERS.flagDebugMode;
    
    handles.parameters = handles.defaultParameters;
    
    % no groundtruth
    set(handles.poplistGroundtruthCellClassSelector, 'Enable', 'off' );               

    % open matlab pool for parallel processing    
    handles.parameters.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.parameters.flagParallelize && ~handles.parameters.flagPoolOpenedAlready 
        matlabpool open;
    end
    
    % make sure weka is in the path    
    AddWekaClassesToPath()
    
    % Set callbacks
    set(handles.CellStateAnalyzer, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes CellStateAnalyzer wait for user response (see UIRESUME)
    % uiwait(handles.CellStateAnalyzer);

% --- Outputs from this function are returned to the command line.
function varargout = CellStateAnalyzer_OutputFcn(hObject, eventdata, handles)
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
    
    % ask the use to select the oif file for the red/green channel
    [fileName,pathName] = uigetfile( fullfile( pathName, '..', '*.oif; *.oib' ), ...
                                              'Select the data file containing the FUCCI Cell Cycle Reporter' );   
        
    if ~fileName 
        return;
    end
    
    dataFilePath{2} = fullfile( pathName, fileName );        
    handles.history.lastAnalyzedRedGreenDir = pathName;
    
    % load nucleus channel data
    PrettyPrintStepDescription( 'Loading Histone Data' );
    
    hStatusDialog = waitbar(0, 'Loading Histone Data');
    try        
        imageSeriesNucleus = loadIntravitalDataset( dataFilePath{1}  );    
        if (imageSeriesNucleus(1).metadata.voxelSpacing(1)/ imageSeriesNucleus(1).metadata.voxelSpacing(3)) >= 1
            imageSeriesNucleus(1).metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
        end        
        imageSeriesNucleus(1).metadata.channelColors = [ 0 , 0 , 1 ];
        metadata_NucleusChannel = imageSeriesNucleus(1).metadata
    catch err
        fprintf( 'ERROR: could not load nucleus channel data from file %s', dataFilePath{1} );
        return;
    end
    
    % load red/green channel data
    PrettyPrintStepDescription( 'Loading Fucci Cell Cycle Reporter Data' );

    waitbar(0, hStatusDialog, 'Loading Fucci Cell Cycle Reporter Data');
    try

        imageSeriesRedGreen = loadIntravitalDataset( dataFilePath{2}  );    

        if imageSeriesRedGreen(1).metadata.channelExcitationWavelength(1) > imageSeriesRedGreen(1).metadata.channelExcitationWavelength(2)

            % first channel is not green so swap it 
            warning( 'GFP was not the first channel in confocal data. It will be swapped to the first place');
            imageSeriesRedGreen(1).imageData = imageSeriesRedGreen(1).imageData(:,[2,1]);
            imageSeriesRedGreen(1).metadata.channelExcitationWavelength = imageSeriesRedGreen(1).metadata.channelExcitationWavelength([2,1]);

        end

        imageSeriesRedGreen(1).metadata.channelColors = [ 0 1 0; 1 0 0 ];

        metadata_FUCCI = imageSeriesRedGreen(1).metadata
        
    catch err

        err
        fprintf( 'ERROR: could not load red-green channel data from file %s', dataFilePath{2} );
        errordlg( sprintf( 'Could not load red-green channel data from file %s', dataFilePath{2} ) ); 
        return;
    end
       
    % perform checks on data
    if any( imageSeriesRedGreen(1).metadata.numChannels ~= 2  )
        errordlg( sprintf( 'Red-Green channel data is expected to contain 2 channels. The file selected by you contains %d channels', imageSeriesRedGreen(1).metadata.numChannels ) );
        return;
    end

    if any( imageSeriesRedGreen(1).metadata.volSize ~= imageSeriesNucleus(1).metadata.volSize )
        errordlg( sprintf( 'Volume Size of Nucleus data doesnt match with red-green channel data' ) );
        return;
    end

    if (imageSeriesRedGreen(1).metadata.voxelSpacing(1)/ imageSeriesRedGreen(1).metadata.voxelSpacing(3)) >= 1
        imageSeriesRedGreen(1).metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
    end

    if any( imageSeriesRedGreen(1).metadata.voxelSpacing ~= imageSeriesNucleus(1).metadata.voxelSpacing )
        errordlg( sprintf( 'Pixel spacing of nucleus data doesnt match with red-green channel data' ) );
        return;
    end

    % store image data in handles structures    
    handles.flagDataLoaded = true;
    handles.data = [];
    handles.data.dataFilePath = dataFilePath;
    handles.data.metadata = imageSeriesNucleus(1).metadata;
    
    set( handles.CellStateAnalyzer, 'Name', sprintf( 'Cell State Analyzer - %s', handles.data.dataFilePath{1} ) );
    
    handles.data.metadata.channelExcitationWavelength = [ imageSeriesNucleus(1).metadata.channelExcitationWavelength, ...
                                                          imageSeriesRedGreen(1).metadata.channelExcitationWavelength ];

    handles.data.metadata.channelColors = cat( 1 , imageSeriesNucleus(1).metadata.channelColors, ...
                                                   imageSeriesRedGreen(1).metadata.channelColors );             
        
    % correct stage-shift        
    PrettyPrintStepDescription( 'Correcting shift between confocal and two-photon data' );

    waitbar(0, hStatusDialog, 'Correcting shift between confocal and two-photon data');
    [handles.data.imageData, ...
     handles.data.imRegValidMask] = CorrectTwoPhotonConfocalStageShift( imageSeriesNucleus(1).imageData, imageSeriesNucleus(1).metadata.voxelSpacing, ...
                                                                        imageSeriesRedGreen(1).imageData, imageSeriesRedGreen(1).metadata.voxelSpacing, ...
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
    if handles.parameters.segmentation.flagPerformRegionMerging
        regionMergingModelFile = handles.parameters.segmentation.regionMergingModelFile;
    end
    
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints] = segmentCellsInIntravitalData( handles.data.imageData{1}, ...
                                                                    handles.data.metadata.voxelSpacing, ...                                                                      
                                                                    'flagParallelize', handles.parameters.flagParallelize, ...
                                                                    'flagDebugMode', handles.parameters.flagDebugMode, ...
                                                                    'cellDiameterRange', handles.parameters.segmentation.cellDiameterRange, ...
                                                                    'thresholdingAlgorithm', 'MinErrorPoissonSliceBySliceLocal', ...
                                                                    'seedPointDetectionAlgorithm', handles.parameters.segmentation.seedPointDetectionAlgorithm, ...
                                                                    'minCellVolume', handles.parameters.segmentation.minCellVolume, ...
                                                                    'flagIgnoreCellsOnXYBorder', handles.parameters.segmentation.flagIgnoreXYBorderCells, ...
                                                                    'roiMask', handles.data.imRegValidMask, ...
                                                                    'minCellROIOverlap', handles.parameters.classification.minCellROIOverlap, ...
                                                                    'regionMergingModelFile', regionMergingModelFile);

    [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.dataDisplay.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints, ones(3,3,3) );
    closeStatusDialog(hStatusDialog);
    
    % compute cell properties
    handles.data.cellStats = ComputeCellProperties( handles );

    % assign all cells to other class initially
    handles.data.cellPatternTypes = handles.defaultCellPatternTypes;
    
    for i = 1:numel(handles.data.cellStats)
        handles.data.cellStats(i).cellPatternId = 1; %Not Annotated
        handles.data.cellStats(i).cellPatternType = handles.data.cellPatternTypes{1};
    end
    
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % set cell pattern types
    set( handles.poplistPredictedCellClassSelector, 'String', handles.data.cellPatternTypes );
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );    
    
    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    set(handles.poplistPredictedCellClassSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % if groundtruth was present remove it
    if isfield( handles.data, 'groundtruth' )
       handles.data = rmfield( handles.data, 'groundtruth');
       set(handles.poplistGroundtruthCellClassSelector, 'Enable', 'off' );           
    end

    % Update handles structure
    guidata(hObject, handles);
   
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

    % run cell cycle state identification if requested
    if handles.parameters.classification.flagPerformCellCycleStateIdentification
        PerformCellCycleStateIdentification(hObject, handles);
    end
    
% --------------------------------------------------------------------
function PerformCellCycleStateIdentification(hObject, handles)

    if ~handles.flagDataLoaded
        return;
    end

    % Ask the user to select the output directory
    modelDir = handles.parameters.classification.cellCycleStateIdentificationModelDir;

    if isempty(modelDir)
        errordlg( 'Could not perform cell cycle identification. Model was not selected. Please go to File->Set Parameters->Cell Cycle State Identification and select the model directory.' );
        return;
    end
    
    modelFileName = [handles.parameters.classification.cellCycleStateIdentificationModelType, '.model'];
    
    if ~exist( fullfile(modelDir, modelFileName), 'file' )
        errordlg( sprintf('Could not perform cell cycle identification. Selected model directory does not contain a model of the type selected. Please make sure the selected model directory contains a file named - %s', modelFileName) );
        return;
    end
    
    % make sure weka is in claspath
    AddWekaClassesToPath();
    
    PrettyPrintStepDescription( 'Identifying cell cycle state of each cell' );
    
    % load model
    cellCycleModelClass = @CellCycleStateClassifier_OneStage;
    cellCycleModel = cellCycleModelClass( fullfile(modelDir, modelFileName) );
    
    % pre-process data
    fprintf('\n>> Pre-processing the image data ... \n');       
    
    imageDataAdjusted = CellCycleStateClassifier.preprocessImageData( handles.data.imageData );
    
    % apply classification model
    numCells = numel(handles.data.cellStats);    
    fprintf('\n>> Applying classification model on each of the %d cells ...\n', numCells);       
    
    predictedClassLabels = cell(numCells, 1);
    classPredictionProbabilities = cell(numCells, 1);  
    cellFeatureStruct = cell( numCells, 1);
    modelClassNameList = cellCycleModel.getClassNameList()   
    
    handles.data.cellPatternTypes = modelClassNameList;
    
    hStatusDlg = waitbar( 0, 'Applying classification model on each cell' );
    
    for cellId = 1:numCells

        [predictedClassLabels{cellId}, ...
         classPredictionProbabilities{cellId}, ...
         cellFeatureStruct{cellId}] = cellCycleModel.predictCell(imageDataAdjusted, ...
                                                                  handles.data.imRegValidMask, ...
                                                                  handles.data.imLabelCellSeg, cellId, ...
                                                                  handles.data.metadata.voxelSpacing); 

        curPatternId = find( strcmpi(handles.data.cellPatternTypes, predictedClassLabels{cellId}) );
        
        handles.data.cellStats(cellId).cellPatternId = curPatternId; 
        handles.data.cellStats(cellId).cellPatternType = handles.data.cellPatternTypes{curPatternId};

        waitbar( cellId/numCells, hStatusDlg );
        
    end
    
    closeStatusDialog(hStatusDlg);
    
    % print class distribution
    fprintf( '\nPredicted Class Distribution: \n' );
    
    tabulate( predictedClassLabels );
    
    fprintf( '\nTotal number of cells: %d\n', numel(handles.data.cellStats) );
    
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % set cell pattern types
    set( handles.poplistPredictedCellClassSelector, 'String', handles.data.cellPatternTypes );
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );    
    
    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    set(handles.poplistPredictedCellClassSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
   
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

function [cellStats] = ComputeCellProperties( handles )

    hStatusDialog = waitbar(0, 'Computing properties of segmented cells');
    
    cellStats = regionprops( handles.data.imLabelCellSeg, ...
                             'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
    
    for i = 1:numel(cellStats)

        % intensity descriptors
        cellPixelIntensities = handles.data.imageData{1}( cellStats(i).PixelIdxList );
        cellStats(i).meanIntensity = mean( cellPixelIntensities );
        cellStats(i).stdIntensity = std( double( cellPixelIntensities ) );
        cellStats(i).minIntensity = min( cellPixelIntensities );
        cellStats(i).maxIntensity = max( cellPixelIntensities );
        
        % volume
        cellStats(i).AreaPhysp = cellStats(i).Area * prod( handles.data.metadata.voxelSpacing );
        
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
    
    if ~handles.flagDataLoaded
        return;
    end

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

    nhist( channelPixelIntensties, ...
           'nolegend', 'noerror', 'pdf', 'smooth', ...
           'minx', 0, 'maxx', maxval, ...
           'color', handles.data.metadata.channelColors );
        
    
    hold off;
    
    legend( {'CFP', 'GFP', 'RFP'}, 'Location', 'best', 'FontSize', 8.0 );
    title( 'Channel Histograms of Cell Pixels', 'FontSize', 11.0 ); 
    
% --------------------------------------------------------------------
function UpdateCellDescriptors(handles)

    if ~handles.flagDataLoaded
        return;
    end
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );
    curCellCentroid = round(curCellStats.Centroid);
    
    % display cell information
    strCellDescription = sprintf( '' );
    if isfield( handles.data, 'groundtruth' )
        curAnnotLabel = handles.data.groundtruth.cellStats(handles.dataDisplay.curCellId).cellPatternType;
        strCellDescription = sprintf( '%s\nGroundtruth Cell Type: %s', strCellDescription, curAnnotLabel);
    end    
    strCellDescription = sprintf( '%s\n\nVolume (cu um): %.2f', strCellDescription, curCellStats.AreaPhysp );
    strCellDescription = sprintf( '%s\n\nMean-std Intensity: [%.2f, %.2f]', strCellDescription, curCellStats.meanIntensity, curCellStats.stdIntensity );
    strCellDescription = sprintf( '%s\n\nIntensity Range: [%d, %d]', strCellDescription, curCellStats.minIntensity, curCellStats.maxIntensity );
    strCellDescription = sprintf( '%s\n\nCentroid: [%d, %d, %d]', strCellDescription, curCellCentroid(1), curCellCentroid(2), curCellCentroid(3) );
                  
    % display difference between pixel intensities of red ang green channels 
    curGreenPixelIntensities = mat2gray(handles.data.imageData{2}( curCellStats.PixelIdxList), handles.dataDisplay.imDisplayRange(2,:));
    curRedPixelIntensities = mat2gray(handles.data.imageData{3}( curCellStats.PixelIdxList), handles.dataDisplay.imDisplayRange(3,:));

    gr_ratio = mat2gray(curGreenPixelIntensities ./ (eps + curRedPixelIntensities));        
    strCellDescription = sprintf( '%s\n\nGreed-Red Median Intensity Ratio: %.2f', ...
                                  strCellDescription, median( gr_ratio ));                          

    ratioOfMedianGreenRedIntensity = median(curGreenPixelIntensities) / median(curRedPixelIntensities);  
    strCellDescription = sprintf( '%s\n\nMedian-Green Median-Red Ratio: %.2f', ...
                                  strCellDescription, ratioOfMedianGreenRedIntensity);                          

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
    
function [ imMaskOverlay ] = genImageRGBMaskOverlay( im, rgbMask, maskAlpha )

    imMaskOverlay = repmat( mat2gray(im), [1,1,3] );
    blnMask = repmat( max( rgbMask, [], 3 ) > 0, [1, 1, 3] );
    imMaskOverlay(blnMask) = (1 - maskAlpha) * imMaskOverlay(blnMask) + maskAlpha * rgbMask(blnMask);
    imMaskOverlay( imMaskOverlay > 1 ) = 1;
    
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 98.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
% --------------------------------------------------------------------
function File_Load_Annotation_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % get annotation file from user
    if isfield( handles.history, 'lastAnnotationDir' )
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnnotationDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    else
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalyzedNucleusDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    end

    if ~fileName 
        return;
    end

    hStatusDialog = waitbar(0, 'Loading selected annotation file ...');
    
    handles.history.lastAnnotationDir = pathName;
    annotationFile = fullfile(pathName, fileName);
    
    % load annotation data from file
    annotationData = load( annotationFile );
    
    % retrieve data needed for this tool
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
        handles.data.imCellSeedPoints = annotationData.imCellSeedPoints;
        handles.data.CellSegColorMap = annotationData.CellSegColorMap;

        handles.dataDisplay.imCellSeedPoints = annotationData.imCellSeedPoints;
        handles.dataDisplay.imCellSegRGBMask = label2rgbND( handles.data.imLabelCellSeg, handles.data.CellSegColorMap );
        
        % cell pattern annotation
        handles.data.groundtruth.cellStats = annotationData.cellStats;
        handles.data.groundtruth.cellPatternTypes = annotationData.cellPatternTypes;
        
        % parameters
        if isfield( annotationData, 'segAlgoParameters' )
            
            fieldsNeeded = { 'cellDiameterRange', ...
                             'minCellVolume', ...
                             'seedPointDetectionAlgorithm', ...
                             'thresholdingAlgorithm' };
                         
            for i = 1:numel(fieldsNeeded)            
                if isfield( annotationData.segAlgoParameters, fieldsNeeded{i} )
                    handles.data.groundtruth.parameters.(fieldsNeeded{i}) = annotationData.segAlgoParameters.(fieldsNeeded{i});
                end                
            end
            
            if isfield( annotationData.segAlgoParameters, 'minObjectSize' )
                handles.data.groundtruth.parameters.minCellVolume = annotationData.segAlgoParameters.minObjectSize;
            end
            
        end        
        
    % compute cell properties    
    handles.data.cellStats = ComputeCellProperties( handles );
    
    % assign all cells to other class initially
    handles.data.cellPatternTypes = handles.defaultCellPatternTypes;
    
    for i = 1:numel(handles.data.cellStats)

        handles.data.cellStats(i).cellPatternId = 1; 
        handles.data.cellStats(i).cellPatternType = handles.data.cellPatternTypes{1};

    end
    
    % data ready for display
    handles.flagDataLoaded = true;
    
    % change window name
    set( handles.CellStateAnalyzer, 'Name', sprintf( 'Cell State Analyzer - %s', handles.data.dataFilePath{1} ) );
    
    % set cell pattern type list and cell class selector list
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );
    set( handles.poplistPredictedCellClassSelector, 'String', {handles.data.cellPatternTypes{:}, 'Any'} );
    
    set(handles.poplistGroundtruthCellClassSelector, 'Enable', 'on' );               
    set(handles.poplistGroundtruthCellClassSelector, 'String', {handles.data.groundtruth.cellPatternTypes{:}, 'Any'} );
    set(handles.poplistGroundtruthCellClassSelector, 'Value', numel(handles.data.groundtruth.cellPatternTypes) + 1 );
    
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    set(handles.poplistPredictedCellClassSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    
    % intialize state of all controls
    set( handles.CellStateAnalyzer, 'Name', sprintf( 'Cell State Analyzer - %s', annotationFile) );
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % close status dialog
    closeStatusDialog(hStatusDialog);

% --------------------------------------------------------------------
function File_SaveAnalysis_Callback(~, eventdata, handles)
% hObject    handle to File_SaveAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % ask user to select the directory in which to save the data    
    if isfield( handles.history, 'lastOutputDir' )
        outputDir = uigetdir( handles.history.lastOutputDir, 'Select output directory');
    else
        outputDir = uigetdir( handles.history.lastAnalyzedNucleusDir, 'Select output directory');
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
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save( fullfile(outputDir, 'CellStateAnalyzer.mat'), '-struct', 'data' );
    closeStatusDialog(h);
    
    clear data;

    % write summary file 
    summary_fid = fopen( fullfile(outputDir, 'cellStateAnalysisSummary.txt'), 'w' );    
    
        % print some information about the dataset
        fprintf( summary_fid, '\n>> Dataset Description:\n' );
        
        fprintf( summary_fid, '\n\tTwo-photon Histone (CFP) data file -- %s\n', handles.data.dataFilePath{1} );
    
        fprintf( summary_fid, '\n\tConfocal FUCCI data file -- %s\n', handles.data.dataFilePath{2} );
        
        fprintf( summary_fid, '\n\tImage Size - [ %s ]\n', sprintf( ' %d ', handles.data.metadata.volSize) );
        fprintf( summary_fid, '\n\tImage Spacing - [ %s ]\n', sprintf( ' %.2f ', handles.data.metadata.voxelSpacing) );
        
        % print analysis summary
        fprintf( summary_fid, '\n>> Analysis Summary:\n' );

        numTotalCells = numel( handles.data.cellStats );
        fprintf( summary_fid, '\n\t%d cells were found by the segmentation algorithm\n', numTotalCells );

        cellPatternStats = { 'Histone (CFP) data', handles.data.dataFilePath{1}, [] };        
        cellPatternStats(end+1,:) = { 'FUCCI data', handles.data.dataFilePath{2}, [] };        
        cellPatternStats(end+1,:) = { 'Image Size', [ '[ ', sprintf( ' %d ', handles.data.metadata.volSize) ,' ]' ], [] };        
        cellPatternStats(end+1,:) = { 'Image Spacing', [ '[ ', sprintf( ' %.2f ', handles.data.metadata.voxelSpacing) ,' ]' ], [] };        
        
        cellPatternStats(end+2,1:2) = { 'Total Cells Detected', numTotalCells };        
        cellPatternStats(end+2,1) = { 'Cell State Distribution' };        
        cellPatternStats(end+1,:) = { 'CellPattern', 'Count', 'Percentage' };
        
        for i = 1:numel(handles.data.cellPatternTypes)+1          

            if i <= numel(handles.data.cellPatternTypes)
                
                curPatternCellIds = find([handles.data.cellStats.cellPatternId] == i);
                numCellsInCurClass = numel( curPatternCellIds );
                fprintf( summary_fid, '\n\t-- %s -- %d (%.2f%%) cells\n', ...
                         handles.data.cellPatternTypes{i}, ...
                         numCellsInCurClass, 100*numCellsInCurClass/numTotalCells );

                % note counts and percentages -- will be written to an excel file
                cellPatternStats{end+1, 1} = handles.data.cellPatternTypes{i}; % pattern type
                cellPatternStats{end, 2} = numCellsInCurClass; % count
                cellPatternStats{end, 3} = numCellsInCurClass * 100.0 / numTotalCells; % percentage
                
            else
                
                % print stats of all cells in this iteration of the for-loop
                curPatternCellIds = 1:numel(handles.data.cellStats); 
                numCellsInCurClass = numel( curPatternCellIds );
                fprintf( summary_fid, '\n\t-- Cell Stats for all %d cells\n', numTotalCells );
                
            end
                
            % get summary info of cells in this class
            if numCellsInCurClass > 2
                
                curClassCellStats = handles.data.cellStats(curPatternCellIds);
                curClassPixelList = ismember( handles.data.imLabelCellSeg, curPatternCellIds );  
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

            end
            
        end    
        
    fclose( summary_fid );
    
    xlswrite( fullfile(outputDir, 'cellCycleStateStats.xlsx'), cellPatternStats );
    
    % ask if the user wants to save annotated cell images
    strQuestion = 'Do you want save images of the cells?';
    button = questdlg( strQuestion, 'Save cell images', 'Yes', 'No', 'No' );
    flagSaveImages = strcmp( button, 'Yes' ); 
    
    h = waitbar(0, 'Saving Images of Cells ... Please Wait' );
    if flagSaveImages

        % create sub-directories for each cell class
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
            curCellOutputDir = fullfile(outputDirImages, curCellPatternType); 
            
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
            imCurCellCropped = double(handles.data.imageData{1}(subinds{1:2}, :));
            imCurCellSegCropped = (handles.data.imLabelCellSeg(subinds{1:2}, :) == cellId);
            
            imCurCellMIP = imresize( mat2gray(max(imCurCellCropped .* imCurCellSegCropped, [], 3)), szOutputImage);

            imwrite( imCurCellMIP, fullfile(curCellOutputDir, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   
                 
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
                     fullfile(curCellOutputDir, sprintf('CellMidSliceHistone_%.3d.png', cellId)), 'png' );   

            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel(:, :, 2:3), channelcolormap(2:3, :) ), ...
                     fullfile(curCellOutputDir, sprintf('CellMidSliceFUCCI_%.3d.png', cellId)), 'png' );   

            imwrite( genMultiChannelOverlay( imCurCellMidSliceAllChannel, channelcolormap ), ...
                     fullfile(curCellOutputDir, sprintf('CellMidSliceAllChannel_%.3d.png', cellId)), 'png' );   

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

    close(handles.CellStateAnalyzer);
        
% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to file_close CellStateAnalyzer.
function CellStateAnalyzer_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CellStateAnalyzer (see GCBO)
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
function File_SetParameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_SetParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function CellStateAnalyzer_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellStateAnalyzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    if isfield(handles, 'history')
        history = handles.history;
        [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
        historyFile = fullfile( pathstr, 'CellStateAnalyzerHistory.mat' );
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

% --- Executes on button press in btnShowPreviousCellInSelectedClass.
function btnShowPreviousCellInSelectedClass_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowPreviousCellInSelectedClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded || numel(handles.data.cellStats) == 0
        return;
    end
    
    curSelPredictedClassId = get(handles.poplistPredictedCellClassSelector, 'Value');
    
    curSelGroundtruthClassId = [];
    if isfield( handles.data, 'groundtruth' )
        curSelGroundtruthClassId = get(handles.poplistGroundtruthCellClassSelector, 'Value');
        if curSelGroundtruthClassId > numel(handles.data.groundtruth.cellPatternTypes)
            curSelGroundtruthClassId = []; % Any type
        end
    end
    
    if isempty( curSelGroundtruthClassId )
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            if handles.dataDisplay.curCellId > 1
                curCellId = handles.dataDisplay.curCellId - 1; % Any predicted type
            else
                return;
            end
        else
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            predictedCellPatternVec = predictedCellPatternVec(handles.dataDisplay.curCellId-1:-1:1);

            indFirstCellInCurClass = find(predictedCellPatternVec == curSelPredictedClassId);

            if isempty( indFirstCellInCurClass )
                return;
            end

            curCellId = handles.dataDisplay.curCellId - indFirstCellInCurClass(1);
        end
        
    else
        
        groundtruthCellPatternVec = [handles.data.groundtruth.cellStats.cellPatternId];
        groundtruthCellPatternVec = groundtruthCellPatternVec(handles.dataDisplay.curCellId-1:-1:1);
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            indFirstCellInCurClass = find(groundtruthCellPatternVec == curSelGroundtruthClassId); %Any predicted type
        else        
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            predictedCellPatternVec = predictedCellPatternVec(handles.dataDisplay.curCellId-1:-1:1);

            indFirstCellInCurClass = find( predictedCellPatternVec == curSelPredictedClassId & ...
                                           groundtruthCellPatternVec == curSelGroundtruthClassId);
        end
        
        if isempty( indFirstCellInCurClass )
            return;
        end
        
        curCellId = handles.dataDisplay.curCellId - indFirstCellInCurClass(1);
        
    end
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
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

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded || numel(handles.data.cellStats) == 0
        return;
    end
    
    curSelPredictedClassId = get(handles.poplistPredictedCellClassSelector, 'Value');
    
    curSelGroundtruthClassId = [];
    if isfield( handles.data, 'groundtruth' )
        curSelGroundtruthClassId = get(handles.poplistGroundtruthCellClassSelector, 'Value');
        if curSelGroundtruthClassId > numel(handles.data.groundtruth.cellPatternTypes)
            curSelGroundtruthClassId = []; % Any type
        end
    end
    
    if isempty( curSelGroundtruthClassId )
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            if handles.dataDisplay.curCellId < numel(handles.data.cellPatternTypes)
                curCellId = handles.dataDisplay.curCellId + 1; % Any predicted type
            else
                return;
            end
        else
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            predictedCellPatternVec = predictedCellPatternVec(handles.dataDisplay.curCellId+1:end);

            indFirstCellInCurClass = find(predictedCellPatternVec == curSelPredictedClassId);

            if isempty( indFirstCellInCurClass )
                return;
            end

            curCellId = handles.dataDisplay.curCellId  + indFirstCellInCurClass(1);
        end
        
    else
        
        groundtruthCellPatternVec = [handles.data.groundtruth.cellStats.cellPatternId];
        groundtruthCellPatternVec = groundtruthCellPatternVec(handles.dataDisplay.curCellId+1:end);
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            indFirstCellInCurClass = find(groundtruthCellPatternVec == curSelGroundtruthClassId); %Any predicted type
        else        
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            predictedCellPatternVec = predictedCellPatternVec(handles.dataDisplay.curCellId+1:end);

            indFirstCellInCurClass = find( predictedCellPatternVec == curSelPredictedClassId & ...
                                           groundtruthCellPatternVec == curSelGroundtruthClassId);
        end
        
        if isempty( indFirstCellInCurClass )
            return;
        end
        
        curCellId = handles.dataDisplay.curCellId  + indFirstCellInCurClass(1);
        
    end
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

% --- Executes on selection change in poplistPredictedCellClassSelector.
function poplistPredictedCellClassSelector_Callback(hObject, eventdata, handles)
% hObject    handle to poplistPredictedCellClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistPredictedCellClassSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistPredictedCellClassSelector

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded || numel(handles.data.cellStats) == 0
        return;
    end
    
    curSelPredictedClassId = get(handles.poplistPredictedCellClassSelector, 'Value');
    
    curSelGroundtruthClassId = [];
    if isfield( handles.data, 'groundtruth' )
        curSelGroundtruthClassId = get(handles.poplistGroundtruthCellClassSelector, 'Value');
        if curSelGroundtruthClassId > numel(handles.data.groundtruth.cellPatternTypes)
            curSelGroundtruthClassId = []; % Any type
        end
    end
    
    if isempty( curSelGroundtruthClassId )
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
             curCellId = 1; % Any predicted type
        else
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];    
            indFirstCellInCurClass = find(predictedCellPatternVec == curSelPredictedClassId);
      
            if isempty( indFirstCellInCurClass )
                return;
            end

            curCellId = indFirstCellInCurClass(1);
        end
        
    else
        
        groundtruthCellPatternVec = [handles.data.groundtruth.cellStats.cellPatternId];
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            indFirstCellInCurClass = find(groundtruthCellPatternVec == curSelGroundtruthClassId); %Any predicted type
        else        
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            indFirstCellInCurClass = find( predictedCellPatternVec == curSelPredictedClassId & ...
                                           groundtruthCellPatternVec == curSelGroundtruthClassId);
        end
        
        if isempty( indFirstCellInCurClass )
            return;
        end
        
        curCellId = indFirstCellInCurClass(1);
        
    end
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
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
function poplistPredictedCellClassSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistPredictedCellClassSelector (see GCBO)
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

    if ~isempty(handles.imarisAppCellSeg) && isobject(handles.imarisAppCellSeg)
        delete( handles.imarisAppCellSeg );
    end       

%     segMaskRGB = label2rgbND(handles.data.imLabelCellSeg);
%     
%     stats = regionprops( bwlabeln( handles.dataDisplay.imCellSeedPoints ), 'Centroid' );
%     cellSeedPointLocations = cat( 1, stats.Centroid );
% 
%     handles.imarisAppCellSeg = Display3DDataAndResultsInImaris( handles.data.imageData{1}, handles.data.metadata.voxelSpacing, ...
%                                                                 'rgbMask', segMaskRGB, ...
%                                                                 'spotLocations', cellSeedPointLocations, ...
%                                                                 'spotRadius', 3);

    % ----
    imageDataPadded = cell(1,3);
    for i = 1:3
        imageDataPadded{i} = padarray( handles.data.imageData{i}, ones(1,3), min(handles.data.imageData{i}(:)) );
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
    stats = regionprops( bwlabeln( handles.dataDisplay.imCellSeedPoints ), 'Centroid' );
    cellSeedPointLocations = cat( 1, stats.Centroid );

    % Display everything in imaris
    handles.imarisAppCellSeg = DisplayMultichannel3DDataInImaris( handles.data.imageData, ...
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
    imCellSeedCropped = padarray( handles.dataDisplay.imCellSeedPoints( subinds{:} ), ones(1,3), 0 );
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

        case { 'Mono_Pre_G1', 'Multi_Pre_G1', 'Mono_G1', 'Multi_G1', 'G1', 'Interphase' } 

            cellColor = [1.0, 0.0, 0.0, 0.5]; % red

        case { 'Mono_S', 'Multi_S', 'S' } 

            cellColor = [1.0, 1.0, 0.0, 0.5]; % yellow

        case { 'Mono_G2', 'Multi_G2', 'G2' } 

           cellColor = [0.0, 1.0, 0.0, 0.5]; % green

        case { 'Mono_Prophase', 'Multi_Prophase', 'M' }     

            cellColor = [0.0, 1.0, 1.0, 0.5]; % orange

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
    
    if ~isempty(handles.imarisAppCellPattern) && isobject(handles.imarisAppCellPattern)
        delete( handles.imarisAppCellPattern );
    end           
    
    imageDataPadded = cell(1,3);
    for i = 1:3
        imageDataPadded{i} = padarray( handles.data.imageData{i}, ones(1,3), min(handles.data.imageData{i}(:)) );
    end
    
    % compute isosurface geometry for cells in each pattern 
    hStatusDlg = waitbar( 0, 'computing surface geometry for cells in each pattern' );
    cellPatternSurfaceObjectList = {};
    cellPatternSurfaceObjectList_Border = {};
    
    borderCellIds = GetBorderCellIds( handles );
    flagIsBorderCell = false( numel(handles.data.cellStats), 1 );
    flagIsBorderCell( borderCellIds ) = true;
    
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
        curPatternIsoSurface_Border.name = [ handles.data.cellPatternTypes{pid}, '_BorderCell' ];
        
        % color
        curPatternIsoSurface.color = MapCellPatternToColor( handles.data.cellPatternTypes{pid} );
        curPatternIsoSurface_Border.color = MapCellPatternToColor( handles.data.cellPatternTypes{pid} );
        
        % compute surface geometry of each cell of current pattern        
        curPatternIsoSurface.surfaces = [];
        curPatternIsoSurface_Border.surfaces = [];
        
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
            if flagIsBorderCell(cid)
                curPatternIsoSurface_Border.surfaces = [ curPatternIsoSurface_Border.surfaces; curCellSurfaceGeometry ];
            else                
                curPatternIsoSurface.surfaces = [ curPatternIsoSurface.surfaces; curCellSurfaceGeometry ];
            end
            
        end
        
        % add cell pattern surface to the list
       cellPatternSurfaceObjectList{end+1} = curPatternIsoSurface;
       cellPatternSurfaceObjectList_Border{end+1} = curPatternIsoSurface_Border;
       
    end    
    closeStatusDialog( hStatusDlg );
    
    cellPatternSurfaceObjectList = cat(2, cellPatternSurfaceObjectList, cellPatternSurfaceObjectList_Border);
    
    % display cell pattern distribution in imaris
    handles.imarisAppCellPattern = DisplayMultichannel3DDataInImaris( imageDataPadded, ...
                                                                      'spacing', handles.data.metadata.voxelSpacing, ...
                                                                      'surfaceObjects', cellPatternSurfaceObjectList, ...
                                                                      'displayRanges', handles.dataDisplay.imDisplayRange, ...
                                                                      'displayColors', handles.data.metadata.channelColors );
    
    % Update handles structure
    guidata(hObject, handles);
    

%---------------------------------------------------------------------
function [borderCellIds] = GetBorderCellIds( handles )

    borderCellIds = setdiff( unique(handles.data.imLabelCellSeg .* ~handles.data.imRegValidMask), 0 );
    
% --------------------------------------------------------------------
function File_Predict_Cell_State_Callback(hObject, eventdata, handles)
% hObject    handle to File_Predict_Cell_State (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~handles.flagDataLoaded
        return;
    end

    PerformCellCycleStateIdentification(hObject, handles);

% --------------------------------------------------------------------
function File_Load_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % get analysis file from user
    if isfield( handles.history, 'lastAnalysisDir' )
        [fileName,pathName] = uigetfile( fullfile( handles.history.lastAnalysisDir, 'CellStateAnalyzer.mat' ), 'Select analysis file' );   
    else
        [fileName,pathName] = uigetfile( 'CellStateAnalyzer.mat', 'Select analysis file' );   
    end

    if ~fileName 
        return;
    end

    hStatusDialog = waitbar(0, 'Loading selected analysis file ...');
    
    handles.history.lastAnalysisDir = pathName;
    analysisFile = fullfile(pathName, fileName);
    
    % load analysis data from file
    analysisData = load( analysisFile );
    
    % extract all the data needed
    handles.data = [];
    
        % basic data
        handles.data.dataFilePath = analysisData.dataFilePath;
        handles.data.metadata = analysisData.metadata;
        handles.data.imageData = analysisData.imageData;
        handles.data.imRegValidMask = analysisData.imRegValidMask;
        
        if ~isfield( handles.data.metadata, 'channelColors' )
            handles.data.metadata.channelColors = [ 0 0 1; 0 1 0; 1 0 0 ];
        end
        handles = ComputeDisplayData(handles);
        
        % segmentation stuff
        handles.data.imLabelCellSeg = analysisData.imLabelCellSeg;        
        handles.data.imCellSeedPoints = analysisData.imCellSeedPoints;        

        if ~isfield( analysisData, 'CellSegColorMap' )
            [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND( handles.data.imLabelCellSeg );
        else
            handles.data.CellSegColorMap = analysisData.CellSegColorMap;        
            handles.dataDisplay.imCellSegRGBMask = label2rgbND( handles.data.imLabelCellSeg, handles.data.CellSegColorMap );
        end
        
        handles.dataDisplay.imCellSeedPoints = imdilate(analysisData.imCellSeedPoints, ones(3,3,3));
        
        % load cell pattern types
        handles.data.cellPatternTypes = analysisData.cellPatternTypes;
        handles.data.cellStats = ComputeCellProperties( handles );
    
        for i = 1:numel(handles.data.cellStats)                
            curCellPatternType = analysisData.cellStats(i).cellPatternType;                
            curCellPatternId = find( strcmpi(handles.data.cellPatternTypes, curCellPatternType) ); 
            handles.data.cellStats(i).cellPatternId = curCellPatternId;  
            handles.data.cellStats(i).cellPatternType = curCellPatternType;
        end
            
    % data ready for display
    handles.flagDataLoaded = true;        

    % set cell pattern type list and cell class selector list
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );
    set( handles.poplistPredictedCellClassSelector, 'String', handles.data.cellPatternTypes );
    set(handles.poplistGroundtruthCellClassSelector, 'Enable', 'off' );               

    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);
    set(handles.poplistPredictedCellClassSelector, 'Value', handles.data.cellStats(handles.dataDisplay.curCellId).cellPatternId);

    % intialize state of all controls
    set( handles.CellStateAnalyzer, 'Name', sprintf( 'Cell State Analyzer - %s', analysisFile ) );

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % close status dialog
    closeStatusDialog(hStatusDialog);    
    
% --- Executes on selection change in poplistGroundtruthCellClassSelector.
function poplistGroundtruthCellClassSelector_Callback(hObject, eventdata, handles)
% hObject    handle to poplistGroundtruthCellClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistGroundtruthCellClassSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistGroundtruthCellClassSelector

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded || numel(handles.data.cellStats) == 0
        return;
    end
    
    if ~isfield( handles.data, 'groundtruth' )
        return;
    end
     
    curSelPredictedClassId = get(handles.poplistPredictedCellClassSelector, 'Value');
    
    curSelGroundtruthClassId = get(handles.poplistGroundtruthCellClassSelector, 'Value');
    if curSelGroundtruthClassId > numel(handles.data.groundtruth.cellPatternTypes)
        curSelGroundtruthClassId = []; % Any groundtruth type
    end
    
    if isempty( curSelGroundtruthClassId )
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
             curCellId = 1; % Any predicted type
        else
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];    
            indFirstCellInCurClass = find(predictedCellPatternVec == curSelPredictedClassId);
      
            if isempty( indFirstCellInCurClass )
                return;
            end

            curCellId = indFirstCellInCurClass(1);
        end
        
    else
        
        groundtruthCellPatternVec = [handles.data.groundtruth.cellStats.cellPatternId];
        
        if curSelPredictedClassId > numel(handles.data.cellPatternTypes) 
            indFirstCellInCurClass = find(groundtruthCellPatternVec == curSelGroundtruthClassId); %Any predicted type
        else        
            predictedCellPatternVec = [handles.data.cellStats.cellPatternId];
            indFirstCellInCurClass = find( predictedCellPatternVec == curSelPredictedClassId & ...
                                           groundtruthCellPatternVec == curSelGroundtruthClassId);
        end
        
        if isempty( indFirstCellInCurClass )
            return;
        end
        
        curCellId = indFirstCellInCurClass(1);
        
    end
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
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
function poplistGroundtruthCellClassSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistGroundtruthCellClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_SetParameters_NucleiSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to File_SetParameters_NucleiSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters.segmentation = NucleiSegmentationParametersGUI( handles.parameters.segmentation );

    % Update handles structure
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function File_SetParametets_CellCycleStateIdentification_Callback(hObject, eventdata, handles)
% hObject    handle to File_SetParametets_CellCycleStateIdentification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters.classification = CellCycleStateIdentificationParametersGUI( handles.parameters.classification );

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function CellStateAnalyzer_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to CellStateAnalyzer (see GCBO)
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
