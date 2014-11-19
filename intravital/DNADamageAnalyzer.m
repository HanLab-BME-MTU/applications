function varargout = DNADamageAnalyzer(varargin)
%DNADAMAGEANALYZER M-file for DNADamageAnalyzer.fig
%      DNADAMAGEANALYZER, by itself, creates a new DNADAMAGEANALYZER or raises the existing
%      singleton*.
%
%      H = DNADAMAGEANALYZER returns the handle to a new DNADAMAGEANALYZER or the handle to
%      the existing singleton*.
%
%      DNADAMAGEANALYZER('Property','Value',...) creates a new DNADAMAGEANALYZER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to DNADamageAnalyzer_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DNADAMAGEANALYZER('CALLBACK') and DNADAMAGEANALYZER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DNADAMAGEANALYZER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DNADamageAnalyzer

% Last Modified by GUIDE v2.5 18-Aug-2014 13:54:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DNADamageAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @DNADamageAnalyzer_OutputFcn, ...
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

% --- Executes just before DNADamageAnalyzer is made visible.
function DNADamageAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

    p = inputParser;
    p.addParamValue('mode', 'analysis', @(x) (ischar(x) && ismember(x, {'training', 'analysis'})) );
    p.addParamValue('flagParallelize', false, @isscalar);
    p.addParamValue('flagDebugMode', false, @isscalar);
    p.parse( varargin{:} );    
    PARAMETERS = p.Results;

    % Choose default command line output for DNADamageAnalyzer
    handles.output = hObject;

    % make sure weka is in the path    
    AddWekaClassesToPath();
    
    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'DNADamageAnalyzerHistory.mat' );
    if exist( historyFile, 'file' )
       handles.history = load( historyFile );
    else
       handles.history.lastImageFileDir = pathstr;
    end

    % Initialize global variables
    handles.cellDisplaySize = 70; 
    
    handles.flagDataLoaded = false;
    handles.flagUseLOG = get(handles.CheckboxGlobalLog, 'Value');
    handles.flagShowNucleiSegMask = get(handles.CheckboxNucleiSegMask, 'Value');
    handles.flagShowNucleiSeedMask = get(handles.CheckboxNucleiSeedMask, 'Value');
    handles.flagShowFociSegMask = get(handles.CheckboxFociSegMask,'Value');
    handles.flagShowCellBBox = get(handles.CheckboxCellBBox, 'Value');
    handles.flagShowDrugMask = get(handles.CheckboxDrugMask, 'Value');
    handles.flagShowMacrophageMask = get(handles.CheckboxMacrophageMask, 'Value');
    
    handles.imarisAppCellSeg = [];
    handles.imarisAppCellPattern = [];
    handles.imarisAppCellSegCropped = [];

    handles.flagParallelize = PARAMETERS.flagParallelize;
    handles.flagDebugMode = PARAMETERS.flagDebugMode;
    
    handles.drugMaskColor = [0, 1, 0];
    handles.macrophageMaskColor = [1, 0, 0];
    
    % default parameters
    f = rdir( fullfile(fileparts(mfilename('fullpath')), '**', 'regionMerging.model') );
    if ~isempty(f)
        handles.defaultParameters.segmentation = NucleiSegmentationParametersGUI( 'defaultWithRegionMerging' );
        handles.defaultParameters.segmentation.flagPerformRegionMerging = true;
        handles.defaultParameters.segmentation.regionMergingModelFile = f(1).name;
    else
        handles.defaultParameters.segmentation = NucleiSegmentationParametersGUI( 'defaultWithoutRegionMerging' );
    end
    
    handles.defaultParameters.fociDetection = FociDetectionParametersGUI('default');
    handles.parameters = handles.defaultParameters;
    
    % open matlab pool for parallel processing    
    handles.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.flagParallelize && ~handles.flagPoolOpenedAlready 
        matlabpool open;
    end
    
    % Set callbacks
    set(handles.DNADamageAnalyzer, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes DNADamageAnalyzer wait for user response (see UIRESUME)
    % uiwait(handles.DNADamageAnalyzer);

% --- Outputs from this function are returned to the command line.
function varargout = DNADamageAnalyzer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

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

    % ask the use to select the oif file for the nucleus channel
    status = bfCheckJavaPath(1);
    [fileName,pathName] = uigetfile(bfGetFileExtensions, ...
                                    'Select the data file', ...
                                    handles.history.lastImageFileDir);
    
    if ~fileName 
        return;
    end
    
    dataFilePath = fullfile( pathName, fileName );
    handles.history.lastImageFileDir = pathName;
    
    % load image data
    PrettyPrintStepDescription( 'Loading Image Data' );
    
        % initialize bio image reader
        try 
            br = BioImageReader(dataFilePath);    
        catch err
            fprintf( 'ERROR: could not load image data data from file %s', dataFilePath );
            err
            errordlg(err.message, 'Error loading image data from file');
            return;
        end
    
        % ask user to select a series
        if br.getSeriesCount() > 1
            [usel, bsel] = settingsdlg('Description', 'Select series', ...
                                       'seriesId', num2cell(1:br.getSeriesCount()));
            if ~strcmpi(bsel, 'ok')
                return;
            end
            seriesId = usel.seriesId;
        else
            seriesId = 1;
        end
        
        % get series metadata
        metadata = br.getMetadata(seriesId);
        
        % make sure that the dataset is 3D
        if numel(metadata.imageSize) ~= 3 || metadata.imageSize(3) < 2
            errordlg( 'Dataset is not 3D' );
            return;
        end
        
        % ask user to set the image properties
        switch metadata.numChannels

            case 1

                [imageData, metadata] = uiGetBioImageData(dataFilePath, 'seriesId', seriesId);

                if isempty(imageData)
                    return;
                end
                
                metadata.channelId53BP1 = 1;
                metadata.channelIdDrug = [];
                metadata.channelIdMacrophage = [];
                metadata.channelNames = {'53BP1'};
                
            case 2

                [imageData, metadata] = uiGetBioImageData(dataFilePath, ...
                                                          'seriesId', seriesId, ...  
                                                          'channelDescriptionsAndNamesList', ...  
                                                          {{'53BP1 Channel', 'channelId53BP1'}, ...
                                                           {'Macrophage Channel', 'channelIdMacrophage'}});
                
                if isempty(imageData)
                    return;
                end
                
                if numel( unique([metadata.channelId53BP1, metadata.channelIdMacrophage]) ) ~= 3
                    errordlg('Drug, 53BP1, and Macrophage channels should be distinct');
                    return;
                end

                metadata.channelIdDrug = [];
                metadata.channelNames{metadata.channelId53BP1} = '53BP1';
                metadata.channelNames{metadata.channelIdMacrophage} = 'Macrophage';
                
            case 3

                [imageData, metadata] = uiGetBioImageData(dataFilePath, ...
                                                          'seriesId', seriesId, ...  
                                                          'channelDescriptionsAndNamesList', ...  
                                                          {{'Drug Channel', 'channelIdDrug'}, ...
                                                           {'53BP1 Channel', 'channelId53BP1'}, ...
                                                           {'Macrophage Channel', 'channelIdMacrophage'}});
                
                if isempty(imageData)
                    return;
                end
                
                if numel( unique([metadata.channelId53BP1, metadata.channelIdDrug, metadata.channelIdMacrophage]) ) ~= 3
                    errordlg('Drug, 53BP1, and Macrophage channels should be distinct');
                    return;
                end

                metadata.channelNames{metadata.channelId53BP1} = '53BP1';
                metadata.channelNames{metadata.channelIdDrug} = 'Drug';
                metadata.channelNames{metadata.channelIdMacrophage} = 'Macrophage';

            otherwise

                errordlg('Invalid Dataset selected. It must contain three channels - Drug, 53BP1, Macrophages');
                return;
        end        
        
    % store image data in handles structures    
    clear handles.data;
    
    handles.flagDataLoaded = true;
    handles.metadata.dataFilePath = dataFilePath;
    handles.data.metadata = metadata;
    handles.data.imageData = imageData;
    
    % update controls    
    set( handles.DNADamageAnalyzer, 'Name', sprintf( 'DNA Damage Analysis - %s', handles.metadata.dataFilePath ) );
    set( handles.popupChannelSelector, 'String', handles.data.metadata.channelNames);
    
    % pre-compute data needed for display
    handles = ComputeDisplayData( handles );
    
    % Run analysis
    handles = RunAnalysis(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % Update Foci Count Distribution Plot
    UpdateFociCountDistributionPlot(handles);

% --------------------------------------------------------------------    
function [ handles ] = ComputeDisplayData(handles)

    for i = 1:numel( handles.data.imageData )
        handles.data.imageData{i} = double( handles.data.imageData{i} );
        handles.dataDisplay.imDisplayRange(i,:) = ComputeImageDynamicRange( handles.data.imageData{i}, 98.0 );   
        handles.dataDisplay.imageDataLOG{i} = ComputeImageLogTransformForDisplay( handles.data.imageData{i} );
        handles.dataDisplay.imLogDisplayRange(i,:) = [ min(handles.dataDisplay.imageDataLOG{i}(:)), max(handles.dataDisplay.imageDataLOG{i}(:))];   
    end
    
% --------------------------------------------------------------------    
function handles = RunAnalysis(hObject, handles)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    % segment cells
    PrettyPrintStepDescription( 'Running Nuclei Sementation Algorithm' );    
    hStatusDialog = waitbar(0, 'Segmentating Nuclei');
   
    regionMergingModelFile = [];
    if handles.parameters.segmentation.flagPerformRegionMerging
        regionMergingModelFile = handles.parameters.segmentation.regionMergingModelFile;
    end
    
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints, ...
     handles.data.segAlgoParameters ] = segmentCellsInIntravitalData( handles.data.imageData{handles.data.metadata.channelId53BP1}, ...
                                                                      handles.data.metadata.pixelSize, ...                                                                      
                                                                      'flagParallelize', handles.flagParallelize, ...
                                                                      'flagDebugMode', handles.flagDebugMode, ...
                                                                      'cellDiameterRange', handles.parameters.segmentation.cellDiameterRange, ...
                                                                      'thresholdingAlgorithm', 'BackgroudRemovalUsingMorphologicalOpening', ...
                                                                      'minSignalToBackgroundRatio', 2.5, ...
                                                                      'seedPointDetectionAlgorithm', handles.parameters.segmentation.seedPointDetectionAlgorithm, ...
                                                                      'minCellVolume', handles.parameters.segmentation.minCellVolume, ...
                                                                      'flagIgnoreCellsOnXYBorder', handles.parameters.segmentation.flagIgnoreXYBorderCells, ...
                                                                      'regionMergingModelFile', regionMergingModelFile);

    [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.dataDisplay.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints, ones(3,3,3) );
    
    % compute properties of each cell
    handles.data.cellStats = ComputeCellProperties( handles );

    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    handles.dataDisplay.curChannelId = handles.data.metadata.channelId53BP1;

    % update ui controls
    set(handles.popupChannelSelector, 'Value', handles.dataDisplay.curChannelId);
    
    % Detect foci
    waitbar(1/3, hStatusDialog, 'Detecting puncta');
    handles = PerformFociSegmenatation(hObject, handles);
    
    % measure colocalization with other channels
    waitbar(1/3, hStatusDialog, 'Computing Colocalization Measures');
    
    handles = PerformColocalizationAnalysis(hObject, handles);
    
    % close progress bar
    closeStatusDialog(hStatusDialog);

% -----------------------------------------------------------------
function handles = PerformFociSegmenatation(hObject, handles)

    if ~handles.flagDataLoaded
        return;
    end

    PrettyPrintStepDescription( 'Detecting Foci' );

    modelFile = fullfile('Z:\intravital\DNADamageAnalysis\models\fociDetection', ...
                         'train_M16_pre_3to5h_24h', 'fociDetectionFeatures.model');
                     
    [handles.data.fociStats, ...
     handles.data.imFociSeedPoints, ...
     handles.data.imLabelFociSeg, ...
     handles.data.fociDetectionParameters ] = segmentFociInsideNuclei(handles.data.imageData{handles.data.metadata.channelId53BP1}, ...
                                                                      handles.parameters.fociDetection.fociDiameterRange, ...
                                                                      'spacing', handles.data.metadata.pixelSize, ...
                                                                      'roiMask', handles.data.imLabelCellSeg, ...
                                                                      'minSignalToBackgroundRatio', handles.parameters.fociDetection.minSignalToBackgroundRatio, ...
                                                                      'minDogResponse', handles.parameters.fociDetection.minDoGResponse, ...
                                                                      'minHessianEigenRatio21', handles.parameters.fociDetection.minHessianEigenRatio21, ...
                                                                      'minHessianEigenRatio31', handles.parameters.fociDetection.minHessianEigenRatio31, ...
                                                                      'minDistanceToROIBoundary', handles.parameters.fociDetection.minDistanceToROIBoundary, ...
                                                                      'fociDetectionModelFile', handles.parameters.fociDetection.fociDetectionModelFile);
                                                                    
   [handles.dataDisplay.imFociSegRGBMask, handles.data.FociSegColorMap] = label2rgbND( handles.data.imLabelFociSeg );
    handles.dataDisplay.imFociSeedPoints = imdilate(handles.data.imFociSeedPoints, ones(3,3,3));
    
    for cid = 1:numel(handles.data.cellStats)
        
        curCellPixInd = find(handles.data.imLabelCellSeg == cid);
        curCellFoci = unique(handles.data.imFociSeedPoints(curCellPixInd));
        curCellFoci = curCellFoci(curCellFoci > 0);
        handles.data.cellStats(cid).foci = curCellFoci;
        handles.data.cellStats(cid).fociCount = numel(curCellFoci);
        
    end

    % update foci selector list
    set(handles.poplistFociCountSelector, 'String', num2cell(unique([handles.data.cellStats.fociCount])) );
    set(handles.poplistFociCountSelector, 'value', 1);
    
    % Update handles structure
    guidata(hObject, handles);
   
% -----------------------------------------------------------------
function handles = PerformColocalizationAnalysis(hObject, handles)

    if ~handles.flagDataLoaded
        return;
    end

    if handles.data.metadata.numChannels ~= 3
        return;
    end
    
    PrettyPrintStepDescription( 'Computing Colocalization Measures' );

    handles.data.cellColocStats = [];
    
    colocTime = tic;
        
    % colocalization with the macrophage channel
    if ~isempty( handles.data.metadata.channelIdMacrophage )

        % segment macrophage channel
        fprintf('\n>> Segmenting macrophage channel ... \n');

        macrophageSegTime = tic;

        imMacrophage = handles.data.imageData{handles.data.metadata.channelIdMacrophage};
        handles.data.imMacrophageSeg = segmentMacrophages(imMacrophage);
        handles.dataDisplay.imMacrophageSegRGB = convertBinaryMaskToRGBMask(handles.data.imMacrophageSeg, handles.macrophageMaskColor);

        fprintf('\ntook %f seconds\n', toc(macrophageSegTime));

        % compute colocalization measures 
        fprintf('\n>> Computing colocalization measures ... \n');

        macrophageColocTime = tic;
        
        [cellColocStats] = ComputeDNADamageColocalizationMeasures(handles.data.imageData{handles.data.metadata.channelId53BP1}, ...   
                                                                  handles.data.imageData{handles.data.metadata.channelIdMacrophage}, ...
                                                                  handles.data.metadata.pixelSize, ...
                                                                  handles.data.imLabelCellSeg, handles.data.cellStats, ...
                                                                  handles.data.imMacrophageSeg);

        handles.data.cellColocStats.macrophage = cellColocStats;

        fprintf('\ntook %f seconds\n', toc(macrophageColocTime));
        
    end

    % colocalization with the drug channel
    if ~isempty( handles.data.metadata.channelIdDrug )
        
        % segment drug channel
        fprintf('\n>> Segmenting drug channel ... \n');

        drugSegTime = tic;

        imDrug = handles.data.imageData{handles.data.metadata.channelIdDrug};
        handles.data.imDrugSeg = segmentDrugCisplatin(imDrug);
        handles.dataDisplay.imDrugSegRGB = convertBinaryMaskToRGBMask(handles.data.imDrugSeg, handles.drugMaskColor);

        fprintf('\ntook %f seconds\n', toc(drugSegTime));

        % compute colocalization measures with drug channel
        fprintf('\n>> Computing colocalization measures for the drug channel ... \n');

        drugColocTime = tic;
        
        [cellColocStats] = ComputeDNADamageColocalizationMeasures(handles.data.imageData{handles.data.metadata.channelId53BP1}, ...   
                                                                  handles.data.imageData{handles.data.metadata.channelIdDrug}, ... 
                                                                  handles.data.metadata.pixelSize, ...
                                                                  handles.data.imLabelCellSeg, handles.data.cellStats, ...
                                                                  handles.data.imDrugSeg);

        handles.data.cellColocStats.drug = cellColocStats;
        
        fprintf('\ntook %f seconds\n', toc(drugColocTime));
        
    end
    
    fprintf('\ncolocalization analysis took a total of %f seconds\n', toc(colocTime));
    
    % Update handles structure
    guidata(hObject, handles);

% --------------------------------------------------------------------
function [cellStats] = ComputeCellProperties( handles )

    cellStats = regionprops( handles.data.imLabelCellSeg, ...
                             'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
    
    for i = 1:numel(cellStats)

        % intensity descriptors
        cellPixelIntensities = handles.data.imageData{handles.data.metadata.channelId53BP1}( cellStats(i).PixelIdxList );
        cellStats(i).meanIntensity = mean( cellPixelIntensities );
        cellStats(i).stdIntensity = std( double(cellPixelIntensities) );
        cellStats(i).minIntensity = min( cellPixelIntensities );
        cellStats(i).maxIntensity = max( cellPixelIntensities );
        
        % volume
        cellStats(i).AreaPhysp = cellStats(i).Area * prod( handles.data.metadata.pixelSize );
        
    end
    
% --------------------------------------------------------------------
function File_Run_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Run_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    handles = RunAnalysis(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % Update Foci Count Distribution Plot
    UpdateFociCountDistributionPlot(handles);

% --------------------------------------------------------------------
function UpdateCellDisplay(handles)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellSliceId = handles.dataDisplay.curCellSliceId;
    curCellCentroid = round( curCellStats.Centroid );    
    curChannelId = handles.dataDisplay.curChannelId;
    
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

    % display mask overlays on global cross section images
    imGlobalXYMasks = {};
    imGlobalXZMasks = {};
    imGlobalYZMasks = {};
    maskAlphas = [];

    nucleiMaskAlpha = 0.2;
    fociMaskAlpha = 0.6;
    seedMaskAlpha = 0.6;
    
    if handles.flagShowNucleiSegMask

        imGlobalXYMasks{end+1} = squeeze( handles.dataDisplay.imCellSegRGBMask( :, :, curCellSliceId(3), : ) );

        imGlobalXZSegMaskRGB = squeeze( handles.dataDisplay.imCellSegRGBMask(curCellSliceId(1), :, :, :) );
        imGlobalXZMasks{end+1} = permute(imGlobalXZSegMaskRGB, [2,1,3] );

        imGlobalYZMasks{end+1} = squeeze( handles.dataDisplay.imCellSegRGBMask(:, curCellSliceId(2), :, :) );

        maskAlphas(end+1) = nucleiMaskAlpha;

    end

    if handles.flagShowFociSegMask && isfield(handles.dataDisplay, 'imFociSegRGBMask')

        imGlobalXYMasks{end+1} = squeeze( handles.dataDisplay.imFociSegRGBMask( :, :, curCellSliceId(3), : ) );

        imGlobalXZSegMaskRGB = squeeze( handles.dataDisplay.imFociSegRGBMask(curCellSliceId(1), :, :, :) );
        imGlobalXZMasks{end+1} = permute(imGlobalXZSegMaskRGB, [2,1,3] );

        imGlobalYZMasks{end+1} = squeeze( handles.dataDisplay.imFociSegRGBMask(:, curCellSliceId(2), :, :) );

        maskAlphas(end+1) = fociMaskAlpha;

    end
    
    if ~isempty(handles.data.metadata.channelIdDrug) && handles.flagShowDrugMask && isfield(handles.dataDisplay, 'imDrugSegRGB')
        
        imGlobalXYMasks{end+1} = squeeze(handles.dataDisplay.imDrugSegRGB(:, :, curCellSliceId(3), :));

        imGlobalXZSegMaskRGB = squeeze(handles.dataDisplay.imDrugSegRGB(curCellSliceId(1), :, :, :));
        imGlobalXZMasks{end+1} = permute(imGlobalXZSegMaskRGB, [2,1,3] );

        imGlobalYZMasks{end+1} = squeeze(handles.dataDisplay.imDrugSegRGB(:, curCellSliceId(2), :, :));

        maskAlphas(end+1) = 0.5;
        
    end

    if ~isempty(handles.data.metadata.channelIdMacrophage) && handles.flagShowMacrophageMask && isfield(handles.dataDisplay, 'imMacrophageSegRGB')
        
        imGlobalXYMasks{end+1} = squeeze(handles.dataDisplay.imMacrophageSegRGB(:, :, curCellSliceId(3), :));

        imGlobalXZSegMaskRGB = squeeze(handles.dataDisplay.imMacrophageSegRGB(curCellSliceId(1), :, :, :));
        imGlobalXZMasks{end+1} = permute(imGlobalXZSegMaskRGB, [2,1,3] );

        imGlobalYZMasks{end+1} = squeeze(handles.dataDisplay.imMacrophageSegRGB(:, curCellSliceId(2), :, :));

        maskAlphas(end+1) = 0.5;
        
    end
    
    if ~isempty(imGlobalXYMasks)
        
        imGlobalXYDisplay = genImageRGBMaskOverlay( imGlobalXY, imGlobalXYMasks, maskAlphas );
        imGlobalXZDisplay = genImageRGBMaskOverlay( imGlobalXZ, imGlobalXZMasks, maskAlphas );
        imGlobalYZDisplay = genImageRGBMaskOverlay( imGlobalYZ, imGlobalYZMasks, maskAlphas );
        
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
        
        imsize = size( handles.data.imageData{handles.data.metadata.channelId53BP1} );
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
    imsize = size(handles.data.imageData{handles.data.metadata.channelId53BP1});
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
        imCellCropped = mat2gray( handles.dataDisplay.imageDataLOG{curChannelId}(subinds{:}), handles.dataDisplay.imLogDisplayRange(1,:) );        
    else
        imCellCropped = mat2gray( handles.data.imageData{curChannelId}(subinds{:}), handles.dataDisplay.imDisplayRange(1,:) );        
    end
    
    imLabelCellSegCropped = handles.data.imLabelCellSeg( subinds{:} );    
    imCellSegCropped = double( imLabelCellSegCropped == handles.dataDisplay.curCellId );    
    imLabelCellSegCropped(~imCellSegCropped) = 0;
    
    imCellSeedCropped = handles.dataDisplay.imCellSeedPoints( subinds{:} );
    imCellSeedCropped(~imCellSegCropped) = 0;
    imCellSeedCropped(:,:,2:3) = 0; 

    % display the nuclei channel of the extracted cell image
    if handles.flagShowNucleiSegMask || handles.flagShowFociSegMask
        
        imCroppedCellRGBMasks = {};
        maskAlphas = [];
        
        if handles.flagShowNucleiSegMask
           imCroppedCellRGBMasks{end+1} = label2rgbND(imLabelCellSegCropped, handles.data.CellSegColorMap); 
           maskAlphas(end+1) = nucleiMaskAlpha;
        end

        if handles.flagShowFociSegMask && isfield(handles.dataDisplay, 'imFociSegRGBMask')

           imLabelFociSegRGBMaskCropped = squeeze( handles.dataDisplay.imFociSegRGBMask( subinds{:}, : ) );
           imLabelFociSegRGBMaskCropped( ~repmat(imCellSegCropped, [1,1,3]) ) = 0;
           imCroppedCellRGBMasks{end+1} = imLabelFociSegRGBMaskCropped; 
           maskAlphas(end+1) = fociMaskAlpha;
           
        end

        if handles.flagShowNucleiSeedMask
           imCroppedCellRGBMasks{end+1} = imCellSeedCropped; 
           maskAlphas(end+1) = seedMaskAlpha;
        end
        
        imNucleusXYDisplay = genImageRGBMaskOverlay(imCellCropped, imCroppedCellRGBMasks, maskAlphas );
        
    else
        
        imNucleusXYDisplay = repmat( imCellCropped, [1,1,3] );
        
    end

    cla( handles.Axes_Nucleus_XY, 'reset' );
    image( imNucleusXYDisplay, 'Parent', handles.Axes_Nucleus_XY );
    set(  handles.Axes_Nucleus_XY, 'XTickLabel', [], 'YTickLabel', [] ); 
    
    % display Nucleus MIP
    subindsMIP = subinds;
    subindsMIP{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));

    imCellSegCropped = handles.data.imLabelCellSeg( subindsMIP{:} ) == handles.dataDisplay.curCellId;
    imCurCellMIP = mat2gray( max(handles.data.imageData{curChannelId}(subindsMIP{:}) .* imCellSegCropped, [], 3 ) );
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
    strCellDescription = sprintf( '%s\n\nNumber of Puncta: %d', strCellDescription, numel(curCellStats.foci) );
    
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
    
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 99.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
% --------------------------------------------------------------------
function File_Load_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % get annotation file from user
    if isfield( handles.history, 'lastAnalysisFileDir' )
        [fileName,pathName] = uigetfile( fullfile(handles.history.lastAnalysisFileDir, 'DNADamageAnalysis.mat'), 'Select analysis file' );   
    else
        [fileName,pathName] = uigetfile( fullfile(pwd, 'DNADamageAnalysis.mat'), 'Select analysis file' );   
    end

    if ~fileName 
        return;
    end
    
    hStatusDialog = waitbar(0, 'Loading selected analysis file ...');
    
    % load annotation data from file
    analysisFile = fullfile(pathName, fileName);
    analysisData = load( analysisFile );
    handles.history.lastAnalysisFileDir = pathName;
    
    % retrieve data needed for this tool
    handles.data = [];
    
        % basic data
        if iscell(analysisData.dataFilePath)
            handles.metadata.dataFilePath = analysisData.dataFilePath{1};
        else
            handles.metadata.dataFilePath = analysisData.dataFilePath;
        end
        handles.data.metadata = analysisData.metadata;
        handles.data.imageData = analysisData.imageData;
        for i = 1:numel(handles.data.imageData)
            handles.data.imageData{i} = double(handles.data.imageData{i});
        end
        
        % basic display data
        handles = ComputeDisplayData(handles);
    
        % nuclei segmentation stuff
        handles.data.imLabelCellSeg = double( analysisData.imLabelCellSeg );        
        handles.data.imCellSeedPoints = double( analysisData.imCellSeedPoints );        

        if ~isfield( analysisData, 'CellSegColorMap' )
            [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND( handles.data.imLabelCellSeg );
        else
            handles.data.CellSegColorMap = analysisData.CellSegColorMap;        
            handles.dataDisplay.imCellSegRGBMask = label2rgbND( handles.data.imLabelCellSeg, handles.data.CellSegColorMap );
        end
        
        handles.dataDisplay.imCellSeedPoints = imdilate(double(handles.data.imCellSeedPoints), ones(3,3,3));

        if isfield(analysisData, 'segAlgoParameters')
            handles.data.segAlgoParameters = analysisData.segAlgoParameters;
        end
        
        % cell stats
        handles.data.cellStats = ComputeCellProperties( handles );
        
        for i = 1:numel(handles.data.cellStats)
            handles.data.cellStats(i).foci = analysisData.cellStats(i).foci;
            handles.data.cellStats(i).fociCount = numel( handles.data.cellStats(i).foci );
        end
        
        % foci 
        handles.data.imLabelFociSeg = double( analysisData.imLabelFociSeg ); 
        handles.data.imFociSeedPoints = double( analysisData.imFociSeedPoints );
        handles.data.fociStats = analysisData.fociStats;

        handles.dataDisplay.imFociSeedPoints = imdilate(handles.data.imFociSeedPoints, ones(3,3,3));
        if ~isfield( analysisData, 'FociSegColorMap' )
            [handles.dataDisplay.imFociSegRGBMask, handles.data.FociSegColorMap] = label2rgbND( handles.data.imLabelFociSeg );
        else
            handles.data.FociSegColorMap = analysisData.FociSegColorMap;        
            handles.dataDisplay.imFociSegRGBMask = label2rgbND( handles.data.imLabelFociSeg, handles.data.FociSegColorMap );
        end
        
        if isfield(analysisData, 'fociDetectionParameters')
            handles.data.fociDetectionParameters = analysisData.fociDetectionParameters;
        end
        
        % drug 
        if handles.data.metadata.numChannels == 3
           
            handles.data.imDrugSeg = analysisData.imDrugSeg;
            handles.dataDisplay.imDrugSegRGB = convertBinaryMaskToRGBMask(handles.data.imDrugSeg, handles.drugMaskColor);
            
            handles.data.imMacrophageSeg = analysisData.imMacrophageSeg;
            handles.dataDisplay.imMacrophageSegRGB = convertBinaryMaskToRGBMask(handles.data.imMacrophageSeg, handles.macrophageMaskColor);
           
            handles.data.cellColocStats = analysisData.cellColocStats;
            
        end
        
    % data ready for display
    handles.flagDataLoaded = true;
    
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    handles.dataDisplay.curChannelId = handles.data.metadata.channelId53BP1;
    
    % update ui controls
    
        % change window name
        set( handles.DNADamageAnalyzer, 'Name', sprintf( 'DNA Damage Analysis - %s', analysisFile ) );

        % update foci selector list
        set( handles.poplistFociCountSelector, 'String', num2cell(unique([handles.data.cellStats.fociCount])) );
        set(handles.poplistFociCountSelector, 'value', 1);

        % update channel selector popup
        set( handles.popupChannelSelector, 'String', handles.data.metadata.channelNames);            
        set(handles.popupChannelSelector, 'Value', handles.dataDisplay.curChannelId);        
        
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % Update Foci Count Distribution Plot
    UpdateFociCountDistributionPlot(handles);
    
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
    if isfield( handles.history, 'lastAnalysisFileDir' )
        outputDir = uigetdir(handles.history.lastAnalysisFileDir, 'Select annotation output directory');
    else
        outputDir = uigetdir(pwd, 'Select annotation output directory');
    end

    if ~outputDir
        return;
    end

    handles.history.lastAnalysisFileDir = outputDir;
    
    % save data
    [pathstr, name, ext] = fileparts( handles.metadata.dataFilePath );
    outputDir = fullfile(outputDir, name);

    if ~isdir( outputDir )
        mkdir( outputDir );
    end
    
    data = handles.data;
    for i = 1:numel(handles.data.imageData)
        data.imageData{i} = uint16(data.imageData{i});
    end
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save( fullfile(outputDir, 'DNADamageAnalysis.mat'), '-struct', 'data' );
    close(h);
    
    clear data;
    
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


% --- Executes when user attempts to file_close DNADamageAnalyzer.
function DNADamageAnalyzer_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnalyzer (see GCBO)
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

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes on button press in CheckboxNucleiSegMask.
function CheckboxNucleiSegMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxNucleiSegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxNucleiSegMask

    % change mask use mode
    handles.flagShowNucleiSegMask = get(hObject,'Value');
    
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

% --- Executes during object deletion, before destroying properties.
function DNADamageAnalyzer_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnalyzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    if isfield(handles, 'history')
        history = handles.history;
        [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
        historyFile = fullfile( pathstr, 'DNADamageAnalyzerHistory.mat' );
        save( historyFile, '-struct', 'history' );  
    end

    % close matlab pool
    if ~handles.flagPoolOpenedAlready && matlabpool( 'size' ) > 0
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
    
    imsize = size(handles.data.imageData{handles.data.metadata.channelId53BP1});

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
function View_Cell_Analysis_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to View_Cell_Analysis_In_Imaris (see GCBO)
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
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );
    
    % create crop indices
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{handles.data.metadata.channelId53BP1});
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
        
        % get foci
        if curCellStats.fociCount > 0
            
            curCellFociStats = handles.data.fociStats(curCellStats.foci);
            fociSeedPointLocations = ind2submat(size(handles.data.imFociSeedPoints), [curCellFociStats.PixelLocationIndex]);
            
            for i = 1:3
               fociSeedPointLocations(:,i) = fociSeedPointLocations(:,i) - min(subinds{i}) + 1;
            end
            
            fociSeedPointLocations = fociSeedPointLocations(:, [2,1,3]); % format as [x, y, z]
            
            imvis.AddSpots(fociSeedPointLocations, 0, ...
                           'radii', [curCellFociStats.Radius], ... 
                           'name', 'Puncta', 'color', [0, 1, 0]);
                       
        end    
    
    % Update handles structure
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function View_Full_Analysis_In_Imaris_Callback(hObject, eventdata, handles)
% hObject    handle to View_Full_Analysis_In_Imaris (see GCBO)
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
        hStatusDlg = waitbar( 0, 'computing surface geometry for cells in each pattern' );
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

        % get foci
        if ~isempty(handles.data.fociStats)
            
            fociSeedPointLocations = ind2submat(size(handles.data.imFociSeedPoints), [handles.data.fociStats.PixelLocationIndex]);
            fociSeedPointLocations = fociSeedPointLocations(:, [2, 1, 3]); % format as [x, y, z]
            
            imvis.AddSpots(fociSeedPointLocations, 0, ...
                           'radii', [handles.data.fociStats.Radius], ... 
                           'name', 'Puncta', 'color', [0, 1, 0]);
                       
        end    
        
        % drug mask
        if isfield(handles.data, 'imDrugSeg')

            drugGeometry = ImarisDataVisualizer.generateSurfaceFromMask(handles.data.imDrugSeg, 'surfaceQuality', surfaceQuality);
            imvis.AddSurfaces(drugGeometry, ...
                              'name', 'Drug', ...
                              'color', [0, 1, 0]);
            
        end
        
        % macrophage mask
        if isfield(handles.data, 'imMacrophageSeg')
            
            macrophageGeometry = ImarisDataVisualizer.generateSurfaceFromMask(handles.data.imMacrophageSeg, 'surfaceQuality', surfaceQuality);
            imvis.AddSurfaces(macrophageGeometry, ...
                              'name', 'Macrophages', ...
                              'color', [1, 0, 0]);
            
        end
        
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
function DNADamageAnalyzer_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnalyzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{handles.data.metadata.channelId53BP1});
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

    end
        
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);    
    

% --- Executes on button press in CheckboxFociSegMask.
function CheckboxFociSegMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxFociSegMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxFociSegMask

    % change mask use mode
    handles.flagShowFociSegMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --------------------------------------------------------------------
function File_Detect_Foci_Callback(hObject, eventdata, handles)
% hObject    handle to File_Detect_Foci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~handles.flagDataLoaded
        return;
    end

    handles = PerformFociSegmenatation(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
    % Update Foci Count Distribution Plot
    UpdateFociCountDistributionPlot(handles);
    
% --------------------------------------------------------------------
function View_Puncta_Distribution_Callback(hObject, eventdata, handles)
% hObject    handle to View_Puncta_Distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    numCells = numel(handles.data.cellStats); 
    punctaCounts = [handles.data.cellStats.fociCount];
    
    figure();
    [counts, centers] = optimalHistogram(punctaCounts, 'discrete');
    bar(centers, counts/sum(counts), ...
        'linewidth', 2.0, ...
        'barwidth', 0.9, ...
        'EdgeColor', 'k', ...
        'FaceColor', 0.7 * ones(1,3));

    [~, fname, ~] = fileparts(handles.metadata.dataFilePath);
    fname = strrep(fname,'_','\_');
    text(0.5, 0.95, {fname, sprintf('#cells = %d, #puncta= %d', numCells, sum(punctaCounts))}, ...
         'units', 'normalized', 'horizontalalignment', 'center');
    
    title({'Distribution of Number of Puncta In a Cell'}, ...
          'FontSize', 12.0, 'FontWeight', 'bold');
    xlabel('Number of puncta per cell');
    ylabel('Proportion of cells');
    ylim([0, 1]);
    
    grid on;
    
    
% --------------------------------------------------------------------
function UpdateFociCountDistributionPlot(handles)
    
    if ~handles.flagDataLoaded
        return;
    end
    
    prev_axis = gca;
    
    % show histograms of individual channels    
    axes( handles.Axes_FociCount_Distribution );
    cla reset;    

    numCells = numel(handles.data.cellStats); 
    punctaCount = zeros(1, numCells);
    for i = 1:numCells
        punctaCount(i) = numel(handles.data.cellStats(i).foci);
    end

    nhist( punctaCount, 'nolegend', 'noerror', 'pdf', 'smooth', 'int', ...
           'xlabel', 'Number of puncta per cell', ...
           'ylabel', 'Proportion of cells', ...
           'minx', 0, 'maxx', 4, ...
           'fsize', get(handles.Axes_FociCount_Distribution, 'FontSize') );
       
    title('Distribution of Number of Puncta In a Cell', 'FontSize', 12.0);
    grid on;
    ylim([0, 1]);
    
    axes(prev_axis);


% --- Executes on button press in btnShowPreviousCellWithSpecifiedFociCount.
function btnShowPreviousCellWithSpecifiedFociCount_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowPreviousCellWithSpecifiedFociCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    fociCountList = cellstr(get(handles.poplistFociCountSelector,'String'));    
    curSelFociCount = str2double(fociCountList{get(handles.poplistFociCountSelector, 'Value')});
    
    cellFociCounts = [handles.data.cellStats.fociCount];    
    cellFociCounts = cellFociCounts(handles.dataDisplay.curCellId-1:-1:1);

    indFirstCellInCurClass = find(cellFociCounts == curSelFociCount);

    if isempty( indFirstCellInCurClass )
        return;
    end

    curCellId = handles.dataDisplay.curCellId  - indFirstCellInCurClass(1);
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

% --- Executes on button press in btnShowNextCellWithSpecifiedFociCount.
function btnShowNextCellWithSpecifiedFociCount_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowNextCellWithSpecifiedFociCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    fociCountList = cellstr(get(handles.poplistFociCountSelector,'String'));    
    curSelFociCount = str2double(fociCountList{get(handles.poplistFociCountSelector, 'Value')});
    
    cellFociCounts = [handles.data.cellStats.fociCount];    
    cellFociCounts = cellFociCounts(handles.dataDisplay.curCellId+1:end);

    indFirstCellInCurClass = find(cellFociCounts == curSelFociCount);

    if isempty( indFirstCellInCurClass )
        return;
    end

    curCellId = handles.dataDisplay.curCellId  + indFirstCellInCurClass(1);
    
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);

% --- Executes on selection change in poplistFociCountSelector.
function poplistFociCountSelector_Callback(hObject, eventdata, handles)
% hObject    handle to poplistFociCountSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistFociCountSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistFociCountSelector

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    fociCountList = cellstr(get(handles.poplistFociCountSelector,'String'));    
    curSelFociCount = str2double(fociCountList{get(handles.poplistFociCountSelector, 'Value')});
    
    cellFociCounts = [handles.data.cellStats.fociCount];    
    indFirstCellInCurClass = find(cellFociCounts == curSelFociCount);

    curCellId = indFirstCellInCurClass(1);
        
    % set current cell id
    handles.dataDisplay.curCellId = curCellId;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function poplistFociCountSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistFociCountSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckboxNucleiSeedMask.
function CheckboxNucleiSeedMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxNucleiSeedMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxNucleiSeedMask

    % change mask use mode
    handles.flagShowNucleiSeedMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);


% --------------------------------------------------------------------
function File_SetParameters_NucleiSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to File_SetParameters_NucleiSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters.segmentation = NucleiSegmentationParametersGUI( handles.parameters.segmentation );
    
    % Update handles structure
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function File_SetParameters_FociDetection_Callback(hObject, eventdata, handles)
% hObject    handle to File_SetParameters_FociDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters.fociDetection = FociDetectionParametersGUI( handles.parameters.fociDetection );
    
    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes on selection change in popupChannelSelector.
function popupChannelSelector_Callback(hObject, eventdata, handles)
% hObject    handle to popupChannelSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupChannelSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupChannelSelector

    if ~handles.flagDataLoaded
        return;
    end
    
    handles.dataDisplay.curChannelId = get(handles.popupChannelSelector,'Value');
    
    % Update handles structure
    guidata(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes during object creation, after setting all properties.
function popupChannelSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupChannelSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckboxMacrophageMask.
function CheckboxMacrophageMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxMacrophageMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxMacrophageMask

    % change mask use mode
    handles.flagShowMacrophageMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes on button press in CheckboxDrugMask.
function CheckboxDrugMask_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxDrugMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxDrugMask

    % change mask use mode
    handles.flagShowDrugMask = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --------------------------------------------------------------------
function File_PerformColocalizationAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_PerformColocalizationAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~handles.flagDataLoaded
        return;
    end

    handles = PerformColocalizationAnalysis(hObject, handles);

    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
    % Update Cell Descriptors
    UpdateCellDescriptors(handles);


% --------------------------------------------------------------------
function File_Colocalization_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Colocalization_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
