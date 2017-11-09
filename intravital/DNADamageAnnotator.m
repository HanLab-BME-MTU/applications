function varargout = DNADamageAnnotator(varargin)
%DNADAMAGEANNOTATOR M-file for DNADamageAnnotator.fig
%      DNADAMAGEANNOTATOR, by itself, creates a new DNADAMAGEANNOTATOR or raises the existing
%      singleton*.
%
%      H = DNADAMAGEANNOTATOR returns the handle to a new DNADAMAGEANNOTATOR or the handle to
%      the existing singleton*.
%
%      DNADAMAGEANNOTATOR('Property','Value',...) creates a new DNADAMAGEANNOTATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to DNADamageAnnotator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DNADAMAGEANNOTATOR('CALLBACK') and DNADAMAGEANNOTATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DNADAMAGEANNOTATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DNADamageAnnotator

% Last Modified by GUIDE v2.5 27-May-2014 16:45:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DNADamageAnnotator_OpeningFcn, ...
                   'gui_OutputFcn',  @DNADamageAnnotator_OutputFcn, ...
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

% --- Executes just before DNADamageAnnotator is made visible.
function DNADamageAnnotator_OpeningFcn(hObject, eventdata, handles, varargin)
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

    % Choose default command line output for DNADamageAnnotator
    handles.output = hObject;

    % make sure weka is in the path    
    AddWekaClassesToPath();
    
    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'DNADamageAnnotatorHistory.mat' );
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
    handles.flagShowFociBBox = get(handles.CheckboxFociBBox,'Value');
    handles.flagShowCellBBox = get(handles.CheckboxCellBBox, 'Value');

    handles.imarisAppCellSeg = [];
    handles.imarisAppCellPattern = [];
    handles.imarisAppCellSegCropped = [];
    
    % default parameters
    if strcmp(PARAMETERS.mode, 'analysis')
        f = rdir( fullfile(fileparts(mfilename('fullpath')), '**', 'regionMerging.model') );
        if ~isempty(f)
            handles.defaultParameters = NucleiSegmentationParametersGUI( 'defaultWithRegionMerging' );
            handles.defaultParameters.flagPerformRegionMerging = true;
            handles.defaultParameters.regionMergingModelFile = f(1).name;
        else
            handles.defaultParameters = NucleiSegmentationParametersGUI( 'defaultWithoutRegionMerging' );
        end
    else
        handles.defaultParameters = NucleiSegmentationParametersGUI( 'defaultRegionMergingTraining' );
    end
    
    handles.defaultParameters.flagParallelize = PARAMETERS.flagParallelize;
    handles.defaultParameters.flagDebugMode = PARAMETERS.flagDebugMode;
    
    handles.parameters = handles.defaultParameters;
    
    % open matlab pool for parallel processing    
    handles.parameters.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if handles.parameters.flagParallelize && ~handles.parameters.flagPoolOpenedAlready 
        matlabpool open;
    end
    
    % Set callbacks
    set(handles.DNADamageAnnotator, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes DNADamageAnnotator wait for user response (see UIRESUME)
    % uiwait(handles.DNADamageAnnotator);

% --- Outputs from this function are returned to the command line.
function varargout = DNADamageAnnotator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

function FociCountDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to FociCountDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FociCountDisplay as text
%        str2double(get(hObject,'String')) returns contents of FociCountDisplay as a double


% --- Executes during object creation, after setting all properties.
function FociCountDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FociCountDisplay (see GCBO)
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
    [fileName,pathName] = uigetfile( fullfile( handles.history.lastImageFileDir, '*.oif; *.oib' ), ...
                                     'Select the file contaning the nuclear marker channel - histone-2B' );   
    
    if ~fileName 
        return;
    end
    
    dataFilePath = fullfile( pathName, fileName );
    handles.history.lastImageFileDir = pathName;
    
    % load nucleus channel data
    PrettyPrintStepDescription( 'Loading Image Data' );
    
    hStatusDialog = waitbar(0, 'Loading Image Data');
    try        
        imageSeries = loadIntravitalDataset( dataFilePath  );    
        
        if (imageSeries.metadata.voxelSpacing(1)/ imageSeries.metadata.voxelSpacing(3)) >= 1
            imageSeries.metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
        end      
        
    catch err
        fprintf( 'ERROR: could not load image data data from file %s', dataFilePath );
        errordlg( 'Error loading image data from file' );
        closeStatusDialog(hStatusDialog);        
        return;
    end

    if numel(imageSeries(1).metadata.volSize) < 3 || ...
       imageSeries(1).metadata.volSize(3) < 2
   
        errordlg( 'Dataset is not 3D' );
        closeStatusDialog(hStatusDialog);        
        return;
        
    end

    if imageSeries(1).metadata.numChannels > 1

        hDisp = imseriesshow_multichannel( imageSeries(1).imageData(1,:), ...
                                           'spacing', imageSeries(1).metadata.voxelSpacing );  
        
        prompt = sprintf('data contains %d channels.\nEnter histone channel id:', imageSeries(1).metadata.numChannels );
        
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        options.Interpreter = 'none';

        channelId53BP1 = inputdlg( prompt, 'Histone channel selector', 1, {'1'}, options ); 
        
        if isempty(channelId53BP1) || isempty( str2num( channelId53BP1{1} ) )
            close(hDisp);
            closeStatusDialog(hStatusDialog);
            return;
        end

        channelId53BP1 = str2num(channelId53BP1{1});
        
        close( hDisp );
        
    else
        channelId53BP1 = 1;
    end
    
    % store image data in handles structures    
    clear handles.data;
    
    handles.flagDataLoaded = true;
    handles.data.dataFilePath = dataFilePath;
    handles.data.metadata = imageSeries(1).metadata;
    handles.data.imageData = imageSeries(1).imageData;
    handles.data.channelId53BP1 = channelId53BP1;
    
    set( handles.DNADamageAnnotator, 'Name', sprintf( 'DNA Damage Annotation - %s', handles.data.dataFilePath ) );
    
    % pre-compute data needed for display
    handles = ComputeDisplayData( handles );
    
    % close progress bar
    close( hStatusDialog );

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
    if handles.parameters.flagPerformRegionMerging
        regionMergingModelFile = handles.parameters.regionMergingModelFile;
    end
    
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints, ...
     handles.data.segAlgoParameters ] = segmentCellsInIntravitalData( handles.data.imageData{handles.data.channelId53BP1}, ...
                                                                      handles.data.metadata.voxelSpacing, ...                                                                      
                                                                      'flagParallelize', handles.parameters.flagParallelize, ...
                                                                      'flagDebugMode', handles.parameters.flagDebugMode, ...
                                                                      'cellDiameterRange', handles.parameters.cellDiameterRange, ...
                                                                      'thresholdingAlgorithm', 'BackgroudRemovalUsingMorphologicalOpening', ...
                                                                      'seedPointDetectionAlgorithm', handles.parameters.seedPointDetectionAlgorithm, ...
                                                                      'minCellVolume', handles.parameters.minCellVolume, ...
                                                                      'flagIgnoreCellsOnXYBorder', handles.parameters.flagIgnoreXYBorderCells, ...
                                                                      'regionMergingModelFile', regionMergingModelFile);

    [handles.dataDisplay.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.dataDisplay.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints, ones(3,3,3) );
    closeStatusDialog(hStatusDialog);
    
    % compute properties of each cell
    handles.data.cellStats = ComputeCellProperties( handles );
    
    % set current cell id to first cell
    handles.dataDisplay.curCellId = 1;
    handles.dataDisplay.curCellSliceId = round(handles.data.cellStats(handles.dataDisplay.curCellId).Centroid([2, 1, 3]));
    
    % Detect foci
    PerformFociSegmenatation(hObject, handles);
    
    % close progress bar
    closeStatusDialog(hStatusDialog);
    
% --------------------------------------------------------------------
function [cellStats] = ComputeCellProperties( handles )

    hStatusDialog = waitbar(0, 'Computing properties of segmented cells');
    
    cellStats = regionprops( handles.data.imLabelCellSeg, ...
                             'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
    
    for i = 1:numel(cellStats)

        % intensity descriptors
        cellPixelIntensities = handles.data.imageData{handles.data.channelId53BP1}( cellStats(i).PixelIdxList );
        cellStats(i).meanIntensity = mean( cellPixelIntensities );
        cellStats(i).stdIntensity = std( double(cellPixelIntensities) );
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
    
    nucleiMaskAlpha = 0.2;
    fociMaskAlpha = 0.6;
    seedMaskAlpha = 0.6;
    
    % display global cross section images  
    if handles.flagUseLOG               
        imGlobalXY = handles.dataDisplay.imageDataLOG{handles.data.channelId53BP1}( :, :, curCellSliceId(3) );       
        imGlobalXZ = squeeze( handles.dataDisplay.imageDataLOG{handles.data.channelId53BP1}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.dataDisplay.imageDataLOG{handles.data.channelId53BP1}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imLogDisplayRange;
    else        
        imGlobalXY = handles.data.imageData{handles.data.channelId53BP1}( :, :, curCellSliceId(3) );
        imGlobalXZ = squeeze( handles.data.imageData{handles.data.channelId53BP1}( curCellSliceId(1), :, : ) );
        imGlobalYZ = squeeze( handles.data.imageData{handles.data.channelId53BP1}( :, curCellSliceId(2), : ) );
        displayrange = handles.dataDisplay.imDisplayRange;
    end
    imGlobalXY = mat2gray(imGlobalXY, displayrange(1,:) );
    imGlobalXZ = mat2gray(imGlobalXZ', displayrange(1,:) );
    imGlobalYZ = mat2gray(imGlobalYZ, displayrange(1,:) );

    if handles.flagShowNucleiSegMask || handles.flagShowFociSegMask

        imGlobalXYMasks = {};
        imGlobalXZMasks = {};
        imGlobalYZMasks = {};
    
        maskAlphas = [];
        
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
        
        imsize = size( handles.data.imageData{handles.data.channelId53BP1} );
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
    imsize = size(handles.data.imageData{handles.data.channelId53BP1});
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
        imCellCropped = mat2gray( handles.dataDisplay.imageDataLOG{handles.data.channelId53BP1}(subinds{:}), handles.dataDisplay.imLogDisplayRange(1,:) );        
    else
        imCellCropped = mat2gray( handles.data.imageData{handles.data.channelId53BP1}(subinds{:}), handles.dataDisplay.imDisplayRange(1,:) );        
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

        if handles.flagShowFociSegMask

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
    imCurCellMIP = mat2gray( max(handles.data.imageData{handles.data.channelId53BP1}(subindsMIP{:}) .* imCellSegCropped, [], 3 ) );
    imHistoneMIPDisplay = repmat( imCurCellMIP, [1,1,3] );
    
    cla( handles.Axes_Histone_MIP, 'reset' );
    image( imHistoneMIPDisplay, 'Parent', handles.Axes_Histone_MIP );
    set( handles.Axes_Histone_MIP, 'XTickLabel', [], 'YTickLabel', [] );
    
    % highlight currently selected foci
    if handles.flagShowFociBBox && numel(handles.data.fociStats) > 0
        
        curFociStats = handles.data.fociStats(handles.dataDisplay.curFociId);
        imsize = size( handles.data.imageData{handles.data.channelId53BP1} );
        spacing = handles.data.metadata.voxelSpacing;
        curFociLoc = ind2submat(imsize, curFociStats.PixelLocationIndex);
        markerSize = 12.0;
        lineWidth = 1.5;
        
        hold( handles.Axes_Global_XY, 'on' );
            plot( handles.Axes_Global_XY, curFociLoc(2), curFociLoc(1), 'ms', ...
                  'LineWidth', lineWidth, 'MarkerSize', markerSize ); 
        hold( handles.Axes_Global_XY, 'off' );    

        hold( handles.Axes_Global_XZ, 'on' );
            plot( handles.Axes_Global_XZ, curFociLoc(2), curFociLoc(3), 'ms', ...
                  'LineWidth', lineWidth, 'MarkerSize', markerSize ); 
        hold( handles.Axes_Global_XZ, 'off' );    

        hold( handles.Axes_Global_YZ, 'on' );
            plot( handles.Axes_Global_YZ, curFociLoc(3), curFociLoc(1), 'ms', ...
                  'LineWidth', lineWidth, 'MarkerSize', markerSize ); 
        hold( handles.Axes_Global_YZ, 'off' );    

        hold( handles.Axes_Nucleus_XY, 'on' );
            plot( handles.Axes_Nucleus_XY, ...
                  curFociLoc(2) - min(subinds{2}) + 1, curFociLoc(1) - min(subinds{1}) + 1, ...
                  'ms', 'LineWidth', lineWidth, 'MarkerSize', markerSize ); 
        hold( handles.Axes_Nucleus_XY, 'off' );    

        hold( handles.Axes_Histone_MIP, 'on' );
            plot( handles.Axes_Histone_MIP, ...
                  curFociLoc(2) - min(subinds{2}) + 1, curFociLoc(1) - min(subinds{1}) + 1, ...
                  'ms', 'LineWidth', lineWidth, 'MarkerSize', markerSize ); 
        hold( handles.Axes_Histone_MIP, 'off' );    
        
    end
    
    % display cell id    
    strtmp = sprintf('%d / %d', handles.dataDisplay.curFociId, numel(handles.data.fociStats) );
    set(handles.FociCountDisplay, 'String', strtmp);    
    
    % display cell slice id
    strtmp = sprintf('x-slice: %.3d / %.3d\ny-slice: %.3d / %.3d\nz-slice: %.3d / %.3d', ...
                     curCellSliceId(2), imsize(2), ...
                     curCellSliceId(1), imsize(1), ...
                     curCellSliceId(3), imsize(3) );
    set(handles.EditboxSliceId, 'String', strtmp);    

% --------------------------------------------------------------------
function UpdateFociDescriptors(handles)

    curFociStats = handles.data.fociStats( handles.dataDisplay.curFociId );
    propNames = fieldnames(curFociStats);
    strCellDescription = '';
    for i = 1:numel(propNames)    
       
       curFieldName =  propNames{i};
       curField = curFociStats.(propNames{i});
       
       if ismember(curFieldName, {'PixelLocationIndex', 'isValid'})
           continue;
       end
       
       if numel(curField) > 1
           strCellDescription = [strCellDescription, ...
                                 sprintf('\n%s: [%s]\n', curFieldName, sprintf(' %.3f ', curField))];
       else           
           strCellDescription = [strCellDescription, ...
                                 sprintf('\n%s: %.3f\n', curFieldName, curField)];
       end
       
    end
    
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
function File_Load_Annotation_Callback(hObject, eventdata, handles)
% hObject    handle to File_Load_Annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % get annotation file from user
    if isfield( handles.history, 'lastAnalysisFileDir' )
        [fileName,pathName] = uigetfile( fullfile(handles.history.lastAnalysisFileDir, 'DNADamageAnnotation.mat'), 'Select analysis file' );   
    else
        [fileName,pathName] = uigetfile( fullfile(pwd, 'DNADamageAnnotation.mat'), 'Select analysis file' );   
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
            handles.data.dataFilePath = analysisData.dataFilePath{1};
        else
            handles.data.dataFilePath = analysisData.dataFilePath;
        end
        handles.data.metadata = analysisData.metadata;
        handles.data.imageData = analysisData.imageData;
        for i = 1:numel(handles.data.imageData)
            handles.data.imageData{i} = double(handles.data.imageData{i});
        end
        handles.data.channelId53BP1 = analysisData.channelId53BP1;
        
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
        handles.data.unannotatedFociList = analysisData.unannotatedFociList;
        
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
        
    % data ready for display
    handles.flagDataLoaded = true;
    
    % set foci id to next unannotated foci
    if ~isempty(handles.data.unannotatedFociList)
        handles.dataDisplay.curFociId = handles.data.unannotatedFociList(1);
    else
        handles.dataDisplay.curFociId = 1;
    end
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % change window name
    set( handles.DNADamageAnnotator, 'Name', sprintf( 'DNA Damage Annotation - %s', analysisFile ) );

    % update foci selector list
    set(handles.poplistFociClassSelector, 'value', 1);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
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
    [pathstr, name, ext] = fileparts( handles.data.dataFilePath );
    outputDir = fullfile(outputDir, name);

    if ~isdir( outputDir )
        mkdir( outputDir );
    end
    
    data = handles.data;
    for i = 1:numel(handles.data.imageData)
        data.imageData{i} = uint16(data.imageData{i});
    end
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save(fullfile(outputDir, 'DNADamageAnnotation.mat'), '-struct', 'data' );    
    writetable(struct2table(handles.data.fociStats), fullfile(outputDir, 'fociStats.csv'));
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


% --- Executes when user attempts to file_close DNADamageAnnotator.
function DNADamageAnnotator_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnnotator (see GCBO)
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

% --- Executes on button press in ShowPreviousFoci.
function ShowPreviousFoci_Callback(hObject, eventdata, handles)
% hObject    handle to ShowPreviousFoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    % decrement foci id
    handles.dataDisplay.curFociId = max(1, handles.dataDisplay.curFociId - 1);
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
% --- Executes on button press in ShowNextFoci.
function ShowNextFoci_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNextFoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % increment foci id
    handles.dataDisplay.curFociId = min( numel( handles.data.fociStats ), handles.dataDisplay.curFociId + 1 );
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
% --- Executes on button press in ShowLastFoci.
function ShowLastFoci_Callback(hObject, eventdata, handles)
% hObject    handle to ShowLastFoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % set foci id as last one
    handles.dataDisplay.curFociId = numel( handles.data.fociStats );
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
% --- Executes on button press in ShowFirstFoci.
function ShowFirstFoci_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFirstFoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end

    % set foci id
    handles.dataDisplay.curFociId = 1;
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);

    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
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

    handles.parameters = NucleiSegmentationParametersGUI( handles.parameters );

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

% --- Executes during object deletion, before destroying properties.
function DNADamageAnnotator_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    if isfield(handles, 'history')
        history = handles.history;
        [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
        historyFile = fullfile( pathstr, 'DNADamageAnnotatorHistory.mat' );
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
    
    imsize = size(handles.data.imageData{handles.data.channelId53BP1});

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
    
    curCellStats = handles.data.cellStats( handles.dataDisplay.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.cellDisplaySize] );
    
    % create crop indices
    subinds = cell(1,3);
    imsize = size(handles.data.imageData{handles.data.channelId53BP1});
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
                                 'spacing', handles.data.metadata.voxelSpacing);
                             
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
               
            curFociInd = find(curCellStats.foci == handles.dataDisplay.curFociId);    
            imvis.AddSpots(fociSeedPointLocations(curFociInd, :), 0, ...
                           'radii', [curCellFociStats(curFociInd).Radius], ... 
                           'name', 'SelectedPuncta', 'color', [1, 1, 0]);
            
        end    
    
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
    imvis = ImarisDataVisualizer( cat(4, handles.data.imageData{:}), 'spacing', handles.data.metadata.voxelSpacing );
    handles.imarisAppCellSeg = imvis;
    
    hSegmentation = imvis.AddDataContainer();
    hPuncta = imvis.AddDataContainer();
    
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

            % Unannotated
            ind = find( [handles.data.fociStats.isValid] < 0 );
            if numel(ind) > 0
                imvis.AddSpots(fociSeedPointLocations(ind,:), 0, ...
                               'radii', [handles.data.fociStats(ind).Radius], ... 
                               'hContainer', hPuncta, ...
                               'name', 'UnannotatedPuncta', 'color', [0, 0, 1]);
            end
            
            % Bad_Detection
            ind = find( [handles.data.fociStats.isValid] == 0 );
            if numel(ind) > 0
                imvis.AddSpots(fociSeedPointLocations(ind,:), 0, ...
                               'radii', [handles.data.fociStats(ind).Radius], ... 
                               'hContainer', hPuncta, ...
                               'name', 'Bad_Detection', 'color', [1, 0, 0]);
            end
            
            % Good_Detection
            ind = find( [handles.data.fociStats.isValid] > 0 );
            if numel(ind) > 0
                imvis.AddSpots(fociSeedPointLocations(ind,:), 0, ...
                               'radii', [handles.data.fociStats(ind).Radius], ... 
                               'hContainer', hPuncta, ...
                               'name', 'Good_Detection', 'color', [0, 1, 0]);
            end
        end    
        
    % Update handles structure
    guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function FociDescriptorsTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FociDescriptorsTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function LabelCellDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelCellDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function DNADamageAnnotator_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnnotator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{handles.data.channelId53BP1});
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

    handles.dataDisplay.curCellSliceId = newCellSliceId;
    pixelPos = num2cell( pixelPos );
    indViewSlicer = { ':', ':', ':' };
    indViewSlicer{viewId} = handles.dataDisplay.curCellSliceId(viewId);
    
    % get cell id under mouse click
    imCurLabelCellSeg = squeeze( handles.data.imLabelCellSeg( indViewSlicer{:} ) );
    
    curCellId = imCurLabelCellSeg( pixelPos{:} );
    if curCellId >= 1
        
        % set cell id
        handles.dataDisplay.curCellId = curCellId;
        
    end

    % get cell id under mouse click
    imCurLabelFociSeg = squeeze( handles.data.imLabelFociSeg( indViewSlicer{:} ) );
    
    curFociId = imCurLabelFociSeg( pixelPos{:} );
    if curFociId >= 1
        
        % set foci id
        handles.dataDisplay.curFociId = curFociId;

        % update cell pattern listbox
        set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
        
    end
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);    
    
% -----------------------------------------------------------------
function PerformFociSegmenatation(hObject, handles)

    if ~handles.flagDataLoaded
        return;
    end

    PrettyPrintStepDescription( 'Detecting Foci' );
    
    [handles.data.fociStats, ...
     handles.data.imFociSeedPoints, ...
     handles.data.imLabelFociSeg, ...
     handles.data.fociDetectionParameters ] = segmentFociInsideNuclei( handles.data.imageData{handles.data.channelId53BP1}, ...
                                                                       min(handles.data.metadata.voxelSpacing) * [3, 7], ...                      
                                                                       'spacing', handles.data.metadata.voxelSpacing, ...
                                                                       'roiMask', handles.data.imLabelCellSeg);
    
    [handles.dataDisplay.imFociSegRGBMask, handles.data.FociSegColorMap] = label2rgbND( handles.data.imLabelFociSeg );
    handles.dataDisplay.imFociSeedPoints = imdilate(handles.data.imFociSeedPoints, ones(3,3,3));
    
    for cid = 1:numel(handles.data.cellStats)
        
        curCellPixInd = find(handles.data.imLabelCellSeg == cid);
        curCellFoci = unique(handles.data.imFociSeedPoints(curCellPixInd));
        curCellFoci = curCellFoci(curCellFoci > 0);
        handles.data.cellStats(cid).foci = curCellFoci;
        handles.data.cellStats(cid).fociCount = numel(curCellFoci);
        
    end

    % add all foci to unannotated list
    handles.data.unannotatedFociList = 1:numel(handles.data.fociStats);
    
    if numel(handles.data.fociStats) > 0
        
        initAnnotationStatus = num2cell(-1 * ones(size(handles.data.fociStats)));
        [handles.data.fociStats.isValid] = initAnnotationStatus{:};
        
        % set foci id
        handles.dataDisplay.curFociId = 1;
        handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
        handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);

        % update cell pattern listbox
        set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);

        % enable Load next unannotated focu button 
        set( handles.ShowNextUnannotatedFoci, 'Enable', 'on' );
        
    else
        handles.dataDisplay.curFociId = 0;
        set( handles.ShowNextUnannotatedFoci, 'Enable', 'off' );
    end
    
    % update foci selector list
    set(handles.poplistFociClassSelector, 'value', 1);
    
    % Update handles structure
    guidata(hObject, handles);
   
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
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

    PerformFociSegmenatation(hObject, handles);


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

    [~, fname, ~] = fileparts(handles.data.dataFilePath);
    fname = strrep(fname,'_','\_');
    text(0.5, 0.95, {fname, sprintf('#cells = %d, #puncta= %d', numCells, sum(punctaCounts))}, ...
         'units', 'normalized', 'horizontalalignment', 'center');
    
    title({'Distribution of Number of Puncta In a Cell'}, ...
          'FontSize', 12.0, 'FontWeight', 'bold');
    xlabel('Number of puncta per cell');
    ylabel('Proportion of cells');
    ylim([0, 1]);
    
    grid on;
    
    
% --- Executes on button press in btnShowPreviousFociInSelectedClass.
function btnShowPreviousFociInSelectedClass_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowPreviousFociInSelectedClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    fociAnnotation = [handles.data.fociStats.isValid];
    
    indPrevFociInCurClass = find(fociAnnotation(handles.dataDisplay.curFociId-1:-1:1) == ...
                                 get(handles.poplistFociClassSelector, 'Value')-2);
    
    if isempty( indPrevFociInCurClass )
        return;        
    end
    
    % navigate to next foci automatically
    handles.dataDisplay.curFociId = handles.dataDisplay.curFociId - indPrevFociInCurClass(1);
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update foci quality listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);

% --- Executes on button press in btnShowNextFociInSelectedClass.
function btnShowNextFociInSelectedClass_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowNextFociInSelectedClass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end
    
    fociAnnotation = [handles.data.fociStats.isValid];
    
    indNextFociInCurClass = find(fociAnnotation(handles.dataDisplay.curFociId+1:end) == ...
                                 get(handles.poplistFociClassSelector, 'Value')-2);
    
    if isempty( indNextFociInCurClass )
        return;        
    end
    
    % navigate to next foci automatically
    handles.dataDisplay.curFociId = handles.dataDisplay.curFociId + indNextFociInCurClass(1);
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update foci quality listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);

% --- Executes on selection change in poplistFociClassSelector.
function poplistFociClassSelector_Callback(hObject, eventdata, handles)
% hObject    handle to poplistFociClassSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistFociClassSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistFociClassSelector

    % first check if data has been loaded and some cells are present
    if ~handles.flagDataLoaded
        return;
    end

    fociAnnotation = [handles.data.fociStats.isValid];
    indFirstFociInClass = find(fociAnnotation == get(hObject,'Value')-2);
    
    if isempty( indFirstFociInClass )
        return;        
    end
    
    % navigate to next foci automatically
    handles.dataDisplay.curFociId = indFirstFociInClass(1);
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function poplistFociClassSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistFociClassSelector (see GCBO)
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


% --- Executes on key press with focus on DNADamageAnnotator or any of its controls.
function DNADamageAnnotator_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to DNADamageAnnotator (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded || isempty(handles.data.fociStats) 
        return;
    end

    if lower(eventdata.Character) == 'f'
        handles.flagShowFociSegMask = ~handles.flagShowFociSegMask;
        set(handles.CheckboxFociSegMask, 'Value', handles.flagShowFociSegMask);
    end
        
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);

% --- Executes on selection change in ListboxFociDetectionQualitySelector.
function ListboxFociDetectionQualitySelector_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxFociDetectionQualitySelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListboxFociDetectionQualitySelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListboxFociDetectionQualitySelector

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end    

    curFociId = handles.dataDisplay.curFociId;
    handles.data.fociStats(curFociId).isValid = get(hObject, 'Value') - 2;
    
    % if pattern type is not equal to 'None' remove it from the list of
    % unannotated foci, if not add it to the list
    if handles.data.fociStats(curFociId).isValid < 0
        handles.data.unannotatedFociList = union(handles.data.unannotatedFociList, curFociId); 
    else
        handles.data.unannotatedFociList = setdiff(handles.data.unannotatedFociList, curFociId); 
    end

    if isempty( handles.data.unannotatedFociList )
        set( handles.ShowNextUnannotatedFoci, 'Enable', 'off' );
    else        
        set( handles.ShowNextUnannotatedFoci, 'Enable', 'on' );
    end

    % navigate to next foci automatically
    handles.dataDisplay.curFociId = min(numel( handles.data.fociStats ), handles.dataDisplay.curFociId + 1 );
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);

    % update foci quality listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);
    
% --- Executes during object creation, after setting all properties.
function ListboxFociDetectionQualitySelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListboxFociDetectionQualitySelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowNextUnannotatedFoci.
function ShowNextUnannotatedFoci_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNextUnannotatedFoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.flagDataLoaded || isempty(handles.data.unannotatedFociList) 
        return;
    end

    % set foci id to next unannotated foci
    handles.dataDisplay.curFociId = handles.data.unannotatedFociList(1);
    handles.dataDisplay.curCellId = handles.data.imLabelCellSeg(handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    handles.dataDisplay.curCellSliceId = ind2submat(size(handles.data.imFociSeedPoints), handles.data.fociStats(handles.dataDisplay.curFociId).PixelLocationIndex);
    
    % update cell pattern listbox
    set(handles.ListboxFociDetectionQualitySelector, 'Value', handles.data.fociStats(handles.dataDisplay.curFociId).isValid + 2);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateFociDescriptors(handles);


% --- Executes on button press in CheckboxFociBBox.
function CheckboxFociBBox_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxFociBBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckboxFociBBox

    % change mask use mode
    handles.flagShowFociBBox = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
