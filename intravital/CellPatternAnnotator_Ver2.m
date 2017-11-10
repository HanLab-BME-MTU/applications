function varargout = CellPatternAnnotator_Ver2(varargin)
%CELLPATTERNANNOTATOR_VER2 M-file for CellPatternAnnotator_Ver2.fig
%      CELLPATTERNANNOTATOR_VER2, by itself, creates a new CELLPATTERNANNOTATOR_VER2 or raises the existing
%      singleton*.
%
%      H = CELLPATTERNANNOTATOR_VER2 returns the handle to a new CELLPATTERNANNOTATOR_VER2 or the handle to
%      the existing singleton*.
%
%      CELLPATTERNANNOTATOR_VER2('Property','Value',...) creates a new CELLPATTERNANNOTATOR_VER2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CellPatternAnnotator_Ver2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CELLPATTERNANNOTATOR_VER2('CALLBACK') and CELLPATTERNANNOTATOR_VER2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CELLPATTERNANNOTATOR_VER2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellPatternAnnotator_Ver2

% Last Modified by GUIDE v2.5 27-Mar-2013 17:20:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellPatternAnnotator_Ver2_OpeningFcn, ...
                   'gui_OutputFcn',  @CellPatternAnnotator_Ver2_OutputFcn, ...
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

% --- Executes just before CellPatternAnnotator_Ver2 is made visible.
function CellPatternAnnotator_Ver2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

    % Choose default command line output for CellPatternAnnotator_Ver2
    handles.output = hObject;

    % load history
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
    if exist( historyFile, 'file' )
       handles.data.history = load( historyFile );
    else
       handles.data.history.lastAnalyzedNucleusDir = pathstr;
       handles.data.history.lastAnalyzedRedGreenDir = pathstr;
    end

    % Initialize global variables
    handles.data.cellDisplaySize = 70;     
    
    handles.data.flagDataLoaded = false;
    handles.data.flagUseLOG = true;
    handles.data.flagShowGlobalSegMask = true;
    handles.data.flagShowCellBBox = true;
    
    handles.data.flagShowLocalMIP = false;
    handles.data.flagShowRedChannel = true;
    handles.data.flagShowGreenChannel = true;
    handles.data.flagShowNucleusChannel = true;    
    
    handles.data.defaultCellPatternTypes = cellstr( get(handles.ListboxCellPatternSelector, 'String') );
    handles.data.cellPatternTypes = cellstr( get(handles.ListboxCellPatternSelector, 'String') );
    
    handles.imarisApp = [];
    handles.imarisAppCellCropped = [];
    handles.imarisAppCellSegCropped = [];
    
    % default parameters
    handles.defaultParameters.cellDiameterRange = [12, 20]; 
    handles.defaultParameters.minCellVolume = 400; 
    handles.defaultParameters.seedPointDetectionAlgorithm = 'AdaptiveMultiscaleLoG'; 
    
    handles.parameters = handles.defaultParameters;
    
    % open matlab pool for parallel processing
    handles.data.flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
    if ~handles.data.flagPoolOpenedAlready 
        matlabpool open;
    end            
    
    % Set callbacks
    set(handles.CellPatternAnnotator, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes CellPatternAnnotator_Ver2 wait for user response (see UIRESUME)
    % uiwait(handles.CellPatternAnnotator_Ver2);

% --- Outputs from this function are returned to the command line.
function varargout = CellPatternAnnotator_Ver2_OutputFcn(hObject, eventdata, handles)
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
    if ~handles.data.flagDataLoaded 
        return;
    end    

    curCellId = handles.data.curCellId;
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
    
    % Update handles structure
    guidata(hObject, handles);
    
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
    [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastAnalyzedNucleusDir, '*.oif' ), 'Select the data file for the nucleus channel' );   
    
    if ~fileName 
        return;
    end
    
    handles.data.dataFilePath{1} = fullfile( pathName, fileName );
    handles.data.history.lastAnalyzedNucleusDir = pathName;
    handles.data.history.lastAnalyzedRedGreenDir = pathName;
    
    % ask the use to select the oif file for the red/green channel
    [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastAnalyzedRedGreenDir, '*.oif' ), 'Select the data file for the red/green channel' );   
        
    if ~fileName 
        handles.data.flagRedGreenChannel = false;        
    else    
        handles.data.dataFilePath{2} = fullfile( pathName, fileName );        
        handles.data.flagRedGreenChannel = true;        
        handles.data.history.lastAnalyzedRedGreenDir = pathName;
    end
    
    % load nucleus channel data
    PrettyPrintStepDescription( 'Loading Nuclear Channel Data' );
    
    hStatusDialog = waitbar(0, 'Loading Nuclear Channel Data');
    try        
        imageSeriesNucleus = loadIntravitalDataset( handles.data.dataFilePath{1}  );    
        if (imageSeriesNucleus(1).metadata.voxelSpacing(1)/ imageSeriesNucleus(1).metadata.voxelSpacing(3)) >= 1
            imageSeriesNucleus(1).metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
        end        
        imageSeriesNucleus(1).metadata.channelColors = [ 0 , 0 , 1 ];
    catch err
        fprintf( 'ERROR: could not load nucleus channel data from file %s', handles.data.dataFilePath{1} );
        return;
    end
    
    % load red/green channel data
    if handles.data.flagRedGreenChannel
    
        PrettyPrintStepDescription( 'Loading Red/Green Channel Data' );
        
        waitbar(0, hStatusDialog, 'Loading Red/Green Channel Data');
        try
            
            imageSeriesRedGreen = loadIntravitalDataset( handles.data.dataFilePath{2}  );    
            
            if imageSeriesRedGreen(1).metadata.channelExcitationWavelength(1) > imageSeriesRedGreen(1).metadata.channelExcitationWavelength(2)
                
                % first channel is not green so swap it 
                warning( 'GFP was not the first channel in confocal data. It will be swapped to the first place');
                imageSeriesRedGreen(1).imageData = imageSeriesRedGreen(1).imageData(:,[2,1]);
                imageSeriesRedGreen(1).metadata.channelExcitationWavelength = imageSeriesRedGreen(1).metadata.channelExcitationWavelength([2,1]);
                
            end
                
            imageSeriesRedGreen(1).metadata.channelColors = [ 0 1 0; 1 0 0 ];
            
        catch err

            err
            fprintf( 'ERROR: could not load nucleus channel data from file %s', handles.data.dataFilePath{1} );
            errordlg( sprintf( 'Could not load nucleus channel data from file %s', handles.data.dataFilePath{1} ) ); 
            return;
        end
        
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
        
        set( handles.LabelAllChannels, 'String', 'All Channel Overlay' );
        set( handles.CheckboxShowNucleusChannel, 'Visible', 'on' );
        set( handles.CheckboxShowRedChannel, 'Visible', 'on' );
        set( handles.CheckboxShowGreenChannel, 'Visible', 'on' );            
        set( handles.CheckboxLocalMIP, 'Visible', 'on' );            
        
    else
        
        set( handles.LabelAllChannels, 'String', 'Nucleus Channel MIP' );
        set( handles.CheckboxShowNucleusChannel, 'Visible', 'off' );
        set( handles.CheckboxShowRedChannel, 'Visible', 'off' );
        set( handles.CheckboxShowGreenChannel, 'Visible', 'off' );
        set( handles.CheckboxLocalMIP, 'Visible', 'off' );      
        
    end
    
    % store image data in handles structures    
    handles.data.flagDataLoaded = true;
    handles.data.metadata = imageSeriesNucleus(1).metadata;
    
    set( handles.CellPatternAnnotator, 'Name', sprintf( 'Cell Pattern Annotator - %s', handles.data.dataFilePath{1} ) );
    
    if handles.data.flagRedGreenChannel

        % merge channel wavelengths
        handles.data.metadata.channelExcitationWavelength = [ imageSeriesNucleus(1).metadata.channelExcitationWavelength, ...
                                                              imageSeriesRedGreen(1).metadata.channelExcitationWavelength ];

        handles.data.metadata.channelColors = cat( 1 , imageSeriesNucleus(1).metadata.channelColors, ...
                                                       imageSeriesRedGreen(1).metadata.channelColors );             
        
        % correct stage-shift        
        PrettyPrintStepDescription( 'Correcting stage-shift between confocal and two-photon data' );
        
        waitbar(0, hStatusDialog, 'Correcting Stage-shift between confocal and two-photon data');
        [handles.data.imageData, ...
         handles.data.imRegValidMask] = CorrectTwoPhotonConfocalStageShift( imageSeriesNucleus(1).imageData, imageSeriesNucleus(1).metadata.voxelSpacing, ...
                                                                            imageSeriesRedGreen(1).imageData, imageSeriesRedGreen(1).metadata.voxelSpacing, ...
                                                                            false );
    
        % compute display ranges and log-transformed images
        for i = 1:3
            handles.data.imDisplayRange(i,:) = ComputeImageDynamicRange( handles.data.imageData{i}, 98.0 );   
            handles.data.imageDataLOG{i} = ComputeImageLogTransformForDisplay( handles.data.imageData{i} );
            handles.data.imLogDisplayRange(i,:) = [ min(handles.data.imageDataLOG{i}(:)), max(handles.data.imageDataLOG{i}(:))];   
        end
    
    else        
        handles.data.imageData{1} = imageSeriesNucleus(1).imageData{1,1};        
        handles.data.imDisplayRange = ComputeImageDynamicRange(handles.data.imageData{1}, 98.0);    
        handles.data.imageDataLOG{1} = ComputeImageLogTransformForDisplay( handles.data.imageData{1} );
        handles.data.imLogDisplayRange = [ min(handles.data.imageDataLOG{1}(:)), max(handles.data.imageDataLOG{1}(:))];    
    end

    % close progress bar
    close( hStatusDialog );

    % Update handles structure
    guidata(hObject, handles);

    % Run analysis
    RunAnalysis(hObject, handles);
    
% --------------------------------------------------------------------    
function RunAnalysis(hObject, handles)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end    

    % segment cells
    PrettyPrintStepDescription( 'Running Cell Sementation Algorithm' );    
    hStatusDialog = waitbar(0, 'Segmentating Cells');
   
    [handles.data.imLabelCellSeg, ...
     handles.data.imCellSeedPoints, ...
     handles.data.segAlgoParameters ] = segmentCellsInIntravitalData( handles.data.imageData{1}, ...
                                                                      handles.data.metadata.voxelSpacing, ...                                                                      
                                                                      'flagParallelize', true, ...
                                                                      'cellDiameterRange', handles.parameters.cellDiameterRange, ...
                                                                      'thresholdingAlgorithm', 'OtsuGlobalSliceBySliceHybrid', ...
                                                                      'seedPointDetectionAlgorithm', handles.parameters.seedPointDetectionAlgorithm, ...
                                                                      'minCellVolume', handles.parameters.minCellVolume );

    [handles.data.imCellSegRGBMask, handles.data.CellSegColorMap] = label2rgbND(handles.data.imLabelCellSeg);
    handles.data.imCellSeedPoints = imdilate( handles.data.imCellSeedPoints, ones(3,3,3) );
    close(hStatusDialog);
    
    % compute properties of each cell
    PrettyPrintStepDescription( 'Computing Properties of Segmentated Cells' );    
    
    hStatusDialog = waitbar(0, 'Computing properties of segmented cells');
    
    handles.data.cellStats = regionprops( handles.data.imLabelCellSeg, ...
                                          'Centroid', 'BoundingBox', 'Area', 'PixelIdxList' );
    
    
    for i = 1:numel(handles.data.cellStats)

        % intensity descriptors
        cellPixelIntensities = handles.data.imageData{1}( handles.data.cellStats(i).PixelIdxList );
        handles.data.cellStats(i).meanIntensity = mean( cellPixelIntensities );
        handles.data.cellStats(i).stdIntensity = std( cellPixelIntensities );
        handles.data.cellStats(i).minIntensity = min( cellPixelIntensities );
        handles.data.cellStats(i).maxIntensity = max( cellPixelIntensities );
        
        % volume
        handles.data.cellStats(i).AreaPhysp = handles.data.cellStats(i).Area * prod( handles.data.metadata.voxelSpacing );
    
        % fit ellipsoid and store its radii
        [yind, xind, zind] = ind2sub( size(handles.data.imageData{1}), handles.data.cellStats(i).PixelIdxList ); 
        ptCell = [xind, yind, zind] .* repmat( handles.data.metadata.voxelSpacing, [numel(xind), 1] );
        ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
        [U, S, V] = svd( (ptCell' * ptCell) / size(ptCell,1) );
        
        handles.data.cellStats(i).ellipsoidRadiusPhysp = zeros(1,3); 
        for j = 1:3
            handles.data.cellStats(i).ellipsoidRadiusPhysp(j) = 2 * sqrt(S(j,j)); % eigen-values are a measure of variance
        end
        
        % update progress
        hStatusDialog = waitbar(i/numel(handles.data.cellStats), hStatusDialog, 'Computing properties of segmented cells');    

    end

    % initialize annotation to none
    for i = 1:numel(handles.data.cellStats)
        
        handles.data.cellStats(i).cellPatternId = 1; %Not Annotated
        handles.data.cellStats(i).cellPatternType = handles.data.cellPatternTypes{1};
    
    end

    % set current cell id to first cell
    handles.data.curCellId = 1;
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));
    
    % add all cells to unannotated list
    handles.data.unannotatedCellList = 1:numel(handles.data.cellStats);
    set( handles.ShowNextUnannotatedCell, 'Enable', 'on' );
    
    % close progress bar
    close(hStatusDialog);
    
    % Update handles structure
    guidata(hObject, handles);
   
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);


% --------------------------------------------------------------------
function File_Run_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to File_Run_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end    

    RunAnalysis(hObject, handles);

% --------------------------------------------------------------------
function UpdateCellDisplay(handles)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end    
    
    curCellStats = handles.data.cellStats( handles.data.curCellId );        
    curCellSliceId = handles.data.curCellSliceId;
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.data.cellDisplaySize] );
    
    % display global MIPs    
    if handles.data.flagUseLOG               
        imGlobalXY = handles.data.imageDataLOG{1}( :, :, curCellSliceId );       
        imGlobalXZ = squeeze( handles.data.imageDataLOG{1}( curCellCentroid(2), :, : ) );
        imGlobalYZ = squeeze( handles.data.imageDataLOG{1}( :, curCellCentroid(1), : ) );
        displayrange = handles.data.imLogDisplayRange;
    else        
        imGlobalXY = handles.data.imageData{1}( :, :, curCellSliceId );
        imGlobalXZ = squeeze( handles.data.imageData{1}( curCellCentroid(2), :, : ) );
        imGlobalYZ = squeeze( handles.data.imageData{1}( :, curCellCentroid(1), : ) );
        displayrange = handles.data.imDisplayRange;
    end
    imGlobalXY = mat2gray(imGlobalXY, displayrange(1,:) );
    imGlobalXZ = mat2gray(imGlobalXZ', displayrange(1,:) );
    imGlobalYZ = mat2gray(imGlobalYZ, displayrange(1,:) );

    if handles.data.flagShowGlobalSegMask 
        
        imGlobalXYSegMaskRGB = squeeze( handles.data.imCellSegRGBMask( :, :, curCellSliceId, : ) );
        
        imGlobalXZSegMaskRGB = squeeze( handles.data.imCellSegRGBMask(curCellCentroid(2), :, :, :) );
        imGlobalXZSegMaskRGB = permute(imGlobalXZSegMaskRGB, [2,1,3] );
        
        imGlobalYZSegMaskRGB = squeeze( handles.data.imCellSegRGBMask(:, curCellCentroid(1), :, :) );

        imGlobalXYDisplay = genImageRGBMaskOverlay( imGlobalXY, imGlobalXYSegMaskRGB, 0.2 );
        imGlobalXZDisplay = genImageRGBMaskOverlay( imGlobalXZ, imGlobalXZSegMaskRGB, 0.2 );
        imGlobalYZDisplay = genImageRGBMaskOverlay( imGlobalYZ, imGlobalYZSegMaskRGB, 0.2 );
        
    else
        
        imGlobalXYDisplay = repmat( imGlobalXY, [1,1,3] );
        imGlobalXZDisplay = repmat( imGlobalXZ, [1,1,3] );
        imGlobalYZDisplay = repmat( imGlobalYZ, [1,1,3] );
        
    end
        
    cla( handles.Axes_Global_XY );
    cla( handles.Axes_Global_XZ );
    cla( handles.Axes_Global_YZ );
    
    image( imGlobalXYDisplay, 'Parent', handles.Axes_Global_XY );           
    image( imGlobalXZDisplay, 'Parent', handles.Axes_Global_XZ );
    image( imGlobalYZDisplay, 'Parent', handles.Axes_Global_YZ );           

    set( handles.Axes_Global_XY, 'XTickLabel', [], 'YTickLabel', [] );
    set( handles.Axes_Global_XZ, 'XTickLabel', [], 'YTickLabel', [] );
    set( handles.Axes_Global_YZ, 'XTickLabel', [], 'YTickLabel', [] );        
    
    % draw bounding box around each cell
    if handles.data.flagShowCellBBox
        
        imsize = size( handles.data.imageData{1} );
        hold( handles.Axes_Global_XY, 'on' );
            w = curCellDisplaySize([1,1]);
            ptCorner = round(curCellCentroid(1:2) - 0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1)-1, 0; w-1 ; 0, w(2)-1; 0, 0 ];        
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(1), 2 ) = imsize(1);            
            plot( handles.Axes_Global_XY, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
        hold( handles.Axes_Global_XY, 'off' );    

        hold( handles.Axes_Global_XZ, 'on' );
            w = [curCellDisplaySize(1), curCellStats.BoundingBox(6)];
            ptCorner = curCellCentroid([1,3]) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(2), 1 ) = imsize(2);
            ptBBox( ptBBox(:,2) > imsize(3), 2 ) = imsize(3);
            plot( handles.Axes_Global_XZ, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
            plot( handles.Axes_Global_XZ, [1, size(imGlobalXZ,2)], [curCellSliceId, curCellSliceId], 'g-' );
        hold( handles.Axes_Global_XZ, 'off' );    

        hold( handles.Axes_Global_YZ, 'on' );
            w = [curCellStats.BoundingBox(6), curCellDisplaySize(1)];
            ptCorner = curCellCentroid([3,2]) - ceil(0.5 * w);
            ptBBox = repmat( ptCorner, 5, 1 ) + [ 0, 0; w(1), 0; w ; 0, w(2); 0, 0 ];
            ptBBox( ptBBox < 1 ) = 1;
            ptBBox( ptBBox(:,1) > imsize(3), 1 ) = imsize(3);
            ptBBox( ptBBox(:,2) > imsize(1), 2 ) = imsize(1);
            plot( handles.Axes_Global_YZ, ptBBox(:,1), ptBBox(:,2), 'r-', 'LineWidth', 2.0 ); 
            plot( handles.Axes_Global_YZ, [curCellSliceId, curCellSliceId], [1, size(imGlobalYZ,1)], 'g-' );
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
    subinds{3} = curCellSliceId;    
    
    if handles.data.flagUseLOG
        imCellCropped = mat2gray( handles.data.imageDataLOG{1}(subinds{:}), handles.data.imLogDisplayRange(1,:) );        
    else
        imCellCropped = mat2gray( handles.data.imageData{1}(subinds{:}), handles.data.imDisplayRange(1,:) );        
    end
    imCellSegCropped = handles.data.imLabelCellSeg( subinds{:} );
    imCellSegCropped = double( imCellSegCropped == handles.data.curCellId );    
    
    imCellSeedCropped = handles.data.imCellSeedPoints( subinds{:} );
    imCellSeedCropped(~imCellSegCropped) = 0;
    imCellSeedCropped(:, :, :, 2:3) = 0;
    
    % display the nuclei channel of the extracted cell image
    if handles.data.flagShowGlobalSegMask 
        
        curCellColor = handles.data.CellSegColorMap( handles.data.curCellId, : );
        imNucleusXYDisplay = genImageMaskOverlay( imCellCropped, {imCellSegCropped, imCellSeedCropped}, [curCellColor; 1 0 0], [0.2, 0.6] );
        
    else
        
        imNucleusXYDisplay = repmat( imCellCropped, [1,1,3] );
        
    end

    cla( handles.Axes_Nucleus_XY );
    image( imNucleusXYDisplay, 'Parent', handles.Axes_Nucleus_XY );
    set(  handles.Axes_Nucleus_XY, 'XTickLabel', [], 'YTickLabel', [] ); 
    
    % display red/green channel 
    if handles.data.flagRedGreenChannel && ~handles.data.flagShowLocalMIP
        
        flagShowChannel = [ handles.data.flagShowNucleusChannel, handles.data.flagShowGreenChannel, handles.data.flagShowRedChannel ];
        
        if ~any(flagShowChannel)
            imMultiChannelOverlay = zeros( [size(imCellCropped), 3] );
        else            
            imCellChannelData = [];
            displaycolors = [];
            for i = 1:3
               if flagShowChannel(i)
                   imCellChannelData = cat(3, imCellChannelData, mat2gray(handles.data.imageData{i}(subinds{:}), handles.data.imDisplayRange(i,:)) );
                   displaycolors = [ displaycolors; handles.data.metadata.channelColors(i,:) ];
               end
            end
            imMultiChannelOverlay = genMultiChannelOverlay(imCellChannelData, displaycolors);            
        end

        imAllChannelXYDisplay = imMultiChannelOverlay;
        
    else
        
        subinds{3} = round(curCellStats.BoundingBox(3):(curCellStats.BoundingBox(3)+curCellStats.BoundingBox(6)-1));

        imCurCellMIP = mat2gray( max( handles.data.imageData{1}( subinds{:} ), [], 3 ) );
        imCurCellSegMIP = max( double( handles.data.imLabelCellSeg( subinds{:} ) == handles.data.curCellId ), [], 3);
        imCurCellMIP(~imCurCellSegMIP) = 0;
        
        if handles.data.flagShowGlobalSegMask 

            curCellColor = handles.data.CellSegColorMap( handles.data.curCellId, : );
            imAllChannelXYDisplay = genImageMaskOverlay( imCurCellMIP, imCurCellSegMIP, curCellColor, 0.2 );
            
        else

            imAllChannelXYDisplay = repmat( imCurCellMIP, [1,1,3] );

        end        
        
    end
    
    cla( handles.Axes_AllChannel_XY );
    image( imAllChannelXYDisplay, 'Parent', handles.Axes_AllChannel_XY );
    set( handles.Axes_AllChannel_XY, 'XTickLabel', [], 'YTickLabel', [] );

    % if image objects were created for the first time -- store their handles
    % display cell id    
    strtmp = sprintf('%d / %d', handles.data.curCellId, numel(handles.data.cellStats) );
    set(handles.CellCountDisplay, 'String', strtmp);    
    
    % display cell slice id
    strtmp = sprintf('Slice: %d / %d', handles.data.curCellSliceId, imsize(3) );
    set(handles.EditboxSliceId, 'String', strtmp);    

    
% --------------------------------------------------------------------
function UpdateRedGreenDifferenceDisplay(handles)
    
    curCellStats = handles.data.cellStats( handles.data.curCellId );   
    
    prev_axis = gca;
    
    % show histograms of individual channels    
    axes( handles.Axes_Channel_Histogram );
    cla reset;    
    hold on;
    
    maxval = 4096;
    for chid = 1:3
        
        channelPixelIntensties{chid} = maxval * mat2gray(handles.data.imageData{chid}( curCellStats.PixelIdxList), ...
                                                       handles.data.imDisplayRange(chid,:));
                                                   
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

    curCellStats = handles.data.cellStats( handles.data.curCellId );
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
    if handles.data.flagRedGreenChannel
        
        curGreenPixelIntensities = mat2gray(handles.data.imageData{2}( curCellStats.PixelIdxList), handles.data.imDisplayRange(2,:));
        curRedPixelIntensities = mat2gray(handles.data.imageData{3}( curCellStats.PixelIdxList), handles.data.imDisplayRange(3,:));

        gr_ratio = mat2gray(curGreenPixelIntensities ./ (eps + curRedPixelIntensities));        
        strCellDescription = sprintf( '%s\n\nGreed-Red Median Intensity Ratio: %.2f', ...
                                      strCellDescription, median( gr_ratio ));                          

        ratioOfMedianGreenRedIntensity = median(curGreenPixelIntensities) / median(curRedPixelIntensities);  
        strCellDescription = sprintf( '%s\n\nMedian-Green Median-Red Ratio: %.2f', ...
                                      strCellDescription, ratioOfMedianGreenRedIntensity);                          

    end
    
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
    if isfield( handles.data.history, 'lastOutputDir' )
        [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastOutputDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    else
        [fileName,pathName] = uigetfile( fullfile( handles.data.history.lastAnalyzedNucleusDir, 'CellPatternAnnotation.mat' ), 'Select annotation file' );   
    end

    if ~fileName 
        return;
    end
    
    annotationFile = fullfile(pathName, fileName);
    
    % load data from file
    handles.data = load( annotationFile );
    
    set( handles.CellPatternAnnotator, 'Name', sprintf( 'Cell Pattern Annotator - %s', handles.data.dataFilePath{1} ) );
    
    if handles.data.flagRedGreenChannel
    
        set( handles.LabelAllChannels, 'String', 'All Channel Overlay' );
        set( handles.CheckboxShowNucleusChannel, 'Visible', 'on' );
        set( handles.CheckboxShowRedChannel, 'Visible', 'on' );
        set( handles.CheckboxShowGreenChannel, 'Visible', 'on' );            
        set( handles.CheckboxLocalMIP, 'Visible', 'on' );            
        
    else
        
        set( handles.LabelAllChannels, 'String', 'Nucleus Channel MIP' );
        set( handles.CheckboxShowNucleusChannel, 'Visible', 'off' );
        set( handles.CheckboxShowRedChannel, 'Visible', 'off' );
        set( handles.CheckboxShowGreenChannel, 'Visible', 'off' );
        set( handles.CheckboxLocalMIP, 'Visible', 'off' );      
        
    end
    
    if ~isfield( handles.data.metadata, 'channelColors' )
        handles.data.metadata.channelColors = [ 0 0 1; 0 1 0; 1 0 0 ];
    end
    
    % set cell pattern type list and cell class selector list
    set( handles.ListboxCellPatternSelector, 'String', handles.data.cellPatternTypes );
    set( handles.poplistCellClassSelector, 'String', handles.data.cellPatternTypes );
    
    % highlight pattern of current cell in listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);

    % Update Cell Descriptors
    UpdateCellDescriptors(handles);
    
% --------------------------------------------------------------------
function File_SaveAnnotation_Callback(~, eventdata, handles)
% hObject    handle to File_SaveAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
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
    if isfield( handles.data.history, 'lastOutputDir' )
        outputDir = uigetdir( handles.data.history.lastOutputDir, 'Select annotation output directory');
    else
        outputDir = uigetdir( handles.data.history.lastAnalyzedNucleusDir, 'Select annotation output directory');
    end

    if ~outputDir
        return;
    end

    handles.data.history.lastOutputDir = outputDir;
    
    % save data
    [pathstr, name, ext] = fileparts( handles.data.dataFilePath{1} );
    outputDir = strtrim( fullfile(outputDir, name) );

    if ~isdir( outputDir )
        mkdir( outputDir );
    end
    
    data = handles.data;
    
    h = waitbar(0, 'Saving Data ... Please Wait' );
    save( fullfile(outputDir, 'CellPatternAnnotation.mat'), '-struct', 'data' );
    close(h);
    
    clear data;

    % write summary file
    summary_fid = fopen( fullfile(outputDir, 'cellPatternAnnotationSummary.txt'), 'w' );    
    
        % print some information about the dataset
        fprintf( summary_fid, '\n>> Dataset Description:\n' );
        
        fprintf( summary_fid, '\n\tTwo-photon data file -- %s\n', handles.data.dataFilePath{1} );
    
        if handles.data.flagRedGreenChannel
            fprintf( summary_fid, '\n\tConfocal data file -- %s\n', handles.data.dataFilePath{2} );
        else            
            fprintf( summary_fid, '\n\tConfocal data file -- Not Specified\n' );
        end
        
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
                                          std( handles.data.imageData{j}(curClassPixelList) ) );
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
            curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.data.cellDisplaySize] );

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
            
            imCurCellCropped = mat2gray( handles.data.imageData{1}( subinds{1:2}, 1:imsize(3) ), handles.data.imDisplayRange(1,:) );
            imCurCellSegCropped = double( handles.data.imLabelCellSeg( subinds{1:2}, 1:imsize(3) ) == cellId );
            
            zROIMask = zeros( size(imCurCellCropped) );
            zROIMask(:, :, subinds{3}) = 1;            
            imCurCellCropped(~zROIMask) = 0;
            
            imCurCellMIP = imresize( genMIPImageForDisplay( imCurCellCropped, imCurCellSegCropped ), [400,400]);

            if handles.data.flagRedGreenChannel 
            
               imCurCellAllChannelMIP = [];
               
               for chid = 1:3
                   imCurChannelCropped = mat2gray( handles.data.imageData{chid}( subinds{1:2}, 1:imsize(3) ), handles.data.imDisplayRange(chid,:) );
                   imCurChannelCropped(~zROIMask) = 0;
                   imCurChannelMIP = imresize( genMIPImageForDisplay( imCurChannelCropped ), [400,400] );
                   imCurCellAllChannelMIP = cat( 3, imCurCellAllChannelMIP, imCurChannelMIP );
               end
               channelcolormap = [ 0, 0, 1; 0, 1, 0; 1, 0, 0 ];
               imCurCellAllChannelMIP = genMultiChannelOverlay( imCurCellAllChannelMIP, channelcolormap );

               imCurCellMIP(end:end+10,:,:) = 1;
               imCurCellMIP = cat(1, imCurCellMIP, imCurCellAllChannelMIP );
               
            end
            
            imwrite( imCurCellMIP, ...
                     fullfile(outputDir, curCellPatternType, sprintf('CellMIP_%.3d.png', cellId)), 'png' );   
                
            waitbar( cellId/numel(handles.data.cellStats), h);
        end
        
    end
    close(h);
    
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


% --- Executes when user attempts to file_close CellPatternAnnotator_Ver2.
function CellPatternAnnotator_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotator_Ver2 (see GCBO)
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
    if ~handles.data.flagDataLoaded 
        return;
    end
    
    % decrement cell id
    handles.data.curCellId = max(1, handles.data.curCellId - 1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));
    
    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
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
    if ~handles.data.flagDataLoaded 
        return;
    end

    % increment cell id
    handles.data.curCellId = min( numel(handles.data.cellStats), handles.data.curCellId + 1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));
    
    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
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
    if ~handles.data.flagDataLoaded 
        return;
    end

    % set cell id as last one
    handles.data.curCellId = numel( handles.data.cellStats );
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));
    
    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
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
    if ~handles.data.flagDataLoaded 
        return;
    end

    % decrement cell id
    handles.data.curCellId = 1;
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));

    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
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
    

% --------------------------------------------------------------------
function File_Set_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to File_Set_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.parameters = CellPatternAnnotatorParametersGUI_Ver1( handles.parameters );

    % Update handles structure
    guidata(hObject, handles);
    
    % Update Cell Visualization
    UpdateCellDisplay(handles);
    
% --- Executes during object deletion, before destroying properties.
function CellPatternAnnotator_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellPatternAnnotator_Ver2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % save history
    history = handles.data.history;
    [pathstr, name, ext] = fileparts( mfilename( 'fullpath' ) );
    historyFile = fullfile( pathstr, 'CellPatternAnnotatorHistory.mat' );
    save( historyFile, '-struct', 'history' );  

%     if ~isempty(handles.imarisApp) && isobject(handles.imarisApp)
%         delete( handles.imarisApp );
%     end       
    
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

% --------------------------------------------------------------------
function FnSliceScroll_Callback(hSrc, eventdata)

    handles = guidata(hSrc);
    
    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end
    
    imsize = size(handles.data.imageData{1});
    curCellSliceId = handles.data.curCellSliceId;

%     curCellId = handles.data.curCellId;
%     curCellMinSlice = round(curCellStats.BoundingBox(3));
%     curCellMaxSlice = round(curCellStats.BoundingBox(3) + curCellStats.BoundingBox(6));    
%     slicePadding = 1;    
%     if eventdata.VerticalScrollCount > 0
%         if curCellSliceId < min( imsize(3), curCellMaxSlice+slicePadding ); 
%             handles.data.curCellSliceId = handles.data.curCellSliceId + 1;
%         end
%     elseif eventdata.VerticalScrollCount < 0
%         if curCellSliceId > max(1, curCellMinSlice-slicePadding)
%             handles.data.curCellSliceId = handles.data.curCellSliceId - 1;
%         end
%     end

    if eventdata.VerticalScrollCount > 0
        if curCellSliceId < imsize(3) 
            handles.data.curCellSliceId = handles.data.curCellSliceId + 1;
        end
    elseif eventdata.VerticalScrollCount < 0
        if curCellSliceId > 1
            handles.data.curCellSliceId = handles.data.curCellSliceId - 1;
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
    handles.data.flagShowNucleusChannel = get(hObject,'Value');
    
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
    handles.data.flagShowGreenChannel = get(hObject,'Value');
    
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
    handles.data.flagShowRedChannel = get(hObject,'Value');
    
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
    handles.data.flagShowLocalMIP = get(hObject,'Value');
    
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
    if ~handles.data.flagDataLoaded || isempty(handles.data.unannotatedCellList) 
        return;
    end

    % set cell id to next unannotated cell
    handles.data.curCellId = handles.data.unannotatedCellList(1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));
    
    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);
    
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
    if ~handles.data.flagDataLoaded
        return;
    end

    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indPrevCellInCurClass = find( cellPatternVec(handles.data.curCellId-1:-1:1) == get(handles.poplistCellClassSelector,'Value') );
    
    if isempty( indPrevCellInCurClass )
        return;
    end
    
    indPrevCellInCurClass = handles.data.curCellId - indPrevCellInCurClass;
    
    % decrement cell id
    handles.data.curCellId = indPrevCellInCurClass(1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));

    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);

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
    if ~handles.data.flagDataLoaded
        return;
    end

    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indNextCellInCurClass = find( cellPatternVec(handles.data.curCellId+1:end) == get(handles.poplistCellClassSelector,'Value') );
        
    if isempty( indNextCellInCurClass )
        return;
    end
    
    indNextCellInCurClass = handles.data.curCellId + indNextCellInCurClass;
    
    % decrement cell id
    handles.data.curCellId = indNextCellInCurClass(1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));

    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);

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
    if ~handles.data.flagDataLoaded
        return;
    end
    
    cellPatternVec = [handles.data.cellStats.cellPatternId];
    
    indFirstCellInCurClass = find( cellPatternVec == get(hObject,'Value') );
    
    if isempty( indFirstCellInCurClass )
        return;
    end
    
    % decrement cell id
    handles.data.curCellId = indFirstCellInCurClass(1);
    handles.data.curCellSliceId = round(handles.data.cellStats(handles.data.curCellId).Centroid(3));

    % update cell pattern listbox
    set(handles.ListboxCellPatternSelector, 'Value', handles.data.cellStats(handles.data.curCellId).cellPatternId);

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

    if ~handles.data.flagDataLoaded
        return;
    end

    if ~isempty(handles.imarisApp) && isobject(handles.imarisApp)
        delete( handles.imarisApp );
    end       

    segMaskRGB = label2rgbND(handles.data.imLabelCellSeg);
    
    stats = regionprops( bwlabeln( handles.data.imCellSeedPoints ), 'Centroid' );
    cellSeedPointLocations = cat( 1, stats.Centroid );

    handles.imarisApp = Display3DDataAndResultsInImaris( handles.data.imageData{1}, handles.data.metadata.voxelSpacing, ...
                                                        'rgbMask', segMaskRGB, ...
                                                        'spotLocations', cellSeedPointLocations, ...
                                                        'spotRadius', 3);

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

    if ~handles.data.flagDataLoaded
        return;
    end
    
    if ~isempty(handles.imarisAppCellSegCropped) && isobject(handles.imarisAppCellSegCropped)
        delete( handles.imarisAppCellSegCropped );
    end       
    
    curCellStats = handles.data.cellStats( handles.data.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.data.cellDisplaySize] );

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
    imCellSegCropped = double( imCellSegCropped == handles.data.curCellId );    
    
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
    curCellDisplayRange = handles.data.imDisplayRange;
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

    if ~handles.data.flagDataLoaded
        return;
    end
    
    curCellStats = handles.data.cellStats( handles.data.curCellId );        
    curCellCentroid = round( curCellStats.Centroid );    
    
    curCellBoundingBox = curCellStats.BoundingBox;
    curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.data.cellDisplaySize] );

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

    if ~handles.data.flagDataLoaded
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
            curCellDisplaySize = max( [curCellBoundingBox(4:5), handles.data.cellDisplaySize] );

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
    close( hStatusDlg );
    
    % display cell pattern distribution in imaris
    handles.imarisAppCellSegCropped = DisplayMultichannel3DDataInImaris( imageDataPadded, ...
                                                                         'spacing', handles.data.metadata.voxelSpacing, ...
                                                                         'surfaceObjects', cellPatternSurfaceObjectList, ...
                                                                         'displayRanges', handles.data.imDisplayRange, ...
                                                                         'displayColors', handles.data.metadata.channelColors );
    
    % Update handles structure
    guidata(hObject, handles);
    
