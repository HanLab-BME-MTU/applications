function varargout = BlobDetectionGUI(varargin)
% BLOBDETECTIONGUI MATLAB code for BlobDetectionGUI.fig
%      BLOBDETECTIONGUI, by itself, creates a new BLOBDETECTIONGUI or raises the existing
%      singleton*.
%
%      H = BLOBDETECTIONGUI returns the handle to a new BLOBDETECTIONGUI or the handle to
%      the existing singleton*.
%
%      BLOBDETECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLOBDETECTIONGUI.M with the given input arguments.
%
%      BLOBDETECTIONGUI('Property','Value',...) creates a new BLOBDETECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BlobDetectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BlobDetectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BlobDetectionGUI

% Last Modified by GUIDE v2.5 06-Dec-2012 16:22:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BlobDetectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BlobDetectionGUI_OutputFcn, ...
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


% --- Executes just before BlobDetectionGUI is made visible.
function BlobDetectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BlobDetectionGUI (see VARARGIN)

% Choose default command line output for BlobDetectionGUI
handles.output = hObject;

% initialize some variables and ui controls
handles.flagShowBlobSizeMask = true;
handles.flagShowBlobSeedMask = true;
handles.curSliceId = 1;

set(handles.figBlobDetectionGUI, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)
set(handles.checkboxShowBlobSizeMask, 'Value', true );
set(handles.checkboxShowBlobSeedMask, 'Value', true );

% set list of detection methods
handles.blobDetectionMethodList = [];

    handles.blobDetectionMethodList(end+1).name = 'FixedScaleLoG';
    handles.blobDetectionMethodList(end).func = @detectBlobsUsingLoG;

    handles.blobDetectionMethodList(end+1).name = 'MultiscaleLoG';
    handles.blobDetectionMethodList(end).func = @detectBlobsUsingMultiscaleLoG;
    
    handles.blobDetectionMethodList(end+1).name = 'AdaptiveMultiscaleLoG';
    handles.blobDetectionMethodList(end).func = @detectBlobsUsingAdaptiveMultiscaleLoG;

    handles.blobDetectionMethodList(end+1).name = 'MultiscaleLoBG';
    handles.blobDetectionMethodList(end).func = @detectBlobsUsingMultiscaleLoBG;

set( handles.poplistBlobDetectionMethod, 'String', { handles.blobDetectionMethodList.name } );    
set( handles.poplistBlobDetectionMethod, 'Value', 1 );    

handles.curBlobDetectionMethodId = get(handles.poplistBlobDetectionMethod,'Value');
handles.curBlobDetectionMethodName = handles.blobDetectionMethodList( handles.curBlobDetectionMethodId ).name;

% get input arguments
if isempty(varargin)

    fgMeanVar = [ 200, 20 ];
    bgMeanVar = [ 180, 20 ];
    blobRadius = 10 * [1, 4];

    im = zeros(200,200);
    fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
    bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

    imsize = size(im);
    foregroundMask = zeros( size(im) );
    [X,Y] = meshgrid(1:imsize(1),1:imsize(2));

    im(:) = random(bgGmObj,numel(im));
    xc = 0.5 * imsize(2);
    yc = 0.5 * imsize(1);
    pts = [ X(:) -  xc, Y(:) - yc ];
    curEllipseInd = find( (pts(:,1).^2 / (blobRadius(2))^2)  + (pts(:,2).^2 / (blobRadius(1))^2) - 1 <= 0 );    
    im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
    foregroundMask( curEllipseInd ) = 1;
    
    handles.data.im = im;
    handles.data.foregroundMask = foregroundMask;
    handles.data.spacing = ones(1,2);
    handles.data.blobDiameterRange = 2 * blobRadius;    
    handles.data.flagBrightBlobs = true;
    handles.flagDataLoaded = true;
   
else

    if numel( varargin ) < 2
        fprintf( '\nUsage: LoGBlobDetectionGUI( im, blobDiameterRange )' );
       error( 'ERROR: improper function call' ); 
    end
    
    im = varargin{1};
    blobDiameterRange = varargin{2};
    
    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'blobDiameterRange', @(x) (numel(x) == 2) );        
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'flagBrightBlobs', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'foregroundMask', [], @(x) (ndims(x) == ndims(im)) );
    p.parse( im, blobDiameterRange, varargin{3:end} );
    
    handles.data.im = p.Results.im;
    handles.data.spacing = p.Results.spacing;
    handles.data.flagBrightBlobs = p.Results.flagBrightBlobs;
    handles.data.blobDiameterRange = p.Results.blobDiameterRange;
    handles.flagDataLoaded = true;
    
    if isempty( p.Results.foregroundMask )        
        handles.data.foregroundMask = double(handles.data.im > thresholdOtsu( im ));
    else
        handles.data.foregroundMask = p.Results.foregroundMask;
    end
    
end

handles.data.displayRange = ComputeImageDynamicRange( handles.data.im, 99.0 );

% compute distance map
handles.data.imForegroundDistMap = bwdistsc( ~handles.data.foregroundMask, handles.data.spacing );

% set blob diameter range uicontrols
set( handles.editCellDiameterMin, 'String',  num2str(handles.data.blobDiameterRange(1)) );
set( handles.editCellDiameterMax, 'String',  num2str(handles.data.blobDiameterRange(2)) );

% Run blob detection algorithm
handles = RunBlobDetectionAlgorithm( handles );

% Update handles structure
guidata(hObject, handles);

% Update display
UpdateDisplay( handles );

% UIWAIT makes BlobDetectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figBlobDetectionGUI);

% --- run blob detection algorithm
function [handles] = RunBlobDetectionAlgorithm(handles)

    if strcmp( handles.curBlobDetectionMethodName, 'FixedScaleLoG' )

        [handles.data.imBlobSeedPoints, handles.data.imFilterResponse ] = detectBlobsUsingLoG( handles.data.im, ...
                                                                                               mean( handles.data.blobDiameterRange ), ...
                                                                                               'spacing', handles.data.spacing, ...
                                                                                               'minBlobDistance', 0.95 * handles.data.blobDiameterRange(1), ...
                                                                                               'flagBrightBlobs', handles.data.flagBrightBlobs );   

        handles.data.imBlobRadius = 0.5 * mean( handles.data.blobDiameterRange ) * ones( size(handles.data.im) );
        
    elseif strcmp( handles.curBlobDetectionMethodName, 'AdaptiveMultiscaleLoG' )

        [handles.data.imBlobSeedPoints, ...
         handles.data.imFilterResponse, ...
         handles.data.imBlobRadius ] = detectBlobsUsingAdaptiveMultiscaleLoG( handles.data.im, ...
                                                                              handles.data.imForegroundDistMap, ...
                                                                              'blobDiameterRange', handles.data.blobDiameterRange, ...
                                                                              'spacing', handles.data.spacing, ...
                                                                              'minBlobDistance', 0.95 * min(handles.data.blobDiameterRange), ...
                                                                              'flagBrightBlobs', handles.data.flagBrightBlobs );   
        
    else
        
        curBlobDetectionMethodId = handles.curBlobDetectionMethodId;
        curBlobDetectionMethod = handles.blobDetectionMethodList( curBlobDetectionMethodId ).func;
        [handles.data.imBlobSeedPoints, ...
         handles.data.imFilterResponse, ...
         handles.data.imBlobRadius ] = curBlobDetectionMethod( handles.data.im, ...
                                                             handles.data.blobDiameterRange, ...
                                                             'spacing', handles.data.spacing, ...
                                                             'minBlobDistance', 0.95 * min(handles.data.blobDiameterRange), ...
                                                             'flagBrightBlobs', handles.data.flagBrightBlobs );   

    end
    
    handles.data.displayRangeFilterResponse = ComputeImageIntensityRange( handles.data.imFilterResponse );
    
    % suppress all maxima with value less than 10% of the maximum value
    maxval = max( handles.data.imBlobSeedPoints(:) );
    handles.data.imBlobSeedPoints( handles.data.imBlobSeedPoints <= 0.1 * maxval ) = 0;
    
    % make a mask showing the size of each blob
    seedPixelInd = find( handles.data.imBlobSeedPoints > 0 );
    handles.data.numSeedPoints = numel( seedPixelInd );
    
    if handles.data.numSeedPoints > 0 
        
        seedPos = cell(1,ndims(handles.data.imBlobSeedPoints));
        [seedPos{:}] = ind2sub( size(handles.data.imBlobSeedPoints), seedPixelInd );
        seedPos = cat( 2, seedPos{[2,1]}, seedPos{3:end} );

        kd = KDTreeSearcher( seedPos * diag(handles.data.spacing) );
        
        pixelPos = cell(1, ndims(handles.data.imBlobSeedPoints));
        [pixelPos{:}] = ind2sub( size(handles.data.imBlobSeedPoints), (1:numel(handles.data.imBlobSeedPoints))' );
        pixelPos = cat( 2, pixelPos{[2,1]}, pixelPos{3:end} );     

        [closestSeedInd,distanceToSeed] = knnsearch(kd, pixelPos * diag(handles.data.spacing));    
        closestSeedPixelInd = seedPixelInd( closestSeedInd );

        handles.data.imBlobMask = zeros( size(handles.data.im) );
        flagIsPixelInSeedVicinity = distanceToSeed <= handles.data.imBlobRadius(closestSeedPixelInd);
        handles.data.imBlobMask( flagIsPixelInSeedVicinity ) = closestSeedInd( flagIsPixelInSeedVicinity );
        handles.data.imBlobMaskRGB = label2rgbND( handles.data.imBlobMask );    
        handles.data.imBlobSeedPointsDilated = imdilate( double(handles.data.imBlobSeedPoints > 0), ones(3,3) );    
        
    end

% ---- 
function [ intensityRange ] = ComputeImageIntensityRange( im )

    intensityRange = [min(im(:)) max(im(:))];
    
% --- Display
function UpdateDisplay( handles )

    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end       
    
    curSliceId = handles.curSliceId;
    im = mat2gray( handles.data.im(:, :, curSliceId), handles.data.displayRange );    
    imFilterResponse = mat2gray( handles.data.imFilterResponse(:, :, curSliceId), handles.data.displayRangeFilterResponse );
    
    if handles.data.numSeedPoints == 0 || (~handles.flagShowBlobSizeMask && ~handles.flagShowBlobSeedMask)

        % display input image
        cla( handles.hAxes_InputImage );
        image( repmat( im, [1,1,3] ), 'Parent', handles.hAxes_InputImage );
        set( handles.hAxes_InputImage, 'XTickLabel', [], 'YTickLabel', [] );

        % display filter response
        cla( handles.hAxes_FilterResponse );
        image( repmat( imFilterResponse, [1,1,3] ), 'Parent', handles.hAxes_FilterResponse );
        set( handles.hAxes_FilterResponse, 'XTickLabel', [], 'YTickLabel', [] );

    else
        
        masks = {};
        maskAlphas = [];
    
        imBlobSeedPointsDilated = handles.data.imBlobSeedPointsDilated(:, :, curSliceId);

        if ndims( handles.data.im ) == 2

            imBlobMaskRGB = handles.data.imBlobMaskRGB;
            imBlobSeedPointsDilatedRGB = imBlobSeedPointsDilated;
            imBlobSeedPointsDilatedRGB(:, :, 2:3) = 0;

        else

            imBlobMaskRGB = squeeze(handles.data.imBlobMaskRGB(:,:,curSliceId,:));
            imBlobSeedPointsDilatedRGB = imBlobSeedPointsDilated;
            imBlobSeedPointsDilatedRGB(:, :, :, 2:3) = 0;

        end
    
        
        if handles.flagShowBlobSizeMask 
            masks{end+1} = imBlobMaskRGB;
            maskAlphas(end+1) = 0.5;
        end

        if handles.flagShowBlobSeedMask 
            masks{end+1} = imBlobSeedPointsDilatedRGB;
            maskAlphas(end+1) = 0.5;
        end
        
        % display input image
        cla( handles.hAxes_InputImage );
        image( genImageRGBMaskOverlay( im, masks, maskAlphas ), ...
               'Parent', handles.hAxes_InputImage );
        set( handles.hAxes_InputImage, 'XTickLabel', [], 'YTickLabel', [] );

        % display filter response
        cla( handles.hAxes_FilterResponse );
        image( genImageRGBMaskOverlay( imFilterResponse, masks, maskAlphas ), ...
               'Parent', handles.hAxes_FilterResponse );
        set( handles.hAxes_FilterResponse, 'XTickLabel', [], 'YTickLabel', [] );
        
    end    
    
    if ndims( handles.data.im ) == 2
        
        % display circle showing scale
        if handles.data.numSeedPoints > 0 && handles.flagShowBlobSeedMask 

            seedInd = find( handles.data.imBlobSeedPoints > 0 );                
            [cy, cx] = ind2sub( size(handles.data.imBlobSeedPoints), seedInd );                

            theta = 0:0.1:(2*pi+0.1);
            cx = cx(:,ones(size(theta)));
            cy = cy(:,ones(size(theta)));
            blobRadius = (0.5 * handles.data.imBlobRadius(seedInd));
            rad = blobRadius(:,ones(size(theta)));
            theta = theta(ones(size(cx,1),1),:);                
            X = cx + cos(theta).* rad;
            Y = cy + sin(theta).* rad;

            X = X / handles.data.spacing(1);
            Y = Y / handles.data.spacing(2);

            axes( handles.hAxes_InputImage );
            hold on;
                line(X', Y', 'Color', 'b', 'LineWidth', 1.5);                
            hold off;

            axes( handles.hAxes_FilterResponse );
            hold on;
                line(X', Y', 'Color', 'b', 'LineWidth', 1.5);                
            hold off;

        end
            
    end
    
% --------------------------------------------------------------------
function FnSliceScroll_Callback(hSrc, eventdata)

    handles = guidata(hSrc);
    
    % first check if data has been loaded
    if ~handles.flagDataLoaded 
        return;
    end        
    
    if ndims(handles.data.im)< 3
        return;
    end
    
    curSliceId = handles.curSliceId;
    
    if eventdata.VerticalScrollCount > 0
        if curSliceId < size(handles.data.im,3) 
            handles.curSliceId = handles.curSliceId + 1;
        end
    elseif eventdata.VerticalScrollCount < 0
        if curSliceId > 1
            handles.curSliceId = handles.curSliceId - 1;
        end
    end
    
    guidata(hSrc, handles);    
    UpdateDisplay(handles);
    

% --- Outputs from this function are returned to the command line.
function varargout = BlobDetectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function editCellDiameterMin_Callback(hObject, eventdata, handles)
% hObject    handle to editCellDiameterMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCellDiameterMin as text
%        str2double(get(hObject,'String')) returns contents of editCellDiameterMin as a double

    strContents = get(hObject,'String');
    
    if isnan( str2double(strContents) ) || str2double(strContents) > handles.data.blobDiameterRange(2)
        set( handles.editCellDiameterMin, 'String',  num2str(handles.data.blobDiameterRange(1)) );
        errordlg( 'minimum cell diameter should be a real number and its value should be less than the maximum cell diameter' );        
        return;
    end
        
    handles.data.blobDiameterRange(1) = str2double( strContents );
    
    if ~handles.flagDataLoaded
        return;
    end        
        
    % Run blob detection algorithm
    handles = RunBlobDetectionAlgorithm( handles );

    % Update handles structure
    guidata(hObject, handles);

    % Update display
    UpdateDisplay( handles );
    

% --- Executes during object creation, after setting all properties.
function editCellDiameterMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCellDiameterMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCellDiameterMax_Callback(hObject, eventdata, handles)
% hObject    handle to editCellDiameterMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCellDiameterMax as text
%        str2double(get(hObject,'String')) returns contents of editCellDiameterMax as a double

    strContents = get(hObject,'String');
    
    if isnan( str2double(strContents) ) || str2double(strContents) < handles.data.blobDiameterRange(1)
        set( handles.editCellDiameterMax, 'String',  num2str(handles.data.blobDiameterRange(2)) );
        errordlg( 'maximum cell diameter should be a real number and its value should be greater than the minimum cell diameter' );        
        return;
    end
        
    handles.data.blobDiameterRange(2) = str2double( strContents );
    
    if ~handles.flagDataLoaded
        return;
    end        
        
    % Run blob detection algorithm
    handles = RunBlobDetectionAlgorithm( handles );

    % Update handles structure
    guidata(hObject, handles);

    % Update display
    UpdateDisplay( handles );

% --- Executes during object creation, after setting all properties.
function editCellDiameterMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCellDiameterMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in poplistBlobDetectionMethod.
function poplistBlobDetectionMethod_Callback(hObject, eventdata, handles)
% hObject    handle to poplistBlobDetectionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplistBlobDetectionMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplistBlobDetectionMethod

    handles.curBlobDetectionMethodId = get(hObject,'Value');
    handles.curBlobDetectionMethodName = handles.blobDetectionMethodList( handles.curBlobDetectionMethodId ).name;
    
    % Run blob detection algorithm
    handles = RunBlobDetectionAlgorithm( handles );

    % Update handles structure
    guidata(hObject, handles);

    % Update display
    UpdateDisplay( handles );
    

% --- Executes during object creation, after setting all properties.
function poplistBlobDetectionMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplistBlobDetectionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figBlobDetectionGUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figBlobDetectionGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkboxShowBlobSizeMask.
function checkboxShowBlobSizeMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowBlobSizeMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkboxShowBlobSizeMask

    handles.flagShowBlobSizeMask = get(hObject,'Value');
    
    guidata(hObject, handles);    
    UpdateDisplay(handles);
    

% --- Executes on button press in checkboxShowBlobSeedMask.
function checkboxShowBlobSeedMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowBlobSeedMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkboxShowBlobSeedMask

    handles.flagShowBlobSeedMask = get(hObject,'Value');
    
    guidata(hObject, handles);    
    UpdateDisplay(handles);
