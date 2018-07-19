function varargout = LoGBlobDetectionGUI(varargin)
% LOGBLOBDETECTIONGUI MATLAB code for LoGBlobDetectionGUI.fig
%      LOGBLOBDETECTIONGUI, by itself, creates a new LOGBLOBDETECTIONGUI or raises the existing
%      singleton*.
%
%      H = LOGBLOBDETECTIONGUI returns the handle to a new LOGBLOBDETECTIONGUI or the handle to
%      the existing singleton*.
%
%      LOGBLOBDETECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOGBLOBDETECTIONGUI.M with the given input arguments.
%
%      LOGBLOBDETECTIONGUI('Property','Value',...) creates a new LOGBLOBDETECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoGBlobDetectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoGBlobDetectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoGBlobDetectionGUI

% Last Modified by GUIDE v2.5 04-Dec-2012 15:35:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoGBlobDetectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LoGBlobDetectionGUI_OutputFcn, ...
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


% --- Executes just before LoGBlobDetectionGUI is made visible.
function LoGBlobDetectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoGBlobDetectionGUI (see VARARGIN)

% Choose default command line output for LoGBlobDetectionGUI
handles.output = hObject;

% get input arguments
if isempty(varargin)

    fgMeanVar = [ 200, 20 ];
    bgMeanVar = [ 180, 20 ];
    blobRadius = 10 * [1, 4];

    im = zeros(200,200);
    fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
    bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

    imsize = size(im);
    [X,Y] = meshgrid(1:imsize(1),1:imsize(2));

    im(:) = random(bgGmObj,numel(im));
    xc = 0.5 * imsize(2);
    yc = 0.5 * imsize(1);
    pts = [ X(:) -  xc, Y(:) - yc ];
    curEllipseInd = find( (pts(:,1).^2 / (blobRadius(2))^2)  + (pts(:,2).^2 / (blobRadius(1))^2) - 1 <= 0 );    
    im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 

    handles.data.im = im;
    handles.data.spacing = ones(1,2);
    handles.data.blobDiameterRange = 2 * blobRadius;    
    handles.data.curBlobDiameter = mean(handles.data.blobDiameterRange);
    handles.data.flagBrightBlobs = true;
    handles.data.flagDataLoaded = true;
   
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
    p.parse( im, blobDiameterRange, varargin{3:end} );
    
    handles.data.im = p.Results.im;
    handles.data.spacing = p.Results.spacing;
    handles.data.flagBrightBlobs = p.Results.flagBrightBlobs;
    handles.data.blobDiameterRange = p.Results.blobDiameterRange;
    handles.data.curBlobDiameter = mean(handles.data.blobDiameterRange);
    handles.data.flagDataLoaded = true;
    
end

handles.data.curSliceId = 1;

set( handles.editCellDiameterMin, 'String',  num2str(handles.data.blobDiameterRange(1)) );
set( handles.editCellDiameterMax, 'String',  num2str(handles.data.blobDiameterRange(2)) );

set( handles.editCellDiameterValue, 'String', num2str(handles.data.curBlobDiameter) );
set( handles.editCellDiameterValue, 'String', num2str(handles.data.curBlobDiameter) );

set( handles.sliderCellDiameter, 'Min',  handles.data.blobDiameterRange(1) );
set( handles.sliderCellDiameter, 'Max',  handles.data.blobDiameterRange(2) );
set( handles.sliderCellDiameter, 'Value', handles.data.curBlobDiameter );

set(handles.figLoGBlobDetectionGUI, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)

% Run blob detection algorithm
handles = RunBlobDetectionAlgorithm( handles );

% Update handles structure
guidata(hObject, handles);

% Update display
UpdateDisplay( handles );

% UIWAIT makes LoGBlobDetectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figLoGBlobDetectionGUI);

% --- run blob detection algorithm
function [handles] = RunBlobDetectionAlgorithm(handles)

    [handles.data.imBlobSeedPoints, handles.data.imFilterResponse ] = detectBlobsUsingLoG( handles.data.im, ...
                                                                                           handles.data.curBlobDiameter, ...
                                                                                           'spacing', handles.data.spacing, ...
                                                                                           'minBlobDistance', 0.95 * min(handles.data.blobDiameterRange), ...
                                                                                           'flagBrightBlobs', handles.data.flagBrightBlobs );   
    
    
    % suppress all maxima with value less than 10% of the maximum value
    maxval = max( handles.data.imBlobSeedPoints(:) );
    handles.data.imBlobSeedPoints( handles.data.imBlobSeedPoints <= 0.1 * maxval ) = 0;
    
    % make an image with a sphere around each 
    [imDistMap, imClosestInd] = bwdist( handles.data.imBlobSeedPoints > 0 );
    handles.data.imBlobSeedScale = zeros( size(handles.data.im) );
    flagIsPixelInSeedVicinity = imDistMap < 0.5 * handles.data.curBlobDiameter;
    handles.data.imBlobSeedScale( flagIsPixelInSeedVicinity ) = imClosestInd( flagIsPixelInSeedVicinity );    
    handles.data.imBlobSeedScaleRGB = label2rgbND( handles.data.imBlobSeedScale );
    
    handles.data.imBlobSeedPointsDilated = imdilate( double(handles.data.imBlobSeedPoints > 0), ones(3,3) );    
    
% --- Display
function UpdateDisplay( handles )

    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end       
    
    curSliceId = handles.data.curSliceId;
    im = handles.data.im(:, :, curSliceId);    
    imBlobSeedPointsDilated = handles.data.imBlobSeedPointsDilated(:, :, curSliceId);
    imFilterResponse = handles.data.imFilterResponse(:, :, curSliceId);
    
    switch ndims( handles.data.im )
        
        case 2

            imBlobSeedScaleRGB = handles.data.imBlobSeedScaleRGB;
            imBlobSeedPointsDilatedRGB = imBlobSeedPointsDilated;
            imBlobSeedPointsDilatedRGB(:, :, 2:3) = 0;
            
            % display input image
            cla( handles.hAxes_InputImage );
            image( genImageRGBMaskOverlay( im, {imBlobSeedScaleRGB, imBlobSeedPointsDilatedRGB} ), ...
                   'Parent', handles.hAxes_InputImage );
            set( handles.hAxes_InputImage, 'XTickLabel', [], 'YTickLabel', [] );

            % display filter response
            cla( handles.hAxes_FilterResponse );
            image( genImageRGBMaskOverlay( imFilterResponse, {imBlobSeedScaleRGB, imBlobSeedPointsDilatedRGB} ), ...
                   'Parent', handles.hAxes_FilterResponse );
            colormap( 'jet' );
            set( handles.hAxes_FilterResponse, 'XTickLabel', [], 'YTickLabel', [] );

            % display circle showing scale
            [cy, cx] = find( handles.data.imBlobSeedPoints > 0 );                

            theta = 0:0.1:(2*pi+0.1);
            cx = cx(:,ones(size(theta)));
            cy = cy(:,ones(size(theta)));
            blobRadius = (0.5 * handles.data.curBlobDiameter);
            rad = blobRadius * ones( size(cx) );
            theta = theta(ones(size(cx,1),1),:);                
            X = cx + cos(theta).* rad;
            Y = cy + sin(theta).* rad;

            axes( handles.hAxes_InputImage );
            hold on;
                line(X', Y', 'Color', 'b', 'LineWidth', 1.5);                
            hold off;

            axes( handles.hAxes_FilterResponse );
            hold on;
                line(X', Y', 'Color', 'b', 'LineWidth', 1.5);                
            hold off;
            
        case 3

            imBlobSeedScaleRGB = squeeze(handles.data.imBlobSeedScaleRGB(:,:,curSliceId,:));
            imBlobSeedPointsDilatedRGB = imBlobSeedPointsDilated;
            imBlobSeedPointsDilatedRGB(:, :, :, 2:3) = 0;
            
            % display input image
            cla( handles.hAxes_InputImage );
            image( genImageRGBMaskOverlay( im, {imBlobSeedScaleRGB, imBlobSeedPointsDilatedRGB} ), ...
                   'Parent', handles.hAxes_InputImage );
            set( handles.hAxes_InputImage, 'XTickLabel', [], 'YTickLabel', [] );

            % display filter response
            cla( handles.hAxes_FilterResponse );
            image( genImageRGBMaskOverlay( imFilterResponse, {imBlobSeedScaleRGB, imBlobSeedPointsDilatedRGB} ), ...
                   'Parent', handles.hAxes_FilterResponse );
            colormap( 'jet' );
            set( handles.hAxes_FilterResponse, 'XTickLabel', [], 'YTickLabel', [] );
            
    end
    
% --------------------------------------------------------------------
function FnSliceScroll_Callback(hSrc, eventdata)

    handles = guidata(hSrc);
    
    % first check if data has been loaded
    if ~handles.data.flagDataLoaded 
        return;
    end        
    
    if ndims(handles.data.im)< 3
        return;
    end
    
    curSliceId = handles.data.curSliceId;
    
    if eventdata.VerticalScrollCount > 0
        if curSliceId < size(handles.data.im,3) 
            handles.data.curSliceId = handles.data.curSliceId + 1;
        end
    elseif eventdata.VerticalScrollCount < 0
        if curSliceId > 1
            handles.data.curSliceId = handles.data.curSliceId - 1;
        end
    end
    
    guidata(hSrc, handles);    
    UpdateDisplay(handles);
    

% --- Outputs from this function are returned to the command line.
function varargout = LoGBlobDetectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderCellDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.data.curBlobDiameter = get(hObject,'Value');
set( handles.editCellDiameterValue, 'String', num2str(handles.data.curBlobDiameter) );

% Run blob detection algorithm
handles = RunBlobDetectionAlgorithm( handles );

% Update handles structure
guidata(hObject, handles);

% Update display
UpdateDisplay( handles );


% --- Executes during object creation, after setting all properties.
function sliderCellDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editCellDiameterMin_Callback(hObject, eventdata, handles)
% hObject    handle to editCellDiameterMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCellDiameterMin as text
%        str2double(get(hObject,'String')) returns contents of editCellDiameterMin as a double


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



function editCellDiameterValue_Callback(hObject, eventdata, handles)
% hObject    handle to editCellDiameterValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCellDiameterValue as text
%        str2double(get(hObject,'String')) returns contents of editCellDiameterValue as a double


% --- Executes during object creation, after setting all properties.
function editCellDiameterValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCellDiameterValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
