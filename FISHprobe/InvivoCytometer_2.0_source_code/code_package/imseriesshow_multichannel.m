function [varargout] = imseriesshow_multichannel(im, varargin ) 
% [varargout] = imseriesshow_multichannel(im, varargin)
% A simple function to show series of 2D images from a volume with multiple channels. 
%
% You can also use it to display a single 2D image with multiple channels.
%
% Input:
% 
%     im - cell array of various channels of an image 
% 
%     Variable arguments (should be specified as param-value pairs)
% 
%       displayColors - allows you to specify the RGB color of each channel
% 
%         should be an (numChannels x 3) matrix where numChannels is the number
%         of channels provided. If you dont specify this, colors will be 
%         assigned internally.            
% 
%       spacing - pixel spacing of the grayscale image
% 
%       displayRanges - intensity display ranges of the each channel
% 
%         specify this if you want to the image channels to be displayed 
%         in a specific intensity range. All intensity values outside this 
%         range will be made equal to the nearest border in the specified
%         intensity range.                           
%    
% % 2D case
% im2d_channel1 = zeros(500,500);
% im2d_channel2 = zeros(500,500);
% im2d_channel3 = zeros(500,500);
% 
% im2d_channel1(100:200,100:200) = rand(101,101);
% im2d_channel2(200:300,200:300) = rand(101,101);
% im2d_channel3(300:400,300:400) = rand(101,101);
% 
% imseriesshow_multichannel( {im2d_channel1, im2d_channel2, im2d_channel3} ); % 2D multi-channel image overlay
% set( gcf, 'Name', '2D Example: Multi-channel imag display' );
%
% % 3D case
% im3d_channel1 = zeros(500,500,10);
% im3d_channel2 = zeros(500,500,10);
% im3d_channel3 = zeros(500,500,10);
% 
% im3d_channel1(100:200,100:200,:) = rand(101,101,10);
% im3d_channel1(200:300,200:300,:) = rand(101,101,10);
% im3d_channel1(300:400,300:400,:) = rand(101,101,10);
%
% imseriesshow_multichannel( {im3d_channel1, im3d_channel2, im3d_channel3} ); % 3D multi-channel image overlay
% set( gcf, 'Name', '3D Example: Multi-channel image display' ); 
%
% imseriesshow_multichannel( {im3d_channel1, im3d_channel2, im3d_channel3}, 'displayColors', [1,1,0; 0,1,1; 1,0,1] ); % 3D multi-channel image overlay
% set( gcf, 'Name', '3D Example: Multi-channel image display with user-specified colors' );
%
% Author: Deepak Roy Chittajallu 
%
%%

% get required parameters
p = inputParser;
p.addRequired( 'im', @(x) ( iscell(x) && ~isempty(x) ) );
p.parse(im);

im = p.Results.im;
volSize = size(im{1});
volSize(3) = size(im{1}, 3);
numChannels = numel(im);

default_displayranges = zeros(numChannels,2);
for i = 1:numChannels
   curChannelVolSize = size(im{i});
   curChannelVolSize(3) = size(im{i}, 3);
   if ~isempty( setxor( curChannelVolSize, volSize ) )
      error('Images corresponding to all channels must of the same size and dimension'); 
   end
   im{i} = double( im{i} );
   %default_displayranges(i,:) = double([ min(im{i}(:)) max(im{i}(:))]);
   default_displayranges(i,:) = ComputeImageDynamicRange( im{i}, 98.0 );
end

cMap(1,:) = [1 0 0];
cMap(2,:) = [0 1 0];
cMap(3,:) = [0 0 1];
cMap(4,:) = [1 0 1];
cMap(5,:) = [0 1 1];
cMap(6,:) = [1 1 0];

% get optional parameters
p.addParamValue( 'spacing', [1,1,1], @(x) ( isnumeric(x) && ~isscalar(x) && numel(x) == ndims(im{1}) ) );
p.addParamValue( 'displayRanges', default_displayranges, @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 2 && size(x,1) == numChannels ) );
p.addParamValue( 'displayColors', cMap(1:numChannels,:), @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 3 && size(x,1) == numChannels ) );
p.parse(im, varargin{:});

spacing = p.Results.spacing;
displayranges = p.Results.displayRanges;
displaycolors = p.Results.displayColors;

%% Data
data.im = im;
data.volSize = volSize;
data.numChannels = numChannels;
data.displayranges = displayranges;
data.displaycolors = displaycolors;

% Slice and planes
data.sliceno = 1;
data.curPlane = 3;
data.plane(1).slice = 1;
data.plane(2).slice = 1;
data.plane(3).slice = 1;

% Planes volume sizes
data.plane(1).volSize = [volSize(3) volSize(1) volSize(2)];
data.plane(2).volSize = [volSize(3) volSize(2) volSize(1)];
data.plane(3).volSize = [volSize(1) volSize(2) volSize(3)];

% Default x-limits and y-limits for each plane
data.plane(1).xlim = [ 0.5 volSize(1)+0.5 ];
data.plane(1).ylim = [ 0.5 volSize(3)+0.5 ];
data.plane(2).xlim = [ 0.5 volSize(2)+0.5 ];
data.plane(2).ylim = [ 0.5 volSize(3)+0.5 ];
data.plane(3).xlim = [ 0.5 volSize(2)+0.5 ];
data.plane(3).ylim = [ 0.5 volSize(1)+0.5 ];

% Data spacing
data.spacingUse = ~ismember( 'spacing', p.UsingDefaults );
data.plane(1).spacing = [spacing(3) spacing(1) 1];
data.plane(2).spacing = [spacing(3) spacing(2) 1];
data.plane(3).spacing = [spacing(2) spacing(1) 1];

% Data log
data.logUse = 0;
data.imLog = cell(1,numChannels);
data.imLogDisplayRanges = zeros( numChannels, 2 );
for i = 1:numChannels     
   imAdjusted = im{i} - min( im{i}(:) );
   ImageIntensityRange = ComputeImageDynamicRange( imAdjusted, 99.0 );
   log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0;
   imAdjusted = log( log_bottom + AdjustImageIntensityRange(imAdjusted, ImageIntensityRange) );
   data.imLogDisplayRanges(i,:) = [ min(imAdjusted(:)), max(imAdjusted(:)) ]; 
   data.imLog{i} = imAdjusted;
end

%% Create UI controls
hMainFigure = figure;

% Display figure toolbar
set(hMainFigure, 'ToolBar', 'figure')

% Default Axis
data.ui.hAxes = axes('Position', [0.0, 0.24, 1.0, 0.74], ...
                'YDir','reverse',...
                'TickDir', 'out', ...
                'XGrid', 'off', ...
                'YGrid', 'off', ...
                'DataAspectRatio', [1 1 1], ...
                'PlotBoxAspectRatio', [volSize(2) volSize(1) 1], ...
                'Visible', 'off');

% Default image
data.ui.hImage = image('BusyAction', 'cancel', ...
           'Parent', data.ui.hAxes, ...
           'Interruptible', 'off');

% Slice navigation controls
data.ui.pbh_dec = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '<<', ...
    'Units', 'normalized', 'Position', [0.3 0.11 0.066 0.05], ...
    'Callback', {@pushFirstSlice_Callback});

data.ui.pbh_dec = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '<', ...
    'Units', 'normalized', 'Position', [0.366 0.11 0.066 0.05], ...
    'Callback', {@pushdec_Callback});

data.ui.eth_sno = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', '0',...
    'Units', 'normalized', 'Position', [0.433 0.11 0.133 0.05]);

data.ui.pbh_inc = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '>', ...
    'Units', 'normalized', 'Position', [0.566 0.11 0.066 0.05], ...
    'Callback', {@pushinc_Callback});

data.ui.pbh_inc = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '>>', ...
    'Units', 'normalized', 'Position', [0.632 0.11 0.066 0.05], ...
    'Callback', {@pushLastSlice_Callback});

% Cursor point info controls
data.ui.eth_xloc = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'X: INV', ...
    'Units', 'normalized', 'Position', [0.3 0.06 0.133 0.05]);

data.ui.eth_yloc = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'Y: INV', ...
    'Units', 'normalized', 'Position', [0.433 0.06 0.133 0.05]);

data.ui.eth_Imval = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'I: INV', ...
    'Units', 'normalized', 'Position', [0.566 0.06 0.133 0.05]);

% Masks selection ui
numLimit = 10;
if numel(data.numChannels) >= numLimit
    sizeBox = 0.1;
else
    sizeBox = 0.06;
end
data.ui.ch_visibility = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.3 0.16 0.4 sizeBox]);
xmask = 5;
ymask = 5;
data.blnShowChannel = ones(1,data.numChannels);
for i = 1:data.numChannels
    if mod(i,numLimit) == 0
        xmask =  xmask + 20;
        ymask = 5;
    end
    % Bool show
    data.ui.channels(i) = uicontrol('Style','checkbox',...
        'Units', 'pixels', 'Position',[ymask xmask 15 15],...
        'parent',data.ui.ch_visibility,'HandleVisibility','off',...
        'Callback',{@checkChannels_Callback}, 'Value', 1, ...
        'BackgroundColor', data.displaycolors(i,:));
    ymask = 25 + ymask;
end

% Planes ui
data.ui.ch_planes = uibuttongroup('visible', 'on', 'Units', 'normalized', ...
    'Position', [0.3 0.01 0.4 0.05], ...
    'SelectionChangeFcn', {@planes_Callback});

if ndims(data.im{1}) == 3
    otherPlaneEnableState = 'on';
else
    otherPlaneEnableState = 'off';
end

data.ui.planeSP = uicontrol('Style', 'Radio', 'String', 'YZ',...
    'Units', 'normalized', 'Position', [0.05 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off','Enable',otherPlaneEnableState);
data.ui.planeCP = uicontrol('Style', 'Radio', 'String', 'XZ',...
    'Units', 'normalized', 'Position', [0.25 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off','Enable',otherPlaneEnableState);
data.ui.planeTP = uicontrol('Style', 'Radio', 'String', 'XY',...
    'Units', 'normalized', 'Position', [0.45 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', 1);
data.ui.planeLog = uicontrol('Style', 'checkbox', 'String', 'log',...
    'Units', 'normalized', 'Position', [0.65 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', data.logUse, ...
    'Callback',{@checkLog_Callback});
data.ui.planeSpacing = uicontrol('Style', 'checkbox', 'String', 'Spacing',...
    'Units', 'normalized', 'Position', [0.85 0.15 0.10 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', data.spacingUse, ...
    'Callback',{@checkSpacing_Callback});

% Set callbacks
set(hMainFigure, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)
set(hMainFigure, 'WindowButtonMotionFcn', @FnMainFig_MouseMotionFunc);

% Pan and Zoom callback
hZoom = zoom(hMainFigure);
set(hZoom, 'ActionPostCallback', @postPanZoom_Callback);
hPan = pan(hMainFigure);
set(hPan, 'ActionPostCallback', @postPanZoom_Callback);

% Define default output and return it if it is requested by users
mOutputArgs{1} = hMainFigure;
if nargout > 0
    [varargout{1:nargout}] = mOutputArgs{:};
end
guidata(hMainFigure,data);
imsliceshow(data);
% close(gcbf);

end

%% Show slice
function imsliceshow(data)

% Compute rgb image and display it
cdata = getRGBSlice(data);
set(data.ui.hImage, 'cdata', cdata);

% Slice number
volSize = data.plane(data.curPlane).volSize;
strtmp = sprintf('%d / %d', data.sliceno, volSize(3));
set(data.ui.eth_sno, 'String', strtmp);

end

%% Compute a rgb slice
function [ slice_rgb ] = getRGBSlice(data)

if ~any( data.blnShowChannel )
    slice_rgb = zeros( [ data.volSize(1:2), 3 ] );
    return;
end

if data.logUse      
    im = data.imLog;
    displayranges = data.imLogDisplayRanges;    
else
    im = data.im;
    displayranges = data.displayranges;
end

% Get the slice, depending on the plane showed
slice = [];
for i = 1:data.numChannels    
    if ~data.blnShowChannel(i)
        continue;
    end        
    if data.curPlane == 1
        curSlice = squeeze(im{i}(:,data.sliceno,:))';
    elseif data.curPlane == 2
        curSlice = squeeze(im{i}(data.sliceno,:,:))';
    else
        curSlice = squeeze(im{i}(:,:,data.sliceno));
    end    
    
    curSlice = mat2gray( curSlice, displayranges(i,:) );    
    slice = cat( 3, slice, curSlice );
end

slice_rgb = blend_multichannel_image( slice, data.displaycolors(data.blnShowChannel>0, :) );

end

%% Check/Uncheck channel visibility
function checkChannels_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Update value for each mask
for i = 1:data.numChannels
    data.blnShowChannel(i) = get(data.ui.channels(i), 'Value');
end

guidata(hSrc, data);
imsliceshow(data);

end


%% Get real coordinates, depending on the plane
function [ ind coord ] = getRealCoord(data, coord)

if data.curPlane == 1
    coord = [coord(2) coord(3) coord(1)];
elseif data.curPlane == 2
    coord = [coord(3) coord(2) coord(1)];
end
ind = sub2ind(data.volSize, coord(1), coord(2), coord(3));

end

%% Zoom post callback
function postPanZoom_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Save x-limits and y-limits
data.plane(data.curPlane).xlim = get(data.ui.hAxes, 'XLim');
data.plane(data.curPlane).ylim = get(data.ui.hAxes, 'YLim');

% Check the plot ratio is good
volSize = data.plane(data.curPlane).volSize;
set(data.ui.hAxes, 'PlotBoxAspectRatio', [volSize(2) volSize(1) 1]);

% Spacing
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

guidata(hSrc, data);

end

%% Check/Uncheck Spacing use
function checkSpacing_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Modify spacing
data.spacingUse = get(data.ui.planeSpacing, 'Value');
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Check/Uncheck Log use
function checkLog_Callback(hSrc, eventdata) %#ok<INUSD>

    data = guidata(hSrc);

    % Modify spacing
    data.logUse = get(data.ui.planeLog, 'Value');

    guidata(hSrc, data);
    imsliceshow(data);

end


%% Planes
function planes_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Save slice number, x-limits and y-limits
data.plane(data.curPlane).slice = data.sliceno;
data.plane(data.curPlane).xlim = get(data.ui.hAxes, 'XLim');
data.plane(data.curPlane).ylim = get(data.ui.hAxes, 'YLim');

% Change to new plane
if get(data.ui.planeSP, 'Value')        % Sagittal
    data.curPlane = 1;
elseif get(data.ui.planeCP, 'Value')    % Coronal
    data.curPlane = 2;
elseif get(data.ui.planeTP, 'Value')    % Tranverse
    data.curPlane = 3;
end
data.sliceno = data.plane(data.curPlane).slice;

% Plot ratio
volSize = data.plane(data.curPlane).volSize;
set(data.ui.hAxes, 'PlotBoxAspectRatio', [volSize(2) volSize(1) 1]);

% Set zoom
set(data.ui.hAxes, 'XLim', data.plane(data.curPlane).xlim);
set(data.ui.hAxes, 'YLim', data.plane(data.curPlane).ylim);

% Spacing
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

% Save and show slice
guidata(hSrc, data);
imsliceshow(data);

end

%% First Slice
function pushFirstSlice_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);
data.sliceno = 1;

guidata(hSrc, data);
imsliceshow(data);

end

%% Last Slice
function pushLastSlice_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);
volSize = data.plane(data.curPlane).volSize;
data.sliceno = volSize(3);

guidata(hSrc, data);
imsliceshow(data);

end

%% Dec Slice
function pushdec_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

if data.sliceno > 1
    data.sliceno = data.sliceno - 1;
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Inc Slice
function pushinc_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

volSize = data.plane(data.curPlane).volSize;
if data.sliceno < volSize(3)
    data.sliceno = data.sliceno+1;
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Slice scrolling
function FnSliceScroll_Callback(hSrc, eventdata)

data = guidata(hSrc);
volSize = data.plane(data.curPlane).volSize;

if eventdata.VerticalScrollCount > 0
    if data.sliceno < volSize(3)
        data.sliceno = data.sliceno + 1;
    end
elseif eventdata.VerticalScrollCount < 0
    if data.sliceno > 1 
        data.sliceno = data.sliceno - 1;
    end
end

guidata(hSrc, data);
imsliceshow(data);
UpdateCursorPointInfo(data);

end

%% Update cursor point info -- xloc, yloc, int_val
function UpdateCursorPointInfo(data)

cp = get(gca, 'CurrentPoint');
cp = round(cp(1, 1:2));
volSize = data.plane(data.curPlane).volSize;

% Display pointer coordinates and value
if IsPointInsideImage( cp(1,1:2) , data )
    set(data.ui.eth_xloc, 'String', sprintf('X: %d / %d', cp(1,1), volSize(2)));
    set(data.ui.eth_yloc, 'String', sprintf('Y: %d / %d', cp(1,2), volSize(1)));
    
    ind = getRealCoord(data, [cp(2) cp(1) data.sliceno]);
    set(data.ui.eth_Imval,'String',sprintf('I: %.1f' , data.im{1}(ind)));
else
    set(data.ui.eth_xloc,'String',sprintf('X: INV') );
    set(data.ui.eth_yloc,'String',sprintf('Y: INV') );
    set(data.ui.eth_Imval,'String',sprintf('I: INV') );
end

end

%% Figure pointer
function FnMainFig_MouseMotionFunc(hSrc, eventdata) %#ok<INUSD>

% No change if use zoom or pan
hZoom = zoom(hSrc);
hPan = pan(hSrc);
data = guidata(hSrc);

if strcmp(get(hZoom, 'Enable'), 'off') && strcmp(get(hPan, 'Enable'), 'off')
    cp = get(gca, 'CurrentPoint');
    if IsPointInsideImage(cp(1,1:2), data)
        set(hSrc, 'Pointer', 'crosshair');
    else
        set(hSrc, 'Pointer', 'arrow');
    end
end

% Cursor info
UpdateCursorPointInfo(data);

end

%%
function [ blnInside ] = IsPointInsideImage(cp, data)

% Point inside figure limits
volInfLim = ceil([ data.plane(data.curPlane).xlim(1) data.plane(data.curPlane).ylim(1) ]);
volSupLim = floor([ data.plane(data.curPlane).xlim(2) data.plane(data.curPlane).ylim(2) ]);
blnInside = all( cp <= volSupLim ) && all( cp >= volInfLim );

end

%%
function [ imCorrected ] = AdjustImageIntensityRange( im, ImageIntensityRange )
    
    imCorrected = mat2gray( im, ImageIntensityRange ) * range(ImageIntensityRange);
    
end

%%
function [ intensityRange ] = ComputeImageDynamicRange( im, cover_percent )

    [p,x] = hist( double(im), 255 );   
    p = p / sum(p);
    
    min_xlow = [];
    min_xhigh = [];
    min_xdiff = [];
    
    for i = 1:numel(x)
        for j = i+1:numel(x)
    
            if sum( p(i:j) ) < 0.01 * cover_percent
                continue;
            end
            
            if isempty(min_xdiff) || (x(j) - x(i)) < min_xdiff
                min_xlow = x(i);
                min_xhigh = x(j);
                min_xdiff = x(j) - x(i);              
            end
        end
    end
    
    w = 0.5 * (x(2) - x(1));
    intensityRange = [min_xlow-w, min_xhigh+w];
    
end