function [varargout] = imseriesmaskshow(im, inmasks, varargin)
%% [varargout] = imseriesmaskshow(im, inmasks, varargin)
% A simple function to show series of 2D images from a volume and the associated binary segmentation masks.
% You can also use this function to overlay a single 2D image with one or more segmentation masks
%
% Input:
% 
%     im - input grayscale image 
% 
%     inmasks - masks that you want to overlay. 
% 
%       If you have more than one mask you can provide all of them in the 
%       form of a cell array or stack them on top of each other
% 
%     Variable arguments (should be specified as param-value pairs)
% 
%       maskColors - allows you to specify the RGB color of each mask
% 
%         should be an (numMasks x 3) matrix where numMasks is the number
%         of masks provided. If you dont specify this colors will be 
%         assigned internally.            
% 
%       maskAlphas - allows you to specify the transparency of each mask
% 
%         should be a vector of values in range [0,1]. The lenght of the 
%         vector should be equal to the number of masks provided.
%         Value of 0 makes the mask entirely transparent and a value of 1 
%         makes it totally opaque. 
% 
%       spacing - pixel spacing of the grayscale image
% 
%       displayRange - intensity display range of the input grayscale image
% 
%         specify this if you want to the grayscale image
%         to be displayed in a specific intensity range. 
%         All intensity values outside this range will be     
%         made equal to the nearest border in the specified
%         intensity range.                           
%                            
% Examples:
%
% % 2D case
% im = rand(500,500);    
% mask2d = zeros(size(im));
% mask2d(100:400,100:400) = 1;
% 
% mask2d_1 = zeros(size(im));
% mask2d_2 = zeros(size(im));
% mask2d_3 = zeros(size(im));
% 
% mask2d_1(100:200,100:200,:) = 1;
% mask2d_2(200:300,200:300,:) = 1;
% mask2d_3(300:400,300:400,:) = 1;
% 
% imseriesmaskshow(im, mask2d); % single mask overlay
% set( gcf, 'Name', '2D Example: Overlay of Single Mask' );
%
% imseriesmaskshow(im, { mask2d_1, mask2d_2, mask2d_3} ); % multi-mask overlay
% set( gcf, 'Name', '2D Example: Overlay of Multiple Masks provided as a cell array' );
% 
% imseriesmaskshow(im, cat( 3, mask2d_1, mask2d_2, mask2d_3 ) );   
% set( gcf, 'Name', '2D Example: Overlay of Multiple Masks provided as a stack of masks' );
%
% 
% % 3D case
% im = rand(500,500,10);    
% mask = zeros(size(im));
% mask(100:400,100:400,:) = 1;
% 
% mask1 = zeros(size(im));
% mask2 = zeros(size(im));
% mask3 = zeros(size(im));
% 
% mask1(100:200,100:200,:) = 1;
% mask2(200:300,200:300,:) = 1;
% mask3(300:400,300:400,:) = 1;
% 
% imseriesmaskshow(im, mask); % single mask overlay
% set( gcf, 'Name', '3D Example: Overlay of Single Mask' );
%
% imseriesmaskshow(im, { mask1, mask2, mask3}, 'maskColors', [1,0,0; 0,1,0; 0,0,1] ); % multi-mask overlay
% set( gcf, 'Name', '3D Example: Overlay of Multiple Masks provided as a cell array' );
%
% imseriesmaskshow(im, cat( 4, mask1, mask2, mask3 ), 'maskAlphas', [0.2,0.2,0.2] ); % multi-mask overlay  
% set( gcf, 'Name', '3D Example: Overlay of Multiple Masks provided as a stack of masks' );
%
% Author: Deepak Roy Chittajallu 
%
%%

% Get and Validate Required Parameters
p = inputParser;
p.addRequired( 'im', @(x)( isnumeric(x) && ~isscalar(x) && ismember( ndims(x), [2,3] ) ) );
p.addRequired( 'inmasks', @(x) ( ~isempty(x) && ( iscell(x) || ( (islogical(x) || isnumeric(x)) && ismember(ndims(x), [2,3,4]) ) ) ) );
p.parse( im, inmasks );

im = double(p.Results.im);
volSize = size(im);
volSize(3) = size(im,3);

if iscell(inmasks)
    for i = 1:numel( inmasks )
        if ndims(inmasks{i}) ~= ndims(im) || any( size(inmasks{i}) ~= size(im) )
            error('ERROR: each mask in the input cell array must be of the same size and dimension as that of the input image');
        end
        masks(i).im = inmasks{i};
    end	
else     
    inmasksize = size(inmasks);
    if any( volSize(1:ndims(im)) ~= inmasksize(1:ndims(im)) ) && ~ismember( (ndims(inmasksize) - ndims(im)), [0,1] )
        error('ERROR: sizes of the masks and the input image must match. check example usages of masks in the documentation');
    end
    switch ndims(im)
        
        case 2
            
            if ndims(inmasks) == 2 % just one 2D mask given
                masks(1).im = inmasks; 
            else ndims(inmasks) == 3 % user gave a stack of 2D masks in the form of a 3D matrix
                for i = 1:size(inmasks,3)
                    masks(i).im = inmasks(:,:,i);
                end
            end

        case 3
            
            if ndims(inmasks) == 3 % just one 3D mask given
                masks(1).im = inmasks; 
            else ndims(inmasks) == 4 % user gave a stack of 3D masks in the form of a 4D matrix
                for i = 1:size(inmasks,4)
                    masks(i).im = inmasks(:,:,:,i);
                end
            end        
    end
end

numMasks = numel(masks);

% Get and Validate Optional Parameters
cMap(1,:) = [1 0 0];
cMap(2,:) = [0 1 0];
cMap(3,:) = [0 0 1];
cMap(4,:) = [1 0 1];
cMap(5,:) = [0 1 1];
cMap(6,:) = [1 1 0];

defaultDisplayRange = ComputeImageDynamicRange( im, 98.0 );

p.addParamValue( 'spacing', ones(1, ndims(im)), @(x) ( isnumeric(x) && ~isscalar(x) && numel(x) == ndims(im) ) );
p.addParamValue( 'displayRange', defaultDisplayRange, @(x) ( isnumeric(x) && numel(x) == 2 ) );
p.addParamValue( 'maskColors', cMap(1:numMasks,:), @(x) (isnumeric(x) && ndims(x) == 2 && size(x,2) == 3 && size(x,1) == numMasks) ); 
p.addParamValue( 'maskAlphas', 0.5 * ones(numMasks,1), @(x) (isnumeric(x) && (isscalar(x) || numel(x) == numMasks)) );
p.parse( im, inmasks, varargin{:} );

spacing = ones(1,3);
spacing(1:ndims(im)) = p.Results.spacing;
maskColorMap = p.Results.maskColors;
displayrange = p.Results.displayRange;

if isscalar(p.Results.maskAlphas) && numMasks > 1
   maskAlpha =  p.Results.maskAlphas  * ones(numMasks,1);
else
   maskAlpha =  p.Results.maskAlphas;
end
    
numColors = size(maskColorMap, 1);
numAlpha  = length(maskAlpha);

if length(unique([ numMasks numColors numAlpha ])) ~= 1
    error('error: number of mask images, alpha values and colors are not the same');
end

%% Data
data.im = im;
data.volSize = volSize;
data.masks = masks;
data.displayrange = displayrange;
data.cmap = maskColorMap;
data.alpha = maskAlpha;

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
data.imLog = ComputeImageLogTransformForDisplay( data.im );
data.imLogDisplayRange = [ min(data.imLog(:)), max(data.imLog(:)) ];

%% Create UI controls
hMainFigure = figure;

% Display figure toolbar
set(hMainFigure, 'ToolBar', 'figure');

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
data.ui.pbh_dec = uicontrol(hMainFigure,'Style','pushbutton','String','<<',...
    'Units' , 'normalized' , 'Position',[0.3 0.11 0.066 0.05],...
    'Callback',{@pushFirstSlice_Callback});

data.ui.pbh_dec = uicontrol(hMainFigure,'Style','pushbutton','String','<',...
    'Units' , 'normalized' , 'Position',[0.366 0.11 0.066 0.05],...
    'Callback',{@pushdec_Callback});

data.ui.eth_sno = uicontrol(hMainFigure,'Style','edit',...
    'String','0',...
    'Units' , 'normalized' , 'Position',[0.433 0.11 0.133 0.05]);

data.ui.pbh_inc = uicontrol(hMainFigure,'Style','pushbutton','String','>',...
    'Units' , 'normalized' , 'Position',[0.566 0.11 0.066 0.05],...
    'Callback',{@pushinc_Callback});

data.ui.pbh_inc = uicontrol(hMainFigure,'Style','pushbutton','String','>>',...
    'Units' , 'normalized' , 'Position',[0.632 0.11 0.066 0.05],...
    'Callback',{@pushLastSlice_Callback});

% Cursor point info controls
data.ui.eth_xloc = uicontrol(hMainFigure,'Style','edit',...
    'String','X: INV',...
    'Units' , 'normalized' , 'Position',[0.3 0.06 0.133 0.05]);

data.ui.eth_yloc = uicontrol(hMainFigure,'Style','edit',...
    'String','Y: INV',...
    'Units' , 'normalized' , 'Position',[0.433 0.06 0.133 0.05]);

data.ui.eth_Imval = uicontrol(hMainFigure,'Style','edit',...
    'String','I: INV',...
    'Units' , 'normalized' , 'Position',[0.566 0.06 0.133 0.05]);

% Masks selection ui
numLimit = 10;
if numel(masks) >= numLimit
    sizeBox = 0.1;
else
    sizeBox = 0.06;
end
data.ui.ch_masks = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.3 0.16 0.4 sizeBox]);
xmask = 5;
ymask = 5;
for i = 1:size(masks, 2)
    if mod(i,numLimit) == 0
        xmask =  xmask + 20;
        ymask = 5;
    end

    % Bool show
    data.masks(i).show = 1;
    data.ui.masks(i) = uicontrol('Style','checkbox',...
        'Units', 'pixels', 'Position',[ymask xmask 15 15],...
        'parent',data.ui.ch_masks,'HandleVisibility','off',...
        'Callback',{@checkMasks_Callback}, 'Value', 1, ...
        'BackgroundColor', maskColorMap(i,:));
    ymask = 25 + ymask;
end

% Planes ui
data.ui.ch_planes = uibuttongroup('visible', 'on', 'Units' , 'normalized', ...
    'Position', [0.30 0.01 0.4 0.05], ...
    'SelectionChangeFcn',{@planes_Callback});

if ndims(data.im) == 3
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
%close(gcbf);

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
function [ slice ] = getRGBSlice(data)

if data.logUse      
    im = data.imLog;
    displayrange = data.imLogDisplayRange;    
else
    im = data.im;
    displayrange = data.displayrange;
end

% Get the slice, depending on the plane showed
if data.curPlane == 1
    imr = squeeze(im(:,data.sliceno,:))';
elseif data.curPlane == 2
    imr = squeeze(im(data.sliceno,:,:))';
else
    imr = squeeze(im(:,:,data.sliceno));
end

imr = im2uint8( mat2gray(imr, displayrange) );

% To rgb
img = imr;
imb = imr;

% Add Masks to rgb
maxVal = 255;
for i = 1:numel(data.masks)
    if data.masks(i).show
        if data.curPlane == 1
            mask = squeeze(data.masks(i).im(:,data.sliceno,:))';
        elseif data.curPlane == 2
            mask = squeeze(data.masks(i).im(data.sliceno,:,:))';
        else
            mask = squeeze(data.masks(i).im(:,:,data.sliceno));
        end
        mask = logical(mask);

        % RGB
        imr(mask) = uint8( double( (1 - data.alpha(i)) * imr(mask) ) + maxVal * data.alpha(i) * data.cmap(i,1) );
        img(mask) = uint8( double( (1 - data.alpha(i)) * img(mask) ) + maxVal * data.alpha(i) * data.cmap(i,2) );
        imb(mask) = uint8( double( (1 - data.alpha(i)) * imb(mask) ) + maxVal * data.alpha(i) * data.cmap(i,3) );
    end
end

% RGB Image
slice = cat(3, imr, img, imb);

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

%% Zoom callback
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

%% Check/Uncheck Masks
function checkMasks_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Update value for each mask
for i = 1:numel(data.masks)
    data.masks(i).show = get(data.ui.masks(i), 'Value');
end

guidata(hSrc, data);
imsliceshow(data);

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

%% Plane to display
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
    data.sliceno = data.sliceno-1;
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
function FnSliceScroll_Callback(hSrc, eventdata) %#ok<INUSD>

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
if IsPointInsideImage(cp(1,1:2), data)
    set(data.ui.eth_xloc, 'String', sprintf('X: %d / %d', cp(1,1), volSize(2)));
    set(data.ui.eth_yloc, 'String', sprintf('Y: %d / %d', cp(1,2), volSize(1)));
    
    ind = getRealCoord(data, [cp(2) cp(1) data.sliceno]);
    set(data.ui.eth_Imval,'String',sprintf('I: %.3f' , data.im(ind)));
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
        set( hSrc ,'Pointer','crosshair');
    else
        set( hSrc ,'Pointer','arrow');
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

    alpha = 0.5 * (1 - 0.01 * cover_percent);
    intensityRange = quantile(double(im(:)), [alpha, 1 - alpha]);

%     [p,x] = hist( im, 255 );   
%     p = p / sum(p);
%     
%     min_xlow = [];
%     min_xhigh = [];
%     min_xdiff = [];
%     
%     for i = 1:numel(x)
%         for j = i+1:numel(x)
%     
%             if sum( p(i:j) ) < 0.01 * cover_percent
%                 continue;
%             end
%             
%             if isempty(min_xdiff) || (x(j) - x(i)) < min_xdiff
%                 min_xlow = x(i);
%                 min_xhigh = x(j);
%                 min_xdiff = x(j) - x(i);              
%             end
%         end
%     end
%     
%     w = 0.5 * (x(2) - x(1));
%     intensityRange = [min_xlow-w, min_xhigh+w];
    
end

%%
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 98.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
end