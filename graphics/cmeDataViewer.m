%cmeDataViewer(data, varargin) displays movies with associated detection and tracking results.
%
% Inputs:    
%             data : single movie structure returned by loadConditionData.m
%     Trajectories : optional input for selecting 'all' (default) or
%                    'valid' CCS trajectories.
%c
% Notes: Only tracks with at least 5 frames are loaded and displayed.

% Francois Aguet, 2011 (last modified 08/24/2013)

function cmeDataViewer(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('Trajectories', 'all', @(x) isempty(x) || isstruct(x) || any(strcmpi(x, {'all', 'valid'})));
ip.addParamValue('LoadTracks', true, @islogical);
ip.parse(data, varargin{:});

% Handles/settings are stored in 'appdata' of the figure handle
handles.data = data;

% detect number of channels (up to 4)
nCh = length(data.channels);
if nCh>4
    error('Max. 4 channels supported.');
end

handles.nCh = nCh;
% master channel index
handles.mCh = find(strcmp(data.source, data.channels));


fidx = 1;
tcur = 1;
xs = [];
ys = [];

nx = data.imagesize(2);
ny = data.imagesize(1);
nf = data.movieLength;

lcolor = hsv2rgb([0.55 0.5 0.8]);

%===============================================================================
% Setup main GUI window/figure
%===============================================================================
hfig = figure('Units', 'normalized', 'Position', [0.025 0.2 0.95 0.8],...
    'PaperPositionMode', 'auto', 'Toolbar', 'figure',...
    'Color', get(0,'defaultUicontrolBackgroundColor'),...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', [getDirFromPath(getExpDir(data)) filesep getCellDir(data)]);

pos = get(hfig, 'Position'); % [pixels]

%-------------------------------------------------------------------------------
% Control panels at bottom of GUI
%-------------------------------------------------------------------------------
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', '', 'Position', [5 5 650 70]);

uicontrol(ph, 'Style', 'text', 'String', 'Data display: ',...
    'Position', [5 40 90 20], 'HorizontalAlignment', 'left');

maskCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Cell mask',...
    'Position', [5 25 100 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);
% plot on top
frameChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Raw', 'Detections', 'RGB'},...
    'Position', [95 42 100 20], 'Callback', @frameChoice_Callback);

labelCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Channel labels',...
    'Position', [5 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @chlabel_Callback);

detectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Detections',...
    'Position', [200 50 100 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);
trackCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Tracks:', 'Value', true,...
    'Position', [200 30 80 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);
trackChoice = uicontrol('Style', 'popup',...
    'String', {'Category', 'Lifetime', 'EAP Status', 'Object Type', 'Random'},...
    'Position', [280 33 100 20], 'Callback', @trackChoice_Callback);
trackRangeButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Settings',...
    'Position', [280 5 80 20], 'HorizontalAlignment', 'left',...
    'Callback', @trackSettings_Callback);


gapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Gaps',...
    'Position', [390 45 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);
trackEventCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Births/Deaths',...
    'Position', [390 25 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);
eapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'EAP status',...
    'Position', [390 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @updateSlice);



trackButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Select track',...
    'Position', [540 40 100 20], 'HorizontalAlignment', 'left',...
    'Callback', @trackButton_Callback);
statsButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Track statistics',...
    'Position', [540 10 100 20], 'HorizontalAlignment', 'left',...
    'Callback', @statsButton_Callback);


%---------------------
% Tracks
%---------------------
% Track plot panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Plot options', 'Position', [pos(3)-220-150-5-160 5 160 70]);
tplotText = uicontrol(ph, 'Style', 'text', 'String', 'Units: ',...
    'Position', [5 35 60 20], 'HorizontalAlignment', 'left');

tplotUnitChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Seconds', 'Frames'},...
    'Position', [40 40 100 15], 'Callback', {@unitChoice_Callback, hfig});
tplotBackgroundCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Subtract background',...
    'Position', [5 20 150 15], 'HorizontalAlignment', 'left', 'Value', true, 'Callback', @updateTrack);
tplotScaleCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Autoscale',...
    'Position', [5 5 120 15], 'HorizontalAlignment', 'left', 'Value', false, 'Callback', @updateTrack);
handles.tplotPanel = ph;

% Montage panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Montage plot', 'Position', [pos(3)-220-150 5 220 70]);
montageAlignCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Align to track',...
    'Position', [95 38 115 15], 'HorizontalAlignment', 'left', 'Value', true);
montageMarkerCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show markers',...
    'Position', [95 23 115 15], 'HorizontalAlignment', 'left');
montageDetectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show detection',...
    'Position', [95 8 115 15], 'HorizontalAlignment', 'left');
montageButton = uicontrol(ph, 'Style', 'pushbutton','String','Generate',...
    'Units', 'pixels', 'Position', [10 20 80 20],...
    'Callback', @montageButton_Callback);
handles.montagePanel = ph;


% Output panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Output', 'Position', [pos(3)-145 5 140 70]);

handles.printButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Print figures',...
    'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.45],...
    'Callback', {@printButton_Callback, hfig});

handles.movieButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Make movie',...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.45],...
    'Callback', @movieButton_Callback);
handles.outputPanel = ph;


setappdata(hfig, 'handles', handles); % write 'handles' to hfig

%===============================================================================
% Set up frame display
%===============================================================================
% fixed width of the track plots, in pixels
w = 320;

tspace = 20;
bspace = 100;
lspace = 10;
rspace = w+30+50;
spacer = 10; 
handles = setupFrameAxes(hfig, [lspace bspace rspace tspace spacer]);

handles.frameLabel = uicontrol('Style', 'text', 'String', ['Frame ' num2str(fidx)], ...
    'Position', [10 pos(4)-20 100 15], 'HorizontalAlignment', 'left');

% Frame slider
if data.movieLength>1
    handles.frameSlider = uicontrol('Style', 'slider', 'Units', 'pixels',...
        'Value', fidx, 'SliderStep', [1/(nf-1) 0.05], 'Min', 1, 'Max', nf,...
        'Position', [lspace 77 pos(3)-rspace-lspace 18]);
end
% this definition (instead of regular callback) enable continuous sliding
addlistener(handle(handles.frameSlider), 'Value', 'PostSet', @frameSlider_Callback);


%===============================================================================
% Set up track display
%===============================================================================
% track panels: 20 spacer, 110 bottom, 20 top
spacer = 15;
% if ~isempty(tracks)
    h_tot = pos(4) - 140;
    h = min((h_tot-(nCh-1)*spacer)/nCh, 200);
    
    opts = {'Parent', hfig, 'Units', 'pixels', 'Box', 'on'};
    dx = pos(3)-w-30;
    switch nCh
        case 1
            handles.tAxes(1) = axes(opts{:}, 'Position', [dx 120+(h_tot-h) w h]);
        case 2
            handles.tAxes(1) = axes(opts{:}, 'Position', [dx 120+(h_tot-h) w h]);
            handles.tAxes(2) = axes(opts{:}, 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
        case 3
            handles.tAxes(1) = axes(opts{:}, 'Position', [dx 120+(h_tot-h) w h]);
            handles.tAxes(2) = axes(opts{:}, 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
            handles.tAxes(3) = axes(opts{:}, 'Position', [dx 120+(h_tot-3*h-2*spacer) w h]);
        case 4
            handles.tAxes(1) = axes(opts{:}, 'Position', [dx 120+(h_tot-h) w h]);
            handles.tAxes(2) = axes(opts{:}, 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
            handles.tAxes(3) = axes(opts{:}, 'Position', [dx 120+(h_tot-3*h-2*spacer) w h]);
            handles.tAxes(4) = axes(opts{:}, 'Position', [dx 120+(h_tot-4*h-3*spacer) w h]);
    end

    handles.trackLabel = uicontrol('Style', 'text', 'String', 'Track 1',...
        'Units', 'pixels', 'Position', [pos(3)-70 pos(4)-20 100 15], 'HorizontalAlignment', 'left');
    
    handles.trackSlider = uicontrol('Style', 'slider',...
        'Value', 1, 'SliderStep', [0.01 0.05], 'Min', 1, 'Max', 1000,...
        'Position', [pos(3)-24 120 10 h_tot]);
    % this definition (instead of regular callback) enable continuous sliding
    addlistener(handle(handles.trackSlider), 'Value', 'PostSet', @trackSlider_Callback);
% end


handles.fAxes = zeros(nCh,3);
hLegend = zeros(1,nCh);
for c = 1:nCh   
    [handles.fAxes(c,:), hLegend(c)] = setupStackViewer(handles.fPanels(c), [nx ny min(nf,  max(nx,ny)/3)], c==1); 
end
hLegend = hLegend(1);
colormap(gray(256));

setappdata(hfig, 'handles', handles);
set(hfig, 'ResizeFcn', @figResize);

% handles for track plot objects in frames window
hpt = []; % tracks
hpd = []; % detections
hpg = []; % gaps
hps = []; % starts/ends

hst = []; % selected track marker

hms = []; % cell mask

% handles for track plots
% ht = [];

cmap = [];

displayType = 'raw';
pUnitType = 's';


%===============================================================================
% Load movie and associated analysis results
%===============================================================================
% readfct = @(path, i) imread(path, i);
fprintf('Loading frames ... ');
stack = cell(1,nCh);
if ~iscell(data.framePaths{1})
    for c = 1:nCh
        %stack{c} = readtiff(data.framePaths{c});
        stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
        for i = 1:data.movieLength
            stack{c}(:,:,i) = imread(data.framePaths{c}, i);
        end
    end
else
    for c = 1:nCh
        stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
        for i = 1:data.movieLength
            stack{c}(:,:,i) = imread(data.framePaths{c}{i});
        end
    end
end
fprintf('done.\n');

%-------------------------------------------------------------------------------
% Load detection masks
%-------------------------------------------------------------------------------
fprintf('Loading detection masks ... ');
dpath = [data.source 'Detection' filesep 'detection_v2.mat'];
if exist(dpath, 'file')==2
    dmask = zeros(ny,nx,nf, 'uint8');
    if ~iscell(data.framePaths{1})
        for i = 1:nf
            dmask(:,:,i) = imread(data.maskPaths, i);
        end
    else
        for i = 1:data.movieLength
            dmask(:,:,i) = imread(data.maskPaths{i});
        end
    end
else
    dmask = [];
end

if exist([data.source 'Detection' filesep 'cellmask.tif'], 'file')==2
    cellMask = imread([data.source 'Detection' filesep 'cellmask.tif']);
else
    cellMask = [];
end


fprintf('done.\n');
%-------------------------------------------------------------------------------
% Load detection files
%-------------------------------------------------------------------------------
% for c = 1:nCh
detectionFile = [data.channels{1} 'Detection' filesep 'detection_v2.mat'];
if (exist(detectionFile, 'file')==2)
    frameInfo = load(detectionFile);
    frameInfo = frameInfo.frameInfo;
else
    frameInfo = [];
end
% end

%-------------------------------------------------------------------------------
% Load tracks
%-------------------------------------------------------------------------------
fprintf('Loading tracks ... ');
tracks = [];
bgA = [];
maxA = [];
% identify track file
fileList = dir([data.source 'Tracking' filesep 'ProcessedTracks*.mat']);
fileList = {fileList.name};
if numel(fileList)>1
    idx = 0;
    while ~(idx>=1 && idx<=numel(fileList) && round(idx)==idx)
        fprintf('Tracking results found for this data set:\n');
        for i = 1:numel(fileList)
            fprintf('[%d] %s\n', i, fileList{i});
        end
        idx = str2double(input('Please enter the number of the set to load: ', 's'));
    end
    fileName = fileList{idx};
elseif numel(fileList)==1
    fileName = fileList{1};
else
    fileName = [];
end   

if exist([data.source 'Tracking' filesep fileName], 'file')==2 && ip.Results.LoadTracks
    tmp = load([data.source 'Tracking' filesep fileName]);
    tracks = tmp.tracks;
    if isfield(tmp, 'bgA')
        bgA = cellfun(@(i) prctile(i, 95, 2), tmp.bgA, 'unif', 0);
        bgA = [bgA{:}];
    end
    clear tmp;
    cutoff_f = 3; % by default; this is adjustable in the GUI
    tracks = tracks([tracks.lifetime_s] >= data.framerate*cutoff_f);
    [~, sortIdx] = sort([tracks.lifetime_s], 'descend');
    tracks = tracks(sortIdx);
    
    % apply cell mask
    if ~isempty(cellMask)
        nt = numel(tracks);
        x = NaN(1,nt);
        y = NaN(1,nt);
        for t = 1:nt
            x(t) = round(nanmean(tracks(t).x(1,:)));
            y(t) = round(nanmean(tracks(t).y(1,:)));
        end
        idx = sub2ind([ny nx], y, x);
        tracks = tracks(cellMask(idx)==1);
    end
    nt = numel(tracks);
    selIndex = true(1,nt);

    nseg = [tracks.nSeg];
    
    np = sum(nseg);
    X = NaN(nf, np);
    Y = NaN(nf, np);
    G = false(nf, np);
    % for significance values, store vectors
    mvec = [tracks.hval_Ar];
    if isfield(tracks, 'significantVsBackground')
        svec = [tracks.significantVsBackground];
    else
        svec = [];
    end
    fvec = [tracks.f];
    xvec = [tracks.x];
    yvec = [tracks.y];
   
    % vector of start indexes since multiple segments/track
    tidx = cumsum([1 nseg(1:end-1)]);
    
    trackStarts = [tracks.start];
    trackEnds = [tracks.end];
    mu_x = NaN(1,nt);
    mu_y = NaN(1,nt);
    
    for t = 1:nt
        if nseg(t)==1
            X(tracks(t).f, tidx(t)) = tracks(t).x(1,:);
            Y(tracks(t).f, tidx(t)) = tracks(t).y(1,:);
            G(tracks(t).f, tidx(t)) = tracks(t).gapVect;
            mu_x(t) = nanmean(X(:,tidx(t)));
            mu_y(t) = nanmean(Y(:,tidx(t)));           
        else
            sep = find(isnan(tracks(t).t));
            sep = [0 sep numel(tracks(t).f)+1]; %#ok<AGROW>
            for s = 1:tracks(t).nSeg
                sidx = sep(s)+1:sep(s+1)-1;
                X(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).x(1,sidx);
                Y(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).y(1,sidx);
                G(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).gapVect(sidx);
            end
            mu_x(t) = nanmean(nanmean(X(:,tidx(t):tidx(t)+s-1)));
            mu_y(t) = nanmean(nanmean(Y(:,tidx(t):tidx(t)+s-1)));
        end
    end
    
    % segment index
    % Example: [1 1 2 3 4 4 ... ] first two cols are from same track
    idx = diff([tidx size(X,2)+1]);
    idx = arrayfun(@(i) i+zeros(1, idx(i)), 1:numel(idx), 'unif', 0);
    tstruct.idx = [idx{:}];
    tstruct.n = numel(tracks);
    
    % min/max track intensities
    maxA = arrayfun(@(t) max(t.A, [], 2), tracks, 'unif', 0);
    maxA = [maxA{:}];
    maxVal = prctile(maxA, 99, 2);
    da = floor(log10(maxVal));
    % y-axis unit
    yunit = round(maxVal ./ 10.^da) .* 10.^(da-1);
    maxVal = ceil(maxVal./yunit) .* yunit;
end
fprintf('done.\n');

% Settings
if ~isempty(tracks)
    minLft = min([tracks.lifetime_s]); % 3 frames. Display tracks >= 5 framesby default
    maxLft = max([tracks.lifetime_s]);
    minVal = data.framerate*5;%minLft;
    maxVal = maxLft;
    catCheckVal = ones(1,8);
    eapCheckVal = ones(1,3);
    maxIntT = 0;
end





% dynamic range for each channel
dRange = cell(1,nCh);
for c = 1:nCh
    dRange{c} = double([min(stack{c}(:)) max(stack{c}(:))]);
end
% dRange{1} = prctile(double(stack{c}(:)), [1 99]);

hues = getFluorophoreHues(data.markers);
rgbColors = arrayfun(@(x) hsv2rgb([x 1 1]), hues, 'unif', 0);


%===============================================================================
% Set visibility for sliders and checkboxes
%===============================================================================
if ~isempty(tracks)
    set(handles.trackSlider, 'Min', 1);
    set(handles.trackSlider, 'Max', nt);
    set(handles.trackSlider, 'SliderStep', [1/(nt-1) 0.05]);
    setTrackColormap('Category');
    setColorbar('Category');
else
    set(handles.trackSlider, 'Visible', 'off');
    set(handles.trackLabel, 'Visible', 'off');
    set(handles.tAxes, 'Visible', 'off');
    set(trackButton, 'Enable', 'off');
    set(statsButton, 'Enable', 'off');
    set(trackCheckbox, 'Value', false);
    set([trackCheckbox trackChoice trackRangeButton gapCheckbox trackEventCheckbox], 'Enable', 'off');
    set([montageAlignCheckbox montageMarkerCheckbox montageDetectionCheckbox montageButton], 'Enable', 'off');
    %set(handles.montagePanel, 'Visible', 'off');
    
    set(hLegend, 'Visible', 'off');
    set([tplotText tplotUnitChoice tplotBackgroundCheckbox tplotScaleCheckbox], 'Enable', 'off');
end
if isempty(cellMask)
    set(maskCheckbox, 'Enable', 'off');
end

if nCh==1
    set(eapCheckbox, 'Enable', 'off');
end

%===============================================================================
% populate with data, plotting functions are called only here, afterwards change data
%===============================================================================
x = round(nx/2);
y = round(ny/2);
hxy = zeros(1,nCh);
hyz = zeros(1,nCh);
hxz = zeros(1,nCh);
hl = zeros(nCh,4);
for c = 1:nCh
    % x,y view
    hxy(c) = imagesc(stack{c}(:,:,fidx), 'Parent', handles.fAxes(c,1), 'HitTest', 'off');
    hold(handles.fAxes(c,1), 'on');
    set(handles.fAxes(c,1), 'ButtonDownFcn', @click_Callback);
    hl(c,1) = plot(handles.fAxes(c,1), [x x], [0.5 ny+0.5], 'Color', lcolor, 'HitTest', 'off', 'DisplayName', 'FrameMarker');
    hl(c,2) = plot(handles.fAxes(c,1), [0.5 nx+0.5], [y y], 'Color', lcolor, 'HitTest', 'off', 'DisplayName', 'FrameMarker');
    
    % y,z view
    hyz(c) = imagesc(squeeze(stack{c}(:,x,:)), 'Parent', handles.fAxes(c,2), 'HitTest', 'off');
    hold(handles.fAxes(c,2), 'on');
    % line in y,z view
    hl(c,3) = plot(handles.fAxes(c,2), fidx*[1 1], [0.5 ny+0.5], 'Color', lcolor, 'HitTest', 'off');
    hold(handles.fAxes(c,2), 'off');
    
    % x,z view
    hxz(c) = imagesc(squeeze(stack{c}(y,:,:))', 'Parent', handles.fAxes(c,3), 'HitTest', 'off');
    hold(handles.fAxes(c,3), 'on');
    % line in x,z view
    hl(c,4) = plot(handles.fAxes(c,3), [0.5 nx+0.5], fidx*[1 1], 'Color', lcolor, 'HitTest', 'off');
    hold(handles.fAxes(c,3), 'off');
    
    arrayfun(@(i) caxis(i, dRange{c}), handles.fAxes(c,:), 'unif', 0);
end
set(handles.fAxes, 'XTick', [], 'YTick', []);
axis(handles.fAxes(:,1), 'equal');
% this fixes a bug with axis 'equal' that allows panning beyond boundaries
set(handles.fAxes(:,1), 'XLim', [0.5 nx+0.5], 'YLim', [0.5 ny+0.5]);

set(handles.fAxes, 'ButtonDownFcn', @click_Callback);

dx = 0.03;
hChLabel = zeros(1,nCh);
for c = 1:nCh
    hChLabel(c) = text(1-dx*ny/nx, dx, data.markers{c},...
        'Color', rgbColors{c}, 'Units', 'normalized',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
        'Parent', handles.fAxes(c,1), 'HitTest', 'off');
end
if nCh <= 2
    set(hChLabel, 'Visible', 'off');
end

if ~isempty(tracks)
    
    
    % plot current frame marker
    hf = zeros(nCh,1);
    hst = zeros(1,nCh);
    for c = 1:nCh
        hf(c) = plot(handles.tAxes(c), ([fidx fidx]-1)*data.framerate,...
            get(handles.tAxes(c), 'YLim'), '--', 'Color', 0.7*[1 1 1]);
        hst(c) = plot(handles.fAxes(c,1), X(fidx, tstruct.idx==tcur),...
                    Y(fidx, tstruct.idx==tcur), 'ws', 'DisplayName', 'TrackMarker', 'MarkerSize', 12);%*nx/diff(get(handles.fAxes(c,1),'XLim')));
    end
    
    updateTrack();
end

% if ~isempty(ip.Results.Trajectories)
% 
%     % load tracks
%     if ischar(ip.Results.Trajectories)
%         if strcmpi(ip.Results.Trajectories, 'valid');
%             c = 'Ia';
%         else
%             c = 'all';
%         end
%         tracks = loadTracks(data, 'Category', c);
%     else
%         tracks = ip.Results.Trajectories;
%     end
% end

%===============================================================================
% Set listeners
%===============================================================================
set(hfig, 'WindowScrollWheelFcn', @scroll_Callback);
set(hfig, 'KeyPressFcn', @key_Callback);

hpan = pan;
set(hpan,'ActionPreCallback',@panstart);
set(hpan,'ActionPostCallback',@panstop);

hz = zoom;
set(hz, 'ActionPostCallback', @czoom);
% setAxesZoomMotion(hz, handles.tAxes, 'horizontal');
% linkaxes(handles.tAxes, 'x');


%===============================================================================
% Listener/display functions
%===============================================================================
    function click_Callback(varargin)
        switch gca
            case num2cell(handles.fAxes(:,1))
                a = get(gca, 'CurrentPoint');
                xs = round(a(1,1));
                ys = round(a(1,2));
                updateProj(); % required for clicking w/o dragging
                set(gcf, 'WindowButtonMotionFcn', @dragProj, 'WindowButtonUpFcn', @stopDragging);
            case num2cell(handles.fAxes(:,2))
                a = get(gca,'CurrentPoint');
                fidx = round(a(1,1));
                updateSlice();
                set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);
            case num2cell(handles.fAxes(:,3))
                a = get(gca,'CurrentPoint');
                fidx = round(a(1,2));
                updateSlice();
                set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);                
        end
    end

    function dragProj(varargin)
        a = get(gca, 'CurrentPoint');
        xs = round(a(1,1));
        ys = round(a(1,2));
        updateProj();
    end

    function dragSlice(varargin)
        a = get(gca, 'CurrentPoint');
        switch gca
            case num2cell(handles.fAxes(:,2))
                fidx = min(max(1,round(a(1,1))),nf);
            case num2cell(handles.fAxes(:,3))
                fidx = min(max(1,round(a(1,2))),nf);
        end
        set(handles.frameSlider, 'Value', fidx);
        updateSlice();
    end

    function stopDragging(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end


    % scroll through stack slices
    function scroll_Callback(~, eventdata)
        if eventdata.VerticalScrollCount < 0
            if fidx < nf
                fidx = fidx + 1;
                set(handles.frameSlider, 'Value', fidx);
                updateSlice();
            end
        elseif eventdata.VerticalScrollCount > 0
            if fidx > 1
                fidx = fidx - 1;
                set(handles.frameSlider, 'Value', fidx);
                updateSlice();
            end
        end
    end


    function key_Callback(~, eventdata)
        switch eventdata.Key
            case 'leftarrow'
                if fidx > 1
                    fidx = fidx - 1;
                    set(handles.frameSlider, 'Value', fidx);
                    updateSlice();
                end
            case 'rightarrow'
                if fidx < nf
                    fidx = fidx + 1;
                    set(handles.frameSlider, 'Value', fidx);
                    updateSlice();
                end
            case 'downarrow'
                if tcur > find(selIndex, 1, 'first');
                    % this invokes trackSlider_callback
                    set(handles.trackSlider, 'Value', get(handles.trackSlider, 'Value')-1);
                end
            case 'uparrow'
                if tcur < find(selIndex, 1, 'last');
                    set(handles.trackSlider, 'Value', get(handles.trackSlider, 'Value')+1);
                end                        
        end
    end



    function updateSlice(varargin)       
        switch displayType
            case 'raw'                
                for ci = 1:nCh
                    set(hxy(ci), 'CData', stack{ci}(:,:,fidx));
                end
            case 'mask'
                set(hxy(1), 'CData', rgbOverlay(stack{1}(:,:,fidx), dmask(:,:,fidx), [1 0 0], dRange{1}));
                for ci = 2:nCh
                    set(hxy(ci), 'CData', stack{ci}(:,:,fidx));
                end
            case 'RGB'
                rframe = zeros(ny,nx,3,'uint8');
                idxRGB = getRGBindex(data.markers);                
                for ci = 1:nCh
                    rframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(stack{ci}(:,:,fidx)), dRange{ci}));
                end
                set(hxy(1), 'CData', rframe);
        end
        
        set(hl(:,3), 'XData', fidx*[1 1]);
        set(hl(:,4), 'YData', fidx*[1 1]);        
        set(handles.frameLabel, 'String', ['Frame ' num2str(fidx)]);
        
        delete(hms);
        hms = [];
        if ~isempty(cellMask) && get(maskCheckbox, 'Value')
            B = bwboundaries(cellMask);
            for ci = 1:nCh
                hms = [hms; cellfun(@(i) plot(handles.fAxes(ci,1), i(:,2), i(:,1), 'Color', 'r', 'LineWidth', 1), B)]; %#ok<AGROW>
            end
        end
        
        if ~isempty(tracks)
            % update current frame marker in track plots
            for ci = 1:nCh
                set(hf(ci), 'XData', ([fidx fidx]-1)*data.framerate,...
                    'YData', get(handles.tAxes(ci), 'YLim'));
            end
        end
        
        delete(hpt);
        delete(hpg);
        delete(hps);
        hpt = [];
        hpg = [];
        hps = [];

        if ~isempty(tracks) && fidx~=1 && get(trackCheckbox, 'Value') && any(~isnan(X(fidx,:)) & selIndex(tstruct.idx))
            vidx = ~isnan(X(fidx,:)) & selIndex(tstruct.idx);
            delete(hpt);
            set(handles.fAxes(1,1), 'ColorOrder', cmap(tstruct.idx(vidx),:));
            hpt = plot(handles.fAxes(1,1), X(1:fidx,vidx), Y(1:fidx,vidx), 'HitTest', 'off');
            if get(gapCheckbox, 'Value')
                hpg = plot(handles.fAxes(1,1), X(fidx,vidx & G(fidx,:)), Y(fidx,vidx & G(fidx,:)), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
            end
            if get(trackEventCheckbox, 'Value')
                % Births
                bcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackStarts==fidx & selIndex), 'unif', 0);
                bcoord = vertcat(bcoord{:});
                if~isempty(bcoord)
                    hps = plot(handles.fAxes(1,1), bcoord(:,1), bcoord(:,2), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                end
                
                % Deaths
                dcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackEnds==fidx & selIndex), 'unif', 0);
                dcoord = vertcat(dcoord{:});
                if ~isempty(dcoord)
                    hps = [hps; plot(handles.fAxes(1,1), dcoord(:,1), dcoord(:,2), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1)];
                end
            end
            set(hst, 'Visible', 'on', 'XData', X(fidx, tstruct.idx==tcur), 'YData', Y(fidx, tstruct.idx==tcur));
        else
            set(hst, 'Visible', 'off');
        end
        if ~isempty(tracks) && get(eapCheckbox, 'Value')
            for ci = 2:nCh
                sel = fvec==fidx & mvec(ci,:)==1;
                hp1 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([1/3 1 0.9]), 'MarkerSize', 8);
                sel = fvec==fidx & mvec(ci,:)==0 & svec(ci,:)==1;
                hp2 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([0.55 1 0.9]), 'MarkerSize', 8);
                sel = fvec==fidx & mvec(ci,:)==0 & svec(ci,:)==0;
                hp3 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', 0.8*[1 1 1], 'MarkerSize', 8);
                hps = [hps; hp1; hp2; hp3]; %#ok<AGROW>
            end
        end
        
        delete(hpd); % clear previous plots
        hpd = [];
        if get(detectionCheckbox, 'Value') && ~isempty(frameInfo)
            isPSF = frameInfo(fidx).isPSF(1,:)==1;
            if any(isPSF)
                hpd(1) = plot(handles.fAxes(1,1), frameInfo(fidx).x(1,isPSF), frameInfo(fidx).y(1,isPSF), 'o', 'Color', [0 0.6 0], 'MarkerSize', 8);
            end
            if any(~isPSF)
                hpd(2) = plot(handles.fAxes(1,1), frameInfo(fidx).x(1,~isPSF), frameInfo(fidx).y(1,~isPSF), 'o', 'Color', [0.6 0 0], 'MarkerSize', 8);
            end
        end
    end


    function updateProj()
        % plot lines
        set(hl(:,1), 'XData', xs*[1 1]);
        set(hl(:,2), 'YData', ys*[1 1]);
        
        % update data
        xi = min(max(xs,1), nx);
        yi = min(max(ys,1), ny);
        
%         switch displayType
%             case 'RGB'
%                 idxRGB = getRGBindex(data.markers);
%                 tframe = zeros(nf,nx,3,'uint8');
%                 lframe = zeros(ny,nf,3,'uint8');
%                 for ci = 1:nCh
%                     tframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(yi,:,:))'), dRange{ci}));
%                     lframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(:,xi,:))), dRange{ci}));
%                 end
%                 set(hxz(1), 'CData', tframe);
%                 set(hyz(1), 'CData', lframe);
%             %case 'mask'
%             %    set(hxy(1), 'CData', rgbOverlay(stack{1}(:,:,fidx), dmask(:,:,fidx), [1 0 0], dRange{1}));
%             %    for ci = 2:nCh
%             %        set(hxy(ci), 'CData', stack{ci}(:,:,fidx));
%             %    end
%             otherwise
                for ci = 1:nCh
                    set(hyz(ci), 'CData', squeeze(stack{ci}(:,xi,:)));
                    set(hxz(ci), 'CData', squeeze(stack{ci}(yi,:,:))');
                end
%         end
    end


    function frameChoice_Callback(varargin)
        contents = cellstr(get(frameChoice,'String'));
        switch contents{get(frameChoice,'Value')}
            case 'Raw'
                displayType = 'raw';
            case 'RGB'
                displayType = 'RGB';
            case 'Detections'
                displayType = 'mask';
        end
        updateSlice();
        %updateProj();
    end


    function czoom(~, eventdata)
        % identify panel
        ci = handles.fAxes(:,1) == eventdata.Axes;
        if any(ci) %&& nCh>1 % x,y axes zoomed
            XLim = get(handles.fAxes(ci,1), 'XLim');
            YLim = get(handles.fAxes(ci,1), 'YLim');
            set(handles.fAxes(:,1), 'XLim', XLim, 'YLim', YLim);
            set(handles.fAxes(:,2), 'YLim', YLim);
            set(handles.fAxes(:,3), 'XLim', XLim);
        end
        
    end

    % Pan functions
    function panstart(~, eventdata)
        set(hfig, 'WindowButtonMotionFcn', {@dopan, eventdata});
    end

    function panstop(varargin)
        set(hfig, 'WindowButtonMotionFcn', '');
    end

    function dopan(~,~,eventdata)
        % get limits of current axes
        XLim = get(eventdata.Axes, 'XLim');
        YLim = get(eventdata.Axes, 'YLim');
        
        switch find(any(handles.fAxes == eventdata.Axes,1))
            case 1
                set(handles.fAxes(:,1), 'XLim', XLim, 'YLim', YLim);
                set(handles.fAxes(:,2), 'YLim', YLim);
                set(handles.fAxes(:,3), 'XLim', XLim);
            case 2
                set(handles.fAxes(:,[1 2]), 'YLim', YLim);
            case 3
                set(handles.fAxes(:,[1 3]), 'XLim', XLim);
        end
    end


    function frameSlider_Callback(~, eventdata)
        obj = get(eventdata, 'AffectedObject'); % this contains the current, continuous value
        fidx = round(get(obj, 'Value'));
        updateSlice();
    end


    function trackSlider_Callback(~, eventdata)
        obj = get(eventdata, 'AffectedObject');
        t0 = round(get(obj, 'Value'));
        tmp = find(selIndex);
        tcur = tmp(t0);
        updateTrack();
        
        % if track not visible, jump to first frame
        % t = handles.tracks{1}(t);
        % if fidx < t.start || fidx > t.end
        %     fidx = t.start;
        %     % set frame number
        %     set(handles.frameLabel, 'String', ['Frame ' num2str(fidx)]);
        %     % set frame slider
        %     set(handles.frameSlider, 'Value', fidx);
        % end
    end


    function updateTrack(varargin)
        
        set(handles.trackLabel, 'String', ['Track ' num2str(tcur)]);
        % update selected track marker position
        set(hst, 'XData', X(fidx, tstruct.idx==tcur), 'YData', Y(fidx, tstruct.idx==tcur));

        
        itrack = tracks(tcur);
        for ci = 1:nCh
            %cla(handles.tAxes(ci));
            hold(handles.tAxes(ci), 'off');
            if get(tplotBackgroundCheckbox, 'Value')
                bgMode = 'zero';
            else
                bgMode = 'data';
            end
            if strcmpi(pUnitType, 'f')
                itrack.t = itrack.f;
                if ~isempty(itrack.startBuffer)
                    itrack.startBuffer.t = itrack.f(1) - (numel(itrack.startBuffer.t):-1:1);
                    itrack.endBuffer.t = itrack.f(end) + (1:numel(itrack.startBuffer.t));
                end
            end
            topts = {'Handle', handles.tAxes(ci), 'Time', 'Movie', 'BackgroundValue', bgMode};
            if get(tplotScaleCheckbox, 'Value')
                topts = [topts, 'YTick', -yunit(ci):yunit(ci):maxVal(ci)]; %#ok<AGROW>
            end
            
            if ~isempty(bgA) && itrack.catIdx<5
                conf = bgA(ci, itrack.f); 
               if ~isempty(itrack.startBuffer)
                   conf = [bgA(ci, itrack.startBuffer.f) conf]; %#ok<AGROW>
               end
               if ~isempty(itrack.endBuffer)
                   conf = [conf bgA(ci, itrack.endBuffer.f)]; %#ok<AGROW>
               end
               topts = [topts 'BackgroundConfidence', conf]; %#ok<AGROW>
            end
            plotTrack(data, itrack, ci, topts{:});
            hold(handles.tAxes(ci), 'on');
            %         dx = 0.03;
            %         if isfield(sTrack, 'significantSignal')
            %             s = sTrack.significantSignal;
            %             if s(ci)==1
            %                 slabel = 'yes';
            %                 scolor = [0 0.8 0];
            %             else
            %                 slabel = 'no';
            %                 scolor = [0.8 0 0];
            %             end
            %             text(1-dx, 1-dx,...
            %                 ['Significant: ' slabel],...
            %                 'Color', scolor, 'Units', 'normalized',...
            %                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
            %                 'Parent', handles.tAxes(ci));
            %         end
            
            hf(ci) = plot(handles.tAxes(ci), ([fidx fidx]-1)*data.framerate,...
                get(handles.tAxes(ci), 'YLim'), '--', 'Color', 0.7*[1 1 1]);
            
            set(handles.tAxes(1:end-1), 'XTickLabel', []);
            
            if strcmpi(pUnitType, 's')
                xlabel(handles.tAxes(end), 'Time (s)');
            else
                xlabel(handles.tAxes(end), 'Frames');
            end
        end
    end

    function trackChoice_Callback(~,~)
        str = cellstr(get(trackChoice, 'String'));
        str = str{get(trackChoice,'Value')};
        setTrackColormap(str);
        setColorbar(str);
        updateSlice();
    end


    function unitChoice_Callback(varargin)       
        contents = cellstr(get(tplotUnitChoice,'String'));
        switch contents{get(tplotUnitChoice,'Value')}
            case 'Seconds'
                pUnitType = 's';
            case 'Frames'
                pUnitType = 'f';
        end
        updateTrack();
    end


    function setTrackColormap(mode)
        switch mode
            case 'Category'
                cmap = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
                cmap = cmap([tracks.catIdx],:);
            case 'Lifetime'
                lifetimes_f = round([tracks.lifetime_s]/data.framerate);
                df = data.movieLength-round(120/data.framerate);
                dcoord = 0.25/df;
                cmap = [jet(round(120/data.framerate)); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                cmap = cmap(lifetimes_f,:);
            case 'EAP Status'
                cmap = hsv2rgb([0 0 0.8; 0.55 1 0.9; 0.33 1 0.9]); % ns, slave sig., master sig.
                S = [tracks.significantSlave];
                M = [tracks.significantMaster];
                eap = ones(1,nt);
                eap(M(2,:)==1) = 3;
                eap(S(2,:)==1 & M(2,:)==0) = 2;
                cmap = cmap(eap,:);                
            case 'Object Type'
                isCCP = [tracks.isCCP];
                cmap = [0.8 0 0; 0 0.8 0];
                cmap = cmap(isCCP+1,:);
            case 'Random'
                cmap = hsv2rgb([rand(nt,1) ones(nt,2)]);
        end
    end

        
    function chlabel_Callback(~,~)
        if get(labelCheckbox, 'Value') %&& ~isRGB
            set(hChLabel, 'Visible', 'on');
        else
            set(hChLabel, 'Visible', 'off');
        end
    end

    function statsButton_Callback(varargin)
        if ~isempty(tracks)
            plotTrackClasses([tracks.catIdx]);
        end
    end

    function trackSettings_Callback(varargin)
        % open window with settings panel
        tpos = get(hfig, 'Position');
        tpos = [tpos(1)+tpos(3)/2-150 tpos(2)+tpos(4)/2-75 300 245];
        pht = figure('Units', 'pixels', 'Position', tpos,...
            'PaperPositionMode', 'auto', 'Menubar', 'none', 'Toolbar', 'none',...
            'Color', get(0,'defaultUicontrolBackgroundColor'),...
            'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels',...
            'Name', 'Track display settings', 'NumberTitle', 'off');
        
        
        b = 215;
        uicontrol(pht, 'Style', 'text', 'String', 'Max. intensity threshold:',...
            'Position', [5 b 165 20], 'HorizontalAlignment', 'left');
        mitText = uicontrol(pht, 'Style', 'edit', 'String', num2str(maxIntT),...
            'Position', [170 b 60 20], 'HorizontalAlignment', 'right');
        
        % Lifetime selection sliders
        b  = 155;
        uicontrol(pht, 'Style', 'text', 'String', 'Lifetimes:',...
            'Position', [5 b+35 90 20], 'HorizontalAlignment', 'left');
        
        uicontrol(pht, 'Style', 'text', 'String', 'Min.:',...
            'Position', [5 b+18 30 20], 'HorizontalAlignment', 'left');
        minLftSlider = uicontrol(pht, 'Style', 'slider',...
            'Value', minVal, 'SliderStep', data.framerate/(maxLft-minLft-data.framerate)*[1 5], 'Min', minLft, 'Max', maxLft,...
            'Position', [40 b+20 200 18]);
        addlistener(handle(minLftSlider), 'Value', 'PostSet', @minSlider_Callback);
        minTxt = uicontrol(pht, 'Style', 'text', 'String', [num2str(minVal) ' s'],...
            'Position', [240 b+18 30 20], 'HorizontalAlignment', 'left');
        
        uicontrol(pht, 'Style', 'text', 'String', 'Max.:',...
            'Position', [5 b-2 30 20], 'HorizontalAlignment', 'left');
        maxLftSlider = uicontrol(pht, 'Style', 'slider',...
            'Value', maxVal, 'SliderStep', data.framerate/(maxLft-minLft-data.framerate)*[1 5], 'Min', minLft, 'Max', maxLft,...
            'Position', [40 b 200 18]);
        addlistener(handle(maxLftSlider), 'Value', 'PostSet', @maxSlider_Callback);
        maxTxt = uicontrol(pht, 'Style', 'text', 'String', [num2str(maxLft) ' s'],...
            'Position', [240 b-2 30 20], 'HorizontalAlignment', 'left');

        
        % Category selection buttons
        b = 115;
        catCheck = zeros(1,8);
        uicontrol(pht, 'Style', 'text', 'String', 'Single tracks: ',...
            'Position', [5 b+10 90 20], 'HorizontalAlignment', 'left');
        catCheck(1) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Valid',...
            'Position', [5 b 60 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(1));
        catCheck(2) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Faulty',...
            'Position', [65 b 140 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(2));
        catCheck(3) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Cut',...
            'Position', [125 b 80 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(3));
        catCheck(4) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Persistent',...
            'Position', [185 b 90 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(4));
        
        b = 75;
        uicontrol(pht, 'Style', 'text', 'String', 'Compound tracks: ',...
            'Position', [5 b+10 120 20], 'HorizontalAlignment', 'left');
        catCheck(5) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Valid',...
            'Position', [5 b 60 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(5));
        catCheck(6) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Faulty',...
            'Position', [65 b 140 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(6));
        catCheck(7) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Cut',...
            'Position', [125 b 80 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(7));
        catCheck(8) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Persistent',...
            'Position', [185 b 90 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(8));
        
        % EAP status selection buttons
        b = 35;
        eapCheck = zeros(1,3);
        uicontrol(pht, 'Style', 'text', 'String', 'EAP significance: ',...
            'Position', [5 b+10 150 20], 'HorizontalAlignment', 'left');
        eapCheck(1) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Independent',...
            'Position', [5 b 110 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(1));
        eapCheck(2) = uicontrol(pht, 'Style', 'checkbox', 'String', 'M/S',...
            'Position', [110 b 140 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(2));
        eapCheck(3) = uicontrol(pht, 'Style', 'checkbox', 'String', 'N.S.',...
            'Position', [165 b 80 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(3));

        
        uicontrol(pht, 'Style', 'pushbutton', 'String', 'Reset',...
            'Position', [20 5 100 20], 'HorizontalAlignment', 'left',...
            'Callback', @resetButton_Callback);
        uicontrol(pht, 'Style', 'pushbutton', 'String', 'Apply',...
            'Position', [180 5 100 20], 'HorizontalAlignment', 'left',...
            'Callback', @applyButton_Callback);
        
        
        function minSlider_Callback(~, eventdata)
            obj = get(eventdata, 'AffectedObject');
            minVal = round(get(obj, 'Value'));
            if minVal >= maxVal
                minVal = maxVal;
                set(minLftSlider, 'Value', minVal);
            end
            set(minTxt, 'String', [num2str(minVal) ' s']);
        end
        
        function maxSlider_Callback(~, eventdata)
            obj = get(eventdata, 'AffectedObject');
            maxVal = round(get(obj, 'Value'));
            if maxVal <= minVal
                maxVal = minVal;
                set(maxLftSlider, 'Value', maxVal);
            end
            set(maxTxt, 'String', [num2str(maxVal) ' s']);            
        end
        
        function resetButton_Callback(varargin)
            set([catCheck eapCheck], 'Value', true);
            maxVal = maxLft;
            minVal = data.framerate*5;%minLft;
            set(mitText, 'String', '0');
            set(maxLftSlider, 'Value', maxVal);
            set(minLftSlider, 'Value', minVal);
        end
        
        function applyButton_Callback(varargin)
            % update track selection index
            catCheckVal = cell2mat(get(catCheck, 'Value'))==1;
            eapCheckVal = cell2mat(get(eapCheck, 'Value'))==1;
            maxIntT = str2double(get(mitText, 'String'));
            
            S = [tracks.significantSlave];
            M = [tracks.significantMaster];
            % EAP: indep: M(2,:)==1; M/S M(2,:)==0 & S(2,:)==1; n.s. S(2,:)==0
            selIndex = ismember([tracks.catIdx], find(catCheckVal)) & ...
                minVal<=[tracks.lifetime_s] & [tracks.lifetime_s]<=maxVal & ...
                ((eapCheckVal(1) & M(2,:)==1) | ...
                (eapCheckVal(2) & M(2,:)==0 & S(2,:)==1) | ...
                (eapCheckVal(3) & S(2,:)==0)) & maxA(1,:)>=maxIntT;
            
            % update track selection
            tcur = find(selIndex, 1, 'first');
            
            set(handles.trackSlider, 'Min', 1);
            set(handles.trackSlider, 'Max', sum(selIndex));
            set(handles.trackSlider, 'SliderStep', [1 1]/(sum(selIndex)-1));
            set(handles.trackSlider, 'Value', 1);
            
            updateSlice();
            %updateTrack();
            close(pht);
            fprintf('# tracks selected: %d\n', sum(selIndex));
        end
        
    end

    function setColorbar(mode)        
        lfont = {'FontName', 'Helvetica', 'FontSize', 12};
        sfont = {'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'normal'};
        if ~isempty(tracks)
            switch mode
                case 'Lifetime'
                    df = 40;
                    dcoord = 0.25/df;
                    lmap = [jet(120); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', [1 20:20:120 160],...
                        'YTickLabel', [data.framerate 20:20:120 (nf-1)*data.framerate], sfont{:});
                    text(-0.1, 80, 'Lifetime (s)', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                case 'Category'
                    xlabels = {' valid', ' faulty', ' cut', ' persistent',...
                        ' valid', ' faulty', ' cut', ' persistent'};
                    lmap = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                    text(-.1, 2.5, 'Single', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                    text(-.1, 6.5, 'Compound', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                case 'EAP Status'
                    xlabels = {' N.S.', ' Signif. M/S', ' Signif. indep.'};
                    lmap = hsv2rgb([0 0 0.8; 0.55 1 0.9; 0.33 1 0.9]); % ns, slave sig., master sig.
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                case 'Object Type'
                    xlabels = {' Diff. lim.', ' Other'};
                    lmap = [0 0.8 0; 0.8 0 0];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                otherwise
                    cla(hLegend);
                    set(hLegend, 'Visible', 'off');
            end
        end
    end


    function montageButton_Callback(varargin)
        
        % Creates a montage based on the master track
        if ~isempty(tcur)
            fprintf('Generating montage...');
            if get(montageAlignCheckbox, 'Value')
                ref = 'Track';
            else
                ref = 'Frame';
            end
            [istack, xa, ya] = getTrackStack(tcur, 6, ref);
            plotTrackMontage(tracks(tcur), istack, xa, ya, 'Labels', data.markers,...
                'ShowMarkers', get(montageMarkerCheckbox, 'Value')==1,...
                'ShowDetection', get(montageDetectionCheckbox, 'Value')==1);
            fprintf(' done.\n');
        else
            fprintf('Cannot create montage: no track selected.\n');
        end
    end


    function [tstack, xa, ya] = getTrackStack(t, w, reference)
        
        sigma = frameInfo(1).s;
        w = ceil(w*sigma);
        
        % coordinate matrices
        x0 = tracks(t).x;
        y0 = tracks(t).y;
        
        % start and end buffer sizes
        if ~isempty(tracks(t).startBuffer)
            sb = numel(tracks(t).startBuffer.t);
            x0 = [tracks(t).startBuffer.x x0];
            y0 = [tracks(t).startBuffer.y y0];
        else
            sb = 0;
        end
        if ~isempty(tracks(t).endBuffer)
            eb = numel(tracks(t).endBuffer.t);
            x0 = [x0 tracks(t).endBuffer.x];
            y0 = [y0 tracks(t).endBuffer.y];
        else
            eb = 0;
        end
        
        % frame index
        tfi = tracks(t).start-sb:tracks(t).end+eb;
        tnf = length(tfi);
        
        
        if tracks(t).nSeg==1 && strcmpi(reference, 'track') % align frames to track
            xi = round(x0(handles.mCh,:));
            yi = round(y0(handles.mCh,:));
            % ensure that window falls within frame bounds
            x0 = xi - min([xi-1 w]);
            x1 = xi + min([nx-xi w]);
            y0 = yi - min([yi-1 w]);
            y1 = yi + min([ny-yi w]);
            % axes for each frame
            xa = arrayfun(@(i) x0(i):x1(i), 1:tnf, 'unif', 0);
            ya = arrayfun(@(i) y0(i):y1(i), 1:tnf, 'unif', 0);
        else
            % window around track mean
            mu_x = round(nanmean(x0,2));
            mu_y = round(nanmean(y0,2));
            x0 = max(1, min(mu_x)-w);
            x1 = min(data.imagesize(2), max(mu_x)+w);
            y0 = max(1, min(mu_y)-w);
            y1 = min(data.imagesize(1), max(mu_y)+w);
            xa = repmat({x0:x1}, [tnf 1]);
            ya = repmat({y0:y1}, [tnf 1]);
        end
        
        tstack = cell(nCh,tnf);
        for ci = 1:nCh
            for k = 1:tnf
                tstack{ci,k} = stack{ci}(ya{k}, xa{k}, tfi(k));
            end
        end
    end


    function trackButton_Callback(varargin)
        [x0,y0] = ginput(1);
        ci = find(handles.fAxes(:,1)==gca, 1);
        if ~isempty(ci) && ~isempty(tracks)
            % track segments visible in current frame
            cidx = find([tracks.start]<=fidx & fidx<=[tracks.end] & selIndex);
            if ~isempty(cidx)
                % distance to mean of tracks
                d = sqrt((x0-mu_x(cidx)).^2 + (y0-mu_y(cidx)).^2);
                [~,d] = nanmin(d);
                tcur = cidx(d);
                set(handles.trackSlider, 'Value', find(find(selIndex)==tcur)); % calls updateTrack
            end
        end
    end


    function printButton_Callback(varargin)
        fprintf('Printing figures ...');
        
        % Tracks
        if ~isempty(tracks)
            for ch = 1:nCh
                plotTrack(data, tracks(tcur), ch,...
                    'FileName', ['track_' num2str(tcur) '_ch' num2str(ch) '.eps'],...
                    'Visible', 'off', 'DisplayMode', 'Print');
            end
            
            if get(montageAlignCheckbox, 'Value')
                ref = 'Track';
            else
                ref = 'Frame';
            end
            [tstack, xa, ya] = getTrackStack(tcur, 6, ref);
            fpath = [data.source 'Figures' filesep 'track_' num2str(tcur) '_montage.eps'];
                plotTrackMontage(tracks(tcur), tstack, xa, ya, 'Labels', data.markers,...
                    'Visible', 'off', 'epsPath', fpath,...
                    'ShowMarkers', get(montageMarkerCheckbox, 'Value')==1,...
                    'ShowDetection', get(montageDetectionCheckbox, 'Value')==1);
        end
        
        % Frames
        if strcmp(displayType, 'RGB')
            maxCh = 1;
            chLabel = {'RGB'};
        else
            maxCh = nCh;
            chLabel = arrayfun(@(i) ['ch' num2str(i)], 1:nCh, 'unif', 0);
        end
                
        f0 = figure('PaperPositionMode', 'auto', 'Position', [20 20 nx ny], 'Visible', 'off',...
            'DefaultLineLineSmoothing', 'on', 'DefaultPatchLineSmoothing', 'on');
        colormap(gray(256));
        fpath = [data.source 'Figures' filesep];
        for ci = 1:maxCh
            h0 = copyobj(handles.fAxes(ci,1),f0);
            hx = findobj(h0, 'DisplayName', 'FrameMarker');
            delete(hx);
            hx = findobj(h0, 'DisplayName', 'TrackMarker');
            delete(hx);
            hx = findobj(h0, 'LineStyle', '-');
            set(hx, 'LineWidth', 1);
            
            hx = findobj(h0, 'Marker', 'o');
            set(hx, 'MarkerSize', 9, 'LineWidth', 1);
            
            set(h0, 'Position', [0 0 1 1]);
            print(f0, '-depsc2', '-loose', [fpath 'frame_' num2str(fidx) '_' chLabel{ci} '.eps']);
            %print(f0, '-dpng', '-loose', [fpath 'frame_' num2str(fidx) '_' chLabel{ci} '.png']);
            delete(h0);
        end
        close(f0);
        
        fprintf([' done. Figures saved in ' getShortPath(data) 'Figures.\n']);
    end


    function movieButton_Callback(varargin)
        
        fopts = {'Visible', 'off', 'Position', [20 20 nx ny],...
            'InvertHardcopy', 'off', 'PaperUnits', 'Points', 'PaperSize', [nx ny],...
            'PaperPosition', [0 0 nx ny], 'PaperPositionMode', 'auto',...
            'DefaultLineLineSmoothing','on', 'DefaultPatchLineSmoothing','on'};
        
        if strcmp(displayType, 'RGB')
            maxCh = 1;
        else
            maxCh = nCh;
        end
        
        mpath = [data.source 'Movies' filesep];
        fpath = [mpath 'Frames' filesep];
        [~,~] = mkdir(mpath);
        [~,~] = mkdir(fpath);
        
        fmt = ['%0' num2str(ceil(log10(nf))) 'd'];
        
        f0 = figure(fopts{:});
        colormap(gray(256));
        ha = axes('Position', [0 0 1 1]);
        for ci = 1:maxCh
            fprintf('Generating movie frames:     ');
            for fi = 1:nf
                switch displayType
                    case 'raw'
                        imagesc(stack{ci}(:,:,fi), 'Parent', ha);
                    case 'mask'
                        if ci==1
                            imagesc(rgbOverlay(stack{1}(:,:,fi), dmask(:,:,fi), [1 0 0], dRange{1}), 'Parent', ha);
                        else
                            imagesc(stack{ci}(:,:,fi), 'Parent', ha);
                        end
                    case 'RGB'
                        rframe = zeros(ny,nx,3,'uint8');
                        idxRGB = getRGBindex(data.markers);
                        for c2 = 1:nCh
                            rframe(:,:,idxRGB(c2)) = uint8(scaleContrast(double(stack{c2}(:,:,fi)), dRange{c2}));
                        end
                        imagesc(rframe, 'Parent', ha);
                end
                hold(ha, 'on');
                caxis(ha, dRange{ci});
                
                if ~isempty(tracks) && fi~=1 && get(trackCheckbox, 'Value')
                    vidx = ~isnan(X(fi,:));
                    
                    set(ha, 'ColorOrder', cmap(tstruct.idx(vidx),:));
                    
                    plot(ha, X(1:fi,vidx), Y(1:fi,vidx), 'HitTest', 'off');
                    if get(gapCheckbox, 'Value')
                        hpg = plot(ha, X(fi,vidx & G(fi,:)), Y(fi,vidx & G(fi,:)), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
                    end
                    if get(trackEventCheckbox, 'Value')
                        % Births
                        bcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackStarts==fi), 'unif', 0);
                        bcoord = vertcat(bcoord{:});
                        plot(ha, bcoord(:,1), bcoord(:,2), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                        
                        % Deaths
                        dcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackEnds==fi), 'unif', 0);
                        dcoord = vertcat(dcoord{:});
                        plot(ha, dcoord(:,1), dcoord(:,2), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);
                    end
                    
                end
                if ~isempty(tracks) && get(eapCheckbox, 'Value') && ci>1
                    sel = fvec==fi & mvec(ci,:)==1;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([1/3 1 0.9]), 'MarkerSize', 8);
                    sel = fvec==fi & mvec(ci,:)==0 & svec(ci,:)==1;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([0.55 1 0.9]), 'MarkerSize', 8);
                    sel = fvec==fi & mvec(ci,:)==0 & svec(ci,:)==0;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', 0.8*[1 1 1], 'MarkerSize', 8);
                end
                
                if get(detectionCheckbox, 'Value') && ~isempty(frameInfo)
                    isPSF = frameInfo(fi).isPSF(1,:)==1;
                    if any(isPSF)
                        plot(ha, frameInfo(fi).x(1,isPSF), frameInfo(fi).y(1,isPSF), 'o', 'Color', [0 0.6 0], 'MarkerSize', 8);
                    end
                    if any(~isPSF)
                        plot(ha, frameInfo(fi).x(1,~isPSF), frameInfo(fi).y(1,~isPSF), 'o', 'Color', [0.6 0 0], 'MarkerSize', 8);
                    end
                end
                
                axis(ha, 'off');
                print(f0, '-dpng', '-loose', ['-r' num2str(1*72)], [fpath 'frame' num2str(fi, fmt) '_ch' num2str(ci) '.png']);
                %print(h, '-djpeg100', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
                cla(ha);
                fprintf('\b\b\b\b%3d%%', round(100*fi/nf));

            end
            fprintf('\n');
        end
        fprintf(['Frames saved to ' getShortPath(data) 'Movies' filesep 'Frames.\n']);
        close(f0);
        
        % side-by-side frame arrangement
        for fi = 1:nf
            % channel frames
            cpath = arrayfun(@(ci) [fpath 'frame' num2str(fi, fmt) '_ch' num2str(ci) '.png '], 1:nCh, 'unif', 0);
            fname = [fpath 'montage' num2str(fi, fmt) '.png'];
            
            cmd = ['export DYLD_LIBRARY_PATH=""; montage -geometry +3+3+0+0 -background "rgb(255,255,255)" '...
                [cpath{:}] ' -compress lzw ' fname];
            system(cmd);
            cmd = ['export DYLD_LIBRARY_PATH=""; convert ' fname ' -shave 3x3 -depth 8 ' fname];
            system(cmd);
        end
        
        % Generate movie, if on a unix system with ffmpeg
        if isunix && ~system('which ffmpeg >/dev/null 2>&1')
            fprintf('Generating movie ... ');
            %fr = num2str(framerate);
            fr = num2str(15);           
            
            %cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'frame' fmt '_ch' num2str(1) '.png' ' -vf "scale=' num2str(2*floor(nx/2)) ':' num2str(2*floor(ny/2))...
            %    '" -c:v libx264 -crf 22 -pix_fmt yuv420p ' mpath 'Movie_ch1.mp4'];
            cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'montage' fmt '.png' ' -vf "scale=' num2str(2*floor((nCh*nx+(nCh-1)*6)/2)) ':' num2str(2*floor(ny/2))...
                '" -c:v libx264 -crf 22 -pix_fmt yuv420p ' mpath getCellDir(data) '.mp4'];
            system(cmd);
            
            fprintf(' done.\n');
        else
            fprintf('A unix system with ffmpeg installed is required to generate movies automatically.\n');
        end
    end

end



function figResize(src,~)
handles = getappdata(src, 'handles');

pos = get(src, 'Position');


set(handles.frameLabel, 'Position', [20 pos(4)-20, 100 15]);

% tracks
set(handles.tplotPanel, 'Position', [pos(3)-535 5 160 70]);
set(handles.montagePanel, 'Position', [pos(3)-370 5 220 70]);
set(handles.outputPanel, 'Position', [pos(3)-145 5 140 70]);

% spacers:
tspace = 20;
bspace = 100;
lspace = 10;
rspace = 400;
spacer = 10; % space between panels

width = pos(3) - rspace - lspace;
height = pos(4) - bspace - tspace;

set(handles.frameSlider, 'Position', [lspace 77 pos(3)-rspace-lspace 18]);

switch numel(handles.fPanels)
    case 1
        set(handles.fPanels(1), 'Position', [lspace bspace width height]);
    case 2
        %if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            width = (width-spacer)/2;
            set(handles.fPanels(1), 'Position', [lspace bspace width height]);
            set(handles.fPanels(2), 'Position', [lspace+width+spacer bspace width height]);
%         else % vertical
%             height = (height-spacer)/2;
%             set(handles.fPanels(1), 'Position', [lspace bspace+spacer+height width height]);
%             set(handles.fPanels(2), 'Position', [lspace bspace width height]);
%         end
    case 3
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        set(handles.fPanels(1), 'Position', [lspace bspace+spacer+height width height]); % top left
        set(handles.fPanels(2), 'Position', [lspace+width+spacer bspace+height+spacer width height]); % top right
        set(handles.fPanels(3), 'Position', [lspace bspace width height]); % bottom left
    case 4
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        set(handles.fPanels(1), 'Position', [lspace bspace+spacer+height width height]); % top left
        set(handles.fPanels(2), 'Position', [lspace+width+spacer bspace+height+spacer width height]); % top right
        set(handles.fPanels(3), 'Position', [lspace bspace width height]); % bottom left
        set(handles.fPanels(4), 'Position', [lspace+width+spacer bspace width height]); % bottom right
end

spacer = 15;
w = 320;
nCh = numel(handles.tAxes);
h_tot = pos(4) - 140;
h = min((h_tot-(nCh-1)*spacer)/nCh, 200);
dx = pos(3)-w-30;
switch nCh
    case 1
        set(handles.tAxes(1), 'Position', [dx 120+(h_tot-h) w h]);
    case 2
        set(handles.tAxes(1), 'Position', [dx 120+(h_tot-h) w h]);
        set(handles.tAxes(2), 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
    case 3
        set(handles.tAxes(1), 'Position', [dx 120+(h_tot-h) w h]);
        set(handles.tAxes(2), 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
        set(handles.tAxes(3), 'Position', [dx 120+(h_tot-3*h-2*spacer) w h]);
    case 4
        set(handles.tAxes(1), 'Position', [dx 120+(h_tot-h) w h]);
        set(handles.tAxes(2), 'Position', [dx 120+(h_tot-2*h-spacer) w h]);
        set(handles.tAxes(3), 'Position', [dx 120+(h_tot-3*h-2*spacer) w h]);
        set(handles.tAxes(4), 'Position', [dx 120+(h_tot-4*h-3*spacer) w h]);
end
set(handles.trackLabel, 'Position', [pos(3)-70 pos(4)-20 100 15]);
set(handles.trackSlider, 'Position', [pos(3)-24 120 18 h_tot]);

end



function handles = setupFrameAxes(hfig, spos, N)

handles = getappdata(hfig, 'handles');
if nargin<3
    N = handles.nCh;
end

pos = get(gcf, 'Position'); % [pixels]

% spacers: 
lspace = spos(1);
bspace = spos(2);
rspace = spos(3);
tspace = spos(4);
spacer = spos(5); % space between panels

width = pos(3) - rspace - lspace;
height = pos(4) - bspace - tspace;

% reset axes etc.
if isfield(handles, 'fPanels') && ~isempty(handles.fPanels)
    delete(handles.fPanels);
end
handles.fPanels = zeros(1,N);
uiOpts = {'Parent', hfig, 'Units', 'pixels', 'BorderType', 'none'};
switch N
    case 1
        handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]);
    case 2
        %if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            width = (width-spacer)/2;
            handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]);
            handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace+width+spacer bspace width height]);
%         else % vertical
%            height = (height-spacer)/2;
%            handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace+spacer+height width height]);
%            handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]);
%         end
    case 3
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace+spacer+height width height]); % top left
        handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace+width+spacer bspace+height+spacer width height]); % top right
        handles.fPanels(3) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]); % bottom left
    case 4
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace+spacer+height width height]); % top left
        handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace+width+spacer bspace+height+spacer width height]); % top right
        handles.fPanels(3) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]); % bottom left
        handles.fPanels(4) = uipanel(uiOpts{:}, 'Position', [lspace+width+spacer bspace width height]); % bottom right
end
setappdata(hfig, 'handles', handles);
end



function [ha, hl] = setupStackViewer(hf, dims, addLegend)
if nargin<3
    addLegend = false;
end

spc = 6; % spacer, fixed [pixels]

nx = dims(1);
ny = dims(2);
nz = dims(3);
pos = get(hf, 'Position');
w = pos(3);
h = pos(4);

% normalized axes dimensions
fx = (w-spc)/(nx+nz);
fy = (h-spc)/(ny+nz);
f = min(fx,fy);
h = (ny+nz)*f+spc; % [pixels]
w = (nx+nz)*f+spc;

rxy = pos(3)/pos(4);
dx = spc/pos(3);
dy = spc/pos(4);
if rxy > w/h % figure is too wide
    f0 = w/h / rxy;
    left = (1-f0)/2;
    ha(1) = axes('Position', [left+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h], 'Parent', hf);
    ha(2) = axes('Position', [left 0 f0*f*nz/w f*ny/h], 'Parent', hf); % bottom left
    ha(3) = axes('Position', [left+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h], 'Parent', hf);
else
    f0 = h/w * rxy;
    left = 0;
    ha(1) = axes('Position', [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h], 'Parent', hf);
    ha(2) = axes('Position', [0 1-f0 f*nz/w f0*f*ny/h], 'Parent', hf);
    ha(3) = axes('Position', [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h], 'Parent', hf);
end
if addLegend
    lpos = get(ha(3), 'Position');
    lpos([1 3]) = [left+15/pos(3) 15/pos(3)];
    hl = axes('Position', lpos, 'Parent', hf);
else
    hl = NaN;
end

set(hf, 'ResizeFcn', @pResize);

    function pResize(~,~)
        ipos = get(hf, 'Position');
        rxy = ipos(3)/ipos(4);
        dx = spc/ipos(3);
        dy = spc/ipos(4);
        if rxy > w/h % figure is too wide
            f0 = w/h / rxy;
            left = (1-f0)/2;
            set(ha(1), 'Position', [left+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h]);
            set(ha(2), 'Position', [left 0 f0*f*nz/w f*ny/h]);
            set(ha(3), 'Position', [left+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h]);
        else
            f0 = h/w * rxy;
            left = 0;
            set(ha(1), 'Position', [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h]);
            set(ha(2), 'Position', [0 1-f0 f*nz/w f0*f*ny/h]);
            set(ha(3), 'Position', [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h]);
        end
        if ~isnan(hl)
            lpos = get(ha(3), 'Position');
            lpos([1 3]) = [left 15/ipos(3)];
            set(hl, 'Position', lpos);
        end
    end
end
