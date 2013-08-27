%cmeDataViewer(data, varargin) displays movies with associated detection and tracking results.
%
% Inputs:    
%             data : single movie structure returned by loadConditionData.m
%     Trajectories : optional input for selecting 'all' (default) or
%                    'valid' CCS trajectories.
%
% Notes: Only tracks with at least 5 frames are loaded and displayed.

% Francois Aguet, 2011 (last modified 08/24/2013)

function cmeDataViewer2(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('Trajectories', 'all', @(x) isempty(x) || isstruct(x) || any(strcmpi(x, {'all', 'valid'})));
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

%---------------------
% Load movie
%---------------------
handles.stack = cell(1,nCh);
% handles.stack{1} = rand(150,150,25);
% handles.stack{2} = rand(150,150,25);
if ~iscell(data.framePaths{1})
    for c = 1:nCh
        %stack{c} = readtiff(data.framePaths{c});
        handles.stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
        for i = 1:data.movieLength
            handles.stack{c}(:,:,i) = imread(data.framePaths{c}, i);
        end
    end
else
    for c = 1:nCh
        handles.stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
        for i = 1:data.movieLength
            handles.stack{c}(:,:,i) = imread(data.framePaths{c}{i});
        end
    end
end
[ny,nx,nf] = size(handles.stack{1});

%--------------------------------
% Load detection files
%--------------------------------
handles.detection = cell(1,nCh);
for c = 1:nCh
    detectionFile = [data.channels{c} 'Detection' filesep 'detection_v2.mat'];
    if (exist(detectionFile, 'file')==2)
%         load(detectionFile);
%         handles.detection{c} = frameInfo;
    end
end

handles.f = 1;


%---------------------
% Load tracks
%---------------------
tracks = loadTracks(data, 'Category', 'all', 'Mask', false, 'Cutoff_f', 5);
nt = numel(tracks);
nseg = [tracks.nSeg];

np = sum(nseg);
X = NaN(nf, np);
Y = NaN(nf, nt);
% vector of start indexes since multiple segments/track
tidx = cumsum([1 nseg(1:end-1)]);

for t = 1:nt
    if nseg(t)==1
        X(tracks(t).f, tidx(t)) = tracks(t).x(1,:);
        Y(tracks(t).f, tidx(t)) = tracks(t).y(1,:);
    else
        sep = find(isnan(tracks(t).t));
        sep = [0 sep numel(tracks(t).f)+1]; %#ok<AGROW>
        for s = 1:tracks(t).nSeg
            sidx = sep(s)+1:sep(s+1)-1;
            X(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).x(1,sidx);
            Y(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).y(1,sidx);
        end
    end
end

tstruct.X = X;
tstruct.Y = Y;
% index 'label' mask
% [1 1 2 3 4 4 ... ] first two cols are from same track
idx = diff([tidx size(X,2)+1]);
idx = arrayfun(@(i) i+zeros(1, idx(i)), 1:numel(idx), 'unif', 0);
tstruct.idx = [idx{:}];
tstruct.n = numel(tracks);

hpt = []; % handles for track plot objects
cmap = hsv2rgb([rand(tstruct.n,1) ones(tstruct.n,2)]);


handles.tracks = cell(1,nCh);
% handles.colorMap = cell(1,nCh);
% handles.maxLifetime_f = [];
% handles.selectedTrack = [];
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
%     
%     if ~isempty(tracks)
%         handles.tracks{handles.mCh} = tracks;
%         
%         if exist([data.source 'Analysis' filesep 'BackgroundFits.mat'],'file')==2
%             load([data.source 'Analysis' filesep 'BackgroundFits.mat']);
%             handles.bg95 = bg95;
%         else
%             handles.bg95 = [];
%         end
%         
%         % random colormaps for track display
%         for c = 1:nCh
%             nt = numel(handles.tracks{c});
%             handles.colorMap{c} = hsv2rgb([rand(nt,1) ones(nt,2)]);
%         end
%         
%         handles.maxLifetime_f = max([handles.tracks{handles.mCh}.end]-[handles.tracks{handles.mCh}.start]+1);
%         handles.selectedTrack = NaN(1,handles.nCh);
%         handles.selectedTrack(handles.mCh) = 1;
%         handles.f = handles.tracks{handles.mCh}(1).start;
%         
%         % min/max track intensities
%         maxA = arrayfun(@(t) max(t.A, [], 2), handles.tracks{1}, 'UniformOutput', false);
%         maxA = [maxA{:}];
%         handles.maxA = zeros(1,nCh);
%         for c = 1:nCh
%             [f_ecdf, x_ecdf] = ecdf(maxA(c,:));
%             handles.maxA(c) = interp1(f_ecdf, x_ecdf, 0.975);
%         end
%         d = floor(log10(handles.maxA));
%         % y-axis unit
%         handles.yunit = round(handles.maxA ./ 10.^d) .* 10.^(d-1);
%         handles.maxA = ceil(handles.maxA ./ handles.yunit) .* handles.yunit;
%     else
%         handles.maxLifetime_f = [];
%         handles.selectedTrack = [];
%         handles.f = 1;
%     end   
% end
% handles.displayType = 'raw';
% handles.pUnitType = 's';
    

%---------------------
% Setup main window
%---------------------
hfig = figure('Units', 'normalized', 'Position', [0.05 0.2 0.9 0.8],...
    'PaperPositionMode', 'auto', 'Toolbar', 'figure',...
    'Color', get(0,'defaultUicontrolBackgroundColor'),...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', getCellDir(data));

pos = get(hfig, 'Position');

% fixed width of the track plots, in pixels
% w = 350;
% dx = pos(3)-w-50; % space available for frame display

%---------------------
% Frames
%---------------------
% width = pos(3) - w-50-100-50 -50;

handles.frameLabel = uicontrol('Style', 'text', 'String', ['Frame ' num2str(handles.f)], ...
    'Position', [20 pos(4)-20 100 15], 'HorizontalAlignment', 'left');

% Frame slider
% if data.movieLength>1
%     handles.frameSlider = uicontrol('Style', 'slider', 'Units', 'pixels',...
%         'Value', handles.f, 'SliderStep', [1/(data.movieLength-1) 0.05], 'Min', 1, 'Max', data.movieLength,...
%         'Position', [20 60 width 20], 'Callback', {@frameSlider_Callback, hfig});
% end
    
% Main control panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', '', 'Position', [5 5 650 70]);

uicontrol(ph, 'Style', 'text', 'String', 'Display: ',...
    'Position', [5 40 60 20], 'HorizontalAlignment', 'left');
handles.frameChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Raw frames', 'Detection', 'RGB'},...
    'Position', [65 42 120 20], 'Callback', {@frameChoice_Callback, hfig});

handles.detectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Detections',...
    'Position', [200 45 100 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.trackCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Tracks:', 'Value', true,...
    'Position', [200 25 80 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.trackChoice = uicontrol('Style', 'popup',...
    'String', {'Category', 'Lifetime', 'Object Type', 'Random'},...
    'Position', [280 28 100 20], 'Callback', {@trackChoice_Callback, hfig});

handles.gapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Gaps',...
    'Position', [390 45 140 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.trackEventCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Births/Deaths',...
    'Position', [390 25 140 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.eapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'EAP status',...
    'Position', [390 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});

handles.labelCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Channel labels',...
    'Position', [200 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});


handles.trackButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Select track',...
    'Position', [540 40 100 20], 'HorizontalAlignment', 'left',...
    'Callback', {@trackButton_Callback, hfig});
handles.statsButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Track statistics',...
    'Position', [540 10 100 20], 'HorizontalAlignment', 'left',...
    'Callback', {@statsButton_Callback, hfig});


%---------------------
% Tracks
%---------------------
handles.trackLabel = uicontrol('Style', 'text', 'String', 'Track 1',...
    'Units', 'pixels', 'Position', [pos(3)-400 pos(4)-20 100 15], 'HorizontalAlignment', 'left');

handles.trackSlider = uicontrol('Style', 'slider',...
    'Value', 1, 'SliderStep', [1 1], 'Min', 1, 'Max', 100,...
    'Position', [pos(3)-35 110 20 pos(4)-130],...
    'Callback', {@trackSlider_Callback, hfig});

% Track plot panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Plot options', 'Position', [pos(3)-540 5 160 70]);

uicontrol(ph, 'Style', 'text', 'String', 'Units: ',...
    'Position', [5 35 60 20], 'HorizontalAlignment', 'left');
handles.tplotUnitChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Seconds', 'Frames'},...
    'Position', [40 40 100 15], 'Callback', {@unitChoice_Callback, hfig});
handles.tplotBackgroundCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Subtract background',...
    'Position', [5 20 120 15], 'HorizontalAlignment', 'left', 'Value', true, 'Callback', {@refreshTracks_Callback, hfig});
handles.tplotScaleCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Autoscale',...
    'Position', [5 5 120 15], 'HorizontalAlignment', 'left', 'Value', false, 'Callback', {@refreshTracks_Callback, hfig});
handles.tplotPanel = ph;


% Montage panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Montage', 'Position', [pos(3)-400 5 200 70]);
handles.montageAlignCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Align to track',...
    'Position', [5 35 120 15], 'HorizontalAlignment', 'left', 'Value', true);
handles.montageMarkerCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show markers',...
    'Position', [5 20 120 15], 'HorizontalAlignment', 'left');
handles.montageDetectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show detection',...
    'Position', [5 5 120 15], 'HorizontalAlignment', 'left');
handles.montageButton = uicontrol(ph, 'Style', 'pushbutton','String','Generate',...
    'Units', 'pixels', 'Position', [120 15 75 30],...%[.1 .55 .6 .4]
    'Callback', {@montageButton_Callback, hfig});
handles.montagePanel = ph;


% Output panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Output', 'Position', [pos(3)-180 5 140 70]);

handles.printButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Print figures',...
    'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.45],...
    'Callback', {@printButton_Callback, hfig});

handles.movieButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Make movie',...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.45],...
    'Callback', {@movieButton_Callback, hfig});
handles.outputPanel = ph;

setappdata(hfig, 'handles', handles);

%================================


% handles.fAspectRatio = handles.data.imagesize(1) / handles.data.imagesize(2);



% dynamic range for each channel
handles.dRange = cell(1,nCh);
for c = 1:nCh
    handles.dRange{c} = [min(handles.stack{c}(:)) max(handles.stack{c}(:))];
end

% initialize handles
handles.trackMode = 'Category';
handles.hues = getFluorophoreHues(data.markers);
handles.rgbColors = arrayfun(@(x) hsv2rgb([x 1 1]), handles.hues, 'UniformOutput', false);

settings.zoom = 1;
setappdata(hfig, 'settings', settings);


%=================================================
% Set initial values for sliders and checkboxes
%=================================================
if ~isempty([handles.tracks{handles.mCh}]) && length(handles.tracks{handles.mCh}) > 1
    set(handles.trackSlider, 'Min', 1);
    nTracks = length(handles.tracks{handles.mCh});
    set(handles.trackSlider, 'Max', nTracks);
    set(handles.trackSlider, 'SliderStep', [1/(nTracks-1) 0.05]);
else
    set(handles.trackSlider, 'Visible', 'off');
end

if nCh > 2
    set(handles.('labelCheckbox'), 'Value', 1);
end




% %=================================================
% % Generate axes
% %=================================================
% % track panels: 20 spacer, 110 bottom, 30 top
% h_tot = pos(4) - 140;
% h = min((h_tot-(nCh-1)*20)/nCh, 200);
% 
% opts = {'Parent', gcf, 'Units', 'pixels', 'Box', 'on'};
% 
% switch nCh
%     case 1
%         handles.tAxes(1) = axes(opts{:}, 'Position', [dx 110+(h_tot-h) w h]);
%     case 2
%         handles.tAxes(1) = axes(opts{:}, 'Position', [dx 110+(h_tot-h) w h]);
%         handles.tAxes(2) = axes(opts{:}, 'Position', [dx 110+(h_tot-2*h-20) w h]);
%     case 3
%         handles.tAxes(1) = axes(opts{:}, 'Position', [dx 110+(h_tot-h) w h]);
%         handles.tAxes(2) = axes(opts{:}, 'Position', [dx 110+(h_tot-2*h-20) w h]);
%         handles.tAxes(3) = axes(opts{:}, 'Position', [dx 110+(h_tot-3*h-40) w h]);
%     case 4        
%         handles.tAxes(1) = axes(opts{:}, 'Position', [dx 110+(h_tot-h) w h]);
%         handles.tAxes(2) = axes(opts{:}, 'Position', [dx 110+(h_tot-2*h-20) w h]);
%         handles.tAxes(3) = axes(opts{:}, 'Position', [dx 110+(h_tot-2*h-40) w h]);
%         handles.tAxes(4) = axes(opts{:}, 'Position', [dx 110+(h_tot-2*h-60) w h]);
% end
% xlabel('Time (s)');
% 
% % Legend (track lifetime/category)
% handles.cAxes = axes('Parent', gcf, 'Units', 'pixels', 'Position', [dx-100 pos(4)-230 15 200], 'Visible', 'on');

set(hfig, 'ResizeFcn', @figResize);

% % save XLim diff. for zoom reference
% handles.refXLimDiff = data.imagesize(2)-1;

setappdata(hfig, 'handles', handles); % write 'handles' to hfig
handles = setupFrameAxes(hfig);

% temp fix: set up stakviewer axes here
handles.fAxes = zeros(nCh,3);
for c = 1:nCh
    handles.fAxes(c,:) = setupStackViewer(handles.fPanels(c), [nx ny nf]); 
end
colormap(gray(256));
% populate with data, plotting functions are called only here, afterwards change data
x = round(nx/2);
y = round(ny/2);


for c = 1:nCh
    % x,y view
    hxy(c) = imagesc(handles.stack{c}(:,:,handles.f), 'Parent', handles.fAxes(c,1), 'HitTest', 'off');
    hold(handles.fAxes(c,1), 'on');
    set(handles.fAxes(c,1), 'ButtonDownFcn', @click_Callback);
    hl(c,1) = plot(handles.fAxes(c,1), [x x], [0.5 ny+0.5], 'r', 'HitTest', 'off');
    hl(c,2) = plot(handles.fAxes(c,1), [0.5 nx+0.5], [y y], 'r', 'HitTest', 'off');
    
    % y,z view
    hyz(c) = imagesc(squeeze(handles.stack{c}(:,x,:)), 'Parent', handles.fAxes(c,2), 'HitTest', 'off');
    hold(handles.fAxes(c,2), 'on');
    % line in y,z view
    hl(c,3) = plot(handles.fAxes(c,2), handles.f*[1 1], [0.5 ny+0.5], 'r');
    
    % x,z view
    hxz(c) = imagesc(squeeze(handles.stack{c}(y,:,:))', 'Parent', handles.fAxes(c,3), 'HitTest', 'off');
    hold(handles.fAxes(c,3), 'on');
    % line in x,z view
    hl(c,4) = plot(handles.fAxes(c,3), [0.5 nx+0.5], handles.f*[1 1], 'r');
    
    arrayfun(@(i) caxis(i, handles.dRange{c}), handles.fAxes(c,:), 'unif', 0);
end
set(handles.fAxes, 'XTick', [], 'YTick', []);


set(hfig, 'WindowScrollWheelFcn', @scroll_Callback);
set(hfig, 'KeyPressFcn', @key_Callback);


set(zoom, 'ActionPostCallback', @czoom);
hpan = pan;
set(hpan,'ActionPreCallback',@panstart);
set(hpan,'ActionPostCallback',@panstop);

setappdata(hfig, 'handles', handles); % messy

%===========================
% initialize figures/plots
%===========================
% linkaxes(handles.tAxes, 'x'); % calls resize??
% 
% set(hfig, 'KeyPressFcn', @keyListener);
% 
% refreshFrameDisplay(hfig);
% refreshTrackDisplay(hfig);
% 
% setColorbar(hfig, handles.trackMode);
% 
% set(zoom, 'ActionPostCallback', {@zoompostcallback, hfig});

    function click_Callback(varargin)
        updateProj(); % when clicking w/o dragging
        set(gcf, 'WindowButtonMotionFcn', @drag, 'WindowButtonUpFcn', @stopDragging);
    end

    function drag(varargin)
        updateProj();
    end

    function stopDragging(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end


    % scroll through stack slices
    function scroll_Callback(src, eventdata)
        hi = getappdata(src, 'handles');
        if eventdata.VerticalScrollCount < 0
            if hi.f < nf
                hi.f = hi.f + 1;
                setappdata(hfig, 'handles', hi);
                updateSlice();
            end
        elseif eventdata.VerticalScrollCount > 0
            if hi.f > 1
                hi.f = hi.f - 1;
                setappdata(hfig, 'handles', hi);
                updateSlice();
            end
        end
    end


    function key_Callback(~, eventdata)
        hi = getappdata(hfig, 'handles');
        switch eventdata.Key
            case 'uparrow'
                if hi.f > 1
                    hi.f = hi.f - 1;
                    setappdata(hfig, 'handles', hi);
                    updateSlice();
                end
            case 'downarrow'
                if hi.f < nf
                    hi.f = hi.f + 1;
                    setappdata(hfig, 'handles', hi);
                    updateSlice();
                end
        end
    end


    function updateSlice()
        hi = getappdata(hfig, 'handles');
        for c = 1:nCh
            set(hxy(c), 'CData', hi.stack{c}(:,:,hi.f));
        end
        set(hl(:,3), 'XData', hi.f*[1 1]);
        set(hl(:,4), 'YData', hi.f*[1 1]);
        %set(ht, 'String', ['Slice ' num2str(slice)]);
        
        if ~isempty(tstruct)
            vidx = ~isnan(tstruct.X(hi.f,:));
            if hi.f~=1
                delete(hpt);
                set(hi.fAxes(1,1), 'ColorOrder', cmap(tstruct.idx(vidx),:));
                hpt = plot(hi.fAxes(1,1), tstruct.X(1:hi.f,vidx), tstruct.Y(1:hi.f,vidx), 'HitTest', 'off');
            else
                delete(hpt);
                hpt = [];
            end
        end
    end

    function updateProj()
        a = get(gca,'CurrentPoint');
        % plot lines
        set(hl(:,1), 'XData', a(1,1)*[1 1]);
        set(hl(:,2), 'YData', a(1,2)*[1 1]);
        
         % update data
        xi = min(max(round(a(1,1)),1), nx);
        yi = min(max(round(a(1,2)),1), ny);
        for c = 1:nCh
            set(hyz(c), 'CData', squeeze(handles.stack{c}(:,xi,:)));
            set(hxz(c), 'CData', squeeze(handles.stack{c}(yi,:,:))');
        end
    end

    function czoom(~, eventdata)
        % identify panel
        ci = handles.fAxes(:,1) == eventdata.Axes;
        if any(ci) && nCh>1 % x,y axes zoomed
            XLim = get(handles.fAxes(ci,1), 'XLim');
            YLim = get(handles.fAxes(ci,1), 'YLim');
            set(handles.fAxes(~ci,1), 'XLim', XLim, 'YLim', YLim);
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


end


%===================================
% Automatic actions after zoom
%===================================
function zoompostcallback(~, eventdata, hfig)
% 
% settings = getappdata(hfig, 'settings');
% handles = getappdata(hfig, 'handles');
% 
% if ismember(eventdata.Axes, handles.fAxes)
%     settings.zoom = handles.refXLimDiff / diff(get(eventdata.Axes, 'XLim'));
%     for c = 1:length(settings.selectedTrackMarkerID)
%         id = settings.selectedTrackMarkerID(c);
%         if ~isnan(id)
%             set(id, 'MarkerSize', 10*settings.zoom);
%         end
%     end
%     setappdata(hfig, 'settings', settings);
% end
end


function figResize(src,~)
pos = get(src, 'Position');
%w = 350;
w = 10;
dx = pos(3)-w-50;
handles = getappdata(src, 'handles');

set(handles.frameLabel, 'Position', [20 pos(4)-20, 100 15]);

% tracks
set(handles.trackLabel, 'Position', [dx pos(4)-20, 100 15]);
set(handles.trackSlider, 'Position', [pos(3)-35 110 20 pos(4)-140]);
set(handles.tplotPanel, 'Position', [pos(3)-540 5 160 70]);
set(handles.outputPanel, 'Position', [pos(3)-160 5 140 70]);
set(handles.montagePanel, 'Position', [pos(3)-370 5 200 70]);

% spacers: 
tspace = 20;
bspace = 100;
lspace = 10;
rspace = 20; %%%%%%%%%%%%%%
spacer = 10; % space between panels

width = pos(3) - rspace - lspace;
height = pos(4) - bspace - tspace;
switch numel(handles.fPanels)
    case 1
        set(handles.fPanels(1), 'Position', [lspace bspace width height]);
    case 2
        if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            width = (width-spacer)/2;
            set(handles.fPanels(1), 'Position', [lspace bspace width height]);
            set(handles.fPanels(2), 'Position', [lspace+width+spacer bspace width height]);
        else % vertical
           height = (height-spacer)/2;
           set(handles.fPanels(1), 'Position', [lspace bspace+spacer+height width height]);
           set(handles.fPanels(2), 'Position', [lspace bspace width height]);
        end
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
end
% h_tot = pos(4) - 140;
% nx = numel(handles.tAxes);
% h = min((h_tot-(nx-1)*20)/nx, 200);
% 
% switch nx
%     case 1
%         set(handles.tAxes(1), 'Position', [dx 110+(h_tot-h) w h]);
%     case 2
%         set(handles.tAxes(1), 'Position', [dx 110+(h_tot-h) w h]);
%         set(handles.tAxes(2), 'Position', [dx 110+(h_tot-2*h-20) w h]);
%     case 3
%         set(handles.tAxes(1), 'Position', [dx 110+(h_tot-h) w h]);
%         set(handles.tAxes(2), 'Position', [dx 110+(h_tot-2*h-20) w h]);
%         set(handles.tAxes(3), 'Position', [dx 110+(h_tot-3*h-40) w h]);
%     case 4        
%         set(handles.tAxes(1), 'Position', [dx 110+(h_tot-h) w h]);
%         set(handles.tAxes(2), 'Position', [dx 110+(h_tot-2*h-20) w h]);
%         set(handles.tAxes(3), 'Position', [dx 110+(h_tot-3*h-40) w h]);
%         set(handles.tAxes(4), 'Position', [dx 110+(h_tot-4*h-60) w h]);
% end
% 
% set(handles.cAxes, 'Position', [dx-100 pos(4)-230 15 200]);
% 
% % frames
% width = pos(3) - 350-50-100-50 -50;
% set(handles.frameLabel, 'Position', [50 pos(4)-20, 100 15]);
% if isfield(handles, 'frameSlider')
%     set(handles.frameSlider, 'Position', [50 75 width 20]);
% end


function handles = setupFrameAxes(hfig, N)

handles = getappdata(hfig, 'handles');
if nargin<2
    N = handles.nCh;
end

pos = get(gcf, 'Position'); % [pixels]

% spacers: 
tspace = 20;
bspace = 100;
lspace = 10;
rspace = 20;%%%%%%%%%%%%% 400;
spacer = 10; % space between panels

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
        if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            width = (width-spacer)/2;
            handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]);
            handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace+width+spacer bspace width height]);
        else % vertical
           height = (height-spacer)/2;
           handles.fPanels(1) = uipanel(uiOpts{:}, 'Position', [lspace bspace+spacer+height width height]);
           handles.fPanels(2) = uipanel(uiOpts{:}, 'Position', [lspace bspace width height]);
        end
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

%===================================
% Plot frames with overlaid tracks
%===================================
function handles = refreshFrameDisplay(hfig)

handles = getappdata(hfig, 'handles');

handles.fAxes

tmp = handles.stack;

% settings = getappdata(hfig, 'settings');
% 
% % save zoom settings
% XLim = get(handles.fAxes(1), 'XLim');
% YLim = get(handles.fAxes(1), 'YLim');
% 
% f = handles.f;
% isRGB = strcmpi(handles.displayType, 'RGB');
% 
% if isRGB
%     if length(handles.fAxes)>1
%         handles = setupFrameAxes(hfig, 1);
%     end
%     cvec = handles.mCh;
%     
% else 
%     if length(handles.fAxes)~=handles.nCh
%         handles = setupFrameAxes(hfig);
%     end
%     cvec = 1:handles.nCh;
% end
% nAxes = length(cvec);
% 
% markerHandles = NaN(1, nAxes);
% 
% for k = 1:nAxes
%     
%     cla(handles.fAxes(k)); % clear axis content
%     
%     % channel index for RGB display
%     if isRGB
%         cidx = 1:min(handles.nCh,3);
%     else
%         cidx = cvec(k);
%     end
%     
%     if ~isempty(handles.tracks{k})
%         chIdx = k;
%     else
%         chIdx = handles.mCh;
%     end
%     
%     if get(handles.('detectionCheckbox'), 'Value') && ~isempty(handles.detection{k})
%         detection = handles.detection{k}(f);
%     else
%         detection = [];
%     end
%     
%     if get(handles.('trackCheckbox'), 'Value') && ~isempty(handles.tracks{cvec(k)})
%         idx = [handles.tracks{cvec(k)}.start]<=f & f<=[handles.tracks{cvec(k)}.end];
%         plotFrame(handles.data, handles.tracks{cvec(k)}(idx), f, cidx,...
%             'Handle', handles.fAxes(cvec(k)), 'iRange', handles.dRange,...
%             'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
%             'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%             'ShowGaps', get(handles.gapCheckbox, 'Value')==1, 'Detection', detection, 'Colormap', handles.colorMap{cvec(k)}(idx,:));
%     else
%         plotFrame(handles.data, [], f, cidx,...
%             'Handle', handles.fAxes(cvec(k)), 'iRange', handles.dRange,...
%             'Mode', handles.displayType, 'Detection', detection);
%     end
%     
%     hold(handles.fAxes(k), 'on');
%     
%     % plot selected track marker
%     if ~isempty(handles.selectedTrack) && get(handles.('trackCheckbox'), 'Value') 
%         selMask = ~isnan(handles.selectedTrack);
%         t = handles.tracks{selMask}(handles.selectedTrack(selMask));
%         fi = f-t.start+1;
%         if 1 <= fi && fi <= length(t.x)
%             xi = t.x(chIdx,fi);
%             yi = t.y(chIdx,fi);
%             markerHandles(k) = plot(handles.fAxes(k), xi, yi, 'ws', 'MarkerSize', 10*settings.zoom);
%         end
%     end
%     
%     if ~isRGB && get(handles.('labelCheckbox'), 'Value')
%         % plot channel name
%         dx = 0.03;
%         text(1-dx*handles.fAspectRatio, dx,...
%             handles.data.markers{k},...
%             'Color', handles.rgbColors{k}, 'Units', 'normalized',...
%             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
%             'Parent', handles.fAxes(k));
%     end
%     
%     % plot EAP status
%     if ~isRGB && get(handles.('eapCheckbox'), 'Value') &&...
%             isfield(handles.tracks{chIdx}, 'significantSignal') && cvec(k) ~= handles.mCh
%         % all tracks
%         tracks = handles.tracks{chIdx};
%         % tracks visible in current frame
%         idx = [tracks.start]<=f & f<=[tracks.end];
%         tracks = tracks(idx);
%         % EAP status
%         eapIdx = [tracks.significantSignal];
%         eapIdx = eapIdx(k,:);
%         % relative position in track
%         fIdx = f-[tracks.start]+1;
%         x = arrayfun(@(i) tracks(i).x(k,fIdx(i)), 1:length(tracks));
%         y = arrayfun(@(i) tracks(i).y(k,fIdx(i)), 1:length(tracks));
%         
%         plot(handles.fAxes(k), x(eapIdx==1), y(eapIdx==1), 'go', 'MarkerSize', 8);
%         plot(handles.fAxes(k), x(eapIdx==0), y(eapIdx==0), 'ro', 'MarkerSize', 8);
%     end
%     
%     hold(handles.fAxes(k), 'off');
% end 
% 
% settings.selectedTrackMarkerID = markerHandles;
% 
% % write zoom level
% set(handles.fAxes(1), 'XLim', XLim);
% set(handles.fAxes(1), 'YLim', YLim);
% 
% setappdata(hfig, 'settings', settings);
% setappdata(hfig, 'handles', handles);
end



% function setColorbar(hfig, mode)
% handles = getappdata(hfig, 'handles');
% 
% lfont = {'FontName', 'Helvetica', 'FontSize', 13};
% sfont = {'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'normal'};
% if ~isempty(handles.tracks{handles.mCh})
%     switch mode
%         case 'Lifetime'
%             df = 40;
%             dcoord = 0.25/df;
%             cmap = [jet(120); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
%             imagesc(reshape(cmap, [size(cmap,1) 1 3]), 'Parent', handles.cAxes);
%             set(handles.cAxes, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
%                'YTick', [1 20:20:120 160],...
%                'YTickLabel', [handles.data.framerate 20:20:120 handles.data.movieLength*handles.data.framerate], sfont{:});
%             text(-0.5, 80, 'Lifetime (s)', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', handles.cAxes, lfont{:});
%         case 'Category'
%             xlabels = {'valid', 'rej. gaps', 'cut', 'persistent',...
%                 'valid', 'rej. gaps', 'cut', 'persistent'};
%             cmap = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
%             imagesc(reshape(cmap, [size(cmap,1) 1 3]), 'Parent', handles.cAxes);
%             set(handles.cAxes, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
%                 'YTick', 1:8, 'YTickLabel', [], 'TickLength', [0 0]);
%             text(-.5, 2.5, 'Single tracks', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', handles.cAxes, lfont{:});
%             text(-.5, 6.5, 'Compound tracks', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', handles.cAxes, lfont{:});
%             
%             for k = 1:8
%                 text(1.6, k-0.1, xlabels{k}, 'Rotation', 45, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Parent', handles.cAxes, sfont{:});
%             end
%         case 'Object Type'
%             cmap = [0 0.8 0; 0.8 0 0];
%             imagesc(reshape(cmap, [size(cmap,1) 1 3]), 'Parent', handles.cAxes);
%             set(handles.cAxes, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
%                 'YTick', 1:8, 'YTickLabel', [], 'TickLength', [0 0]);
%             text(-.5, 1, 'CCP', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', handles.cAxes, lfont{:});
%             text(-.5, 2, 'Other', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', handles.cAxes, lfont{:});
%         otherwise
%             set(handles.cAxes, 'Visible', 'off');
%             cla(handles.cAxes);
%     end
% end



%=========================
% Plot tracks
%=========================
% function refreshTrackDisplay(hfig)
% 
% handles = getappdata(hfig, 'handles');
% 
% if ~isempty(handles.selectedTrack)
%     
%     for ci = 1:handles.nCh
%         h = handles.tAxes(ci);
%         hold(h, 'off');
%         selMask = ~isnan(handles.selectedTrack);
%         sTrack = handles.tracks{selMask}(handles.selectedTrack(selMask));
%                 
%         if size(sTrack.A, 1)==1
%             cx = 1;
%         else
%             cx = ci;
%         end
%         
%         if get(handles.tplotBackgroundCheckbox, 'Value')
%             bgMode = 'zero';
%         else
%             bgMode = 'data';
%         end
%         if strcmpi(handles.pUnitType, 'f')
%             sTrack.t = sTrack.f;
%             if ~isempty(sTrack.startBuffer)
%                 sTrack.startBuffer.t = sTrack.f(1) - (numel(sTrack.startBuffer.t):-1:1);
%                 sTrack.endBuffer.t = sTrack.f(end) + (1:numel(sTrack.startBuffer.t));
%             end
%         end
%         topts = {'Handle', h, 'Time', 'Movie', 'BackgroundValue', bgMode};
%         if get(handles.tplotScaleCheckbox, 'Value')
%             topts = [topts, 'YTick', -handles.yunit(ci):handles.yunit(ci):handles.maxA(ci)];
%         end
%         if ~isempty(handles.bg95) && sTrack.catIdx<5
%             conf = handles.bg95(cx, sTrack.z(cx,:));
%             if ~isempty(sTrack.startBuffer)
%                 conf = [handles.bg95(cx, sTrack.startBuffer.z(cx,:)) conf];
%             end
%             if ~isempty(sTrack.endBuffer)
%                 conf = [conf handles.bg95(cx, sTrack.endBuffer.z(cx,:))];
%             end
%             topts = [topts 'BackgroundConfidence', conf];
%         end
%         plotTrack(handles.data, sTrack, cx, topts{:});
%         box on;
%                      
%         % plot current frame position
%         ybounds = get(h, 'YLim');
%         plot(h, ([handles.f handles.f]-1)*handles.data.framerate, ybounds, '--', 'Color', 0.7*[1 1 1], 'HandleVisibility', 'off');
%         hold(h, 'off');
%         
%         % display result of classification, if available
%         %if isfield(handles.tracks{1}, 'cStatus')
%         %    cStatus = handles.tracks{1}(handles.selectedTrack(1)).cStatus(2);
%         %    if cStatus == 1
%         %        set(handles.statusLabel, 'String', 'Ch. 2: EAP+');
%         %    else
%         %        set(handles.statusLabel, 'String', 'Ch. 2: EAP-');
%         %    end
%         %end
% 
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
%         
%         if isfield(sTrack, 'corrDisappearance') && ci~= handles.mCh
%             
%             if sTrack.corrDisappearance
%                 slabel = 'yes';
%                 scolor = [0 0.8 0];
%             else
%                 slabel = 'no';
%                 scolor = [0.8 0 0];
%             end
%             
%             text(1-dx, 1-4*dx,...
%                 ['Corr. disappearance: ' slabel],...
%                 'Color', scolor, 'Units', 'normalized',...
%                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
%                 'Parent', handles.tAxes(ci));
%             
%         end
%         
%         
%         %     if ~isRGB && get(handles.('labelCheckbox'), 'Value')
%         %         dx = 0.03;
%         %         text(1-dx*handles.fAspectRatio, dx,...
%         %             handles.data.markers{k},...
%         %             'Color', handles.rgbColors{k}, 'Units', 'normalized',...
%         %             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
%         %             'Parent', handles.fAxes(k));
%         %     end
%         
%         
%     end
%     
%     xlabel(h, ['Time (' handles.pUnitType ')']);
% end
% setappdata(hfig, 'handles', handles);


%========================
% Callback functions
%========================

% function refresh_Callback(~,~,hfig)
% refreshFrameDisplay(hfig);
% 
% function refreshTracks_Callback(~,~,hfig)
% refreshTrackDisplay(hfig);
% 
% 
% 
% % function frameSlider_CreateFcn(hObject, ~, ~)
% % if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
% %     set(hObject,'BackgroundColor',[.9 .9 .9]);
% % end
% 
% 
% function frameSlider_Callback(hObject, ~, hfig)
% f = round(get(hObject, 'value'));
% set(hObject, 'Value', f);
% handles = getappdata(hfig, 'handles');
% set(handles.frameLabel, 'String', ['Frame ' num2str(f)]);
% handles.f = f;
% setappdata(hfig, 'handles', handles);
% refreshFrameDisplay(hfig);
% refreshTrackDisplay(hfig);
% 
% 
% 
% function trackButton_Callback(~, ~, hfig)
% 
% handles = getappdata(hfig, 'handles');
% 
% % set focus for next input
% %axes(handles.fAxes(1)); % linked to axes2, selection possible in both
% % set(figure_handle,'CurrentAxes',handles.fAxes(1))
% [x,y] = ginput(1);
% chIdx = find(gca==handles.fAxes);
% 
% if ~isempty(chIdx) && ~isempty(handles.tracks{chIdx})
%     
%     % track segments visible in current frame
%     f = handles.f;
%     idx = find([handles.tracks{chIdx}.start]<=f & f<=[handles.tracks{chIdx}.end]);
%     if ~isempty(idx)
%         np = arrayfun(@(i) numel(i.t), handles.tracks{chIdx}(idx)); % points in each track
%         nt = numel(idx);
%         
%         maxn = max(np);
%         X = NaN(maxn, nt);
%         Y = NaN(maxn, nt);
%         F = NaN(maxn, nt);
%         
%         for k = 1:nt
%             i = 1:np(k);
%             X(i,k) = handles.tracks{chIdx}(idx(k)).x(chIdx,:);
%             Y(i,k) = handles.tracks{chIdx}(idx(k)).y(chIdx,:);
%             F(i,k) = handles.tracks{chIdx}(idx(k)).f;
%         end
%         
%         X(F~=f) = NaN;
%         Y(F~=f) = NaN;
%         mu_x = nanmean(X,1); % average position for compound tracks
%         mu_y = nanmean(Y,1);
%         
%         % nearest point
%         d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
%         handles.selectedTrack = NaN(1,handles.nCh);
%         handles.selectedTrack(chIdx) = idx(d==nanmin(d));
%         set(handles.trackSlider, 'Value', handles.selectedTrack(chIdx));
%         set(handles.trackLabel, 'String', ['Track ' num2str(handles.selectedTrack(chIdx))]);
%         setappdata(hfig, 'handles', handles);
%         % axis(handles.axes3, [0 handles.data.movieLength 0 1]);
%         refreshFrameDisplay(hfig);
%         refreshTrackDisplay(hfig);
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% function statsButton_Callback(~, ~, hfig)
% handles = getappdata(hfig, 'handles');
% if ~isempty(handles.tracks{handles.mCh})
%     
%     % Categories
%     % Ia)  Single tracks with valid gaps
%     % Ib)  Single tracks with invalid gaps
%     % Ic)  Single tracks cut at beginning or end
%     % Id)  Single tracks, persistent
%     % IIa) Compound tracks with valid gaps
%     % IIb) Compound tracks with invalid gaps
%     % IIc) Compound tracks cut at beginning or end
%     % IId) Compound tracks, persistent
%     v = [handles.tracks{handles.mCh}.catIdx];
%     plotTrackClasses(v);
% end


% --- Executes on button press in montageButton.
function montageButton_Callback(~, ~, hfig)
handles = getappdata(hfig, 'handles');

% Creates a montage based on the master track
if ~isempty(handles.selectedTrack)
    fprintf('Generating montage...');
    %options = get(handles.montageOptions, 'String');
    if get(handles.montageAlignCheckbox, 'Value')
        ref = 'Track';
    else
        ref = 'Frame';
    end
    itrack = handles.tracks{handles.mCh}(handles.selectedTrack(1));
    [stack, xa, ya] = getTrackStack(handles.data, itrack,...
        'WindowWidth', 6, 'Reference', ref);
    plotTrackMontage(itrack, stack, xa, ya, 'Labels', handles.data.markers,...
        'ShowMarkers', get(handles.montageMarkerCheckbox, 'Value')==1,...
        'ShowDetection', get(handles.montageDetectionCheckbox, 'Value')==1);
    fprintf(' done.\n');
else
    fprintf('Cannot create montage: no track selected.\n');
end

end

% % --- Executes on selection change in popupmenu1.
% function frameChoice_Callback(hObject, ~, hfig)
% % hObject    handle to popupmenu1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from popupmenu1
% 
% handles = getappdata(hfig, 'handles');
% 
% contents = cellstr(get(hObject,'String'));
% switch contents{get(hObject,'Value')}
%     case 'Raw frames'
%         handles.displayType = 'raw';
%     case 'RGB'
%         handles.displayType = 'RGB';
%     case 'Detection'
%         handles.displayType = 'mask';
% end
% setappdata(hfig, 'handles', handles);
% refreshFrameDisplay(hfig);



function unitChoice_Callback(hObject, ~, hfig)
handles = getappdata(hfig, 'handles');

contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'Seconds'
        handles.pUnitType = 's';
    case 'Frames'
        handles.pUnitType = 'f';
end
setappdata(hfig, 'handles', handles);
refreshTrackDisplay(hfig);

end

% function trackChoice_Callback(hObject, ~, hfig)
% handles = getappdata(hfig, 'handles');
% contents = cellstr(get(hObject, 'String'));
% handles.trackMode = contents{get(hObject,'Value')};
% setColorbar(hfig, handles.trackMode);
% setappdata(hfig, 'handles', handles);
% refreshFrameDisplay(hfig);



% function trackSlider_Callback(hObject, ~, hfig)
% handles = getappdata(hfig, 'handles');
% 
% t = round(get(hObject, 'value'));
% set(hObject, 'Value', t);
% set(handles.trackLabel, 'String', ['Track ' num2str(t)]);
% 
% selMask = ~isnan(handles.selectedTrack);
% handles.selectedTrack(selMask) = t;
% 
% % if track not visible, jump to first frame
% t = handles.tracks{1}(t);
% if handles.f < t.start || handles.f > t.end
%     handles.f = t.start;
%     % set frame number
%     set(handles.frameLabel, 'String', ['Frame ' num2str(handles.f)]);
%     % set frame slider
%     set(handles.frameSlider, 'Value', handles.f);
% end
% 
% setappdata(hfig, 'handles', handles);
% 
% refreshFrameDisplay(hfig);
% refreshTrackDisplay(hfig);






% % --- Executes on button press in printButton.
% function printButton_Callback(~, ~, hfig)
% % hObject    handle to printButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% fprintf('Printing figures ...');
% handles = getappdata(hfig, 'handles');
% 
% % Tracks
% selMask = ~isnan(handles.selectedTrack);
% if ~isempty(selMask) && ~isempty(handles.tracks{selMask})
%     sTrack = handles.tracks{selMask}(handles.selectedTrack(selMask));
%     for ch = 1:handles.nCh
%         plotTrack(handles.data, sTrack, ch,...
%             'FileName', ['track_' num2str(handles.selectedTrack(selMask)) '_ch' num2str(ch) '.eps'],...
%             'Visible', 'off', 'DisplayMode', 'Print');
%     end
%     
%     if get(handles.montageAlignCheckbox, 'Value')
%         ref = 'Track';
%     else
%         ref = 'Frame';
%     end
%     itrack = handles.tracks{handles.mCh}(handles.selectedTrack(1));
%     [stack, xa, ya] = getTrackStack(handles.data, itrack, 'WindowWidth', 6, 'Reference', ref);
%     fpath = [handles.data.source 'Figures' filesep 'track_' num2str(handles.selectedTrack(1)) '_montage.eps'];
%     plotTrackMontage(itrack, stack, xa, ya, 'Labels', handles.data.markers,...
%         'Visible', 'off', 'epsPath', fpath,...
%         'ShowMarkers', get(handles.montageMarkerCheckbox, 'Value')==1,...
%         'ShowDetection', get(handles.montageDetectionCheckbox, 'Value')==1);
% end
% 
% % Frames
% if strcmp(handles.displayType, 'RGB')
%     if ~isempty(handles.tracks{handles.mCh}) && get(handles.('trackCheckbox'), 'Value')
%         idx = [handles.tracks{handles.mCh}.start]<=handles.f & handles.f<=[handles.tracks{handles.mCh}.end];
%     else
%         idx = [];
%     end
%     plotFrame(handles.data, handles.tracks{handles.mCh}(idx), handles.f, 1:min(handles.nCh,3),...
%         'iRange', handles.dRange,...
%         'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
%         'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%         'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
%         'Colormap', handles.colorMap{handles.mCh}(idx,:), 'Print', 'on', 'Visible', 'off');
% else
%     for c = 1:handles.nCh
%         if get(handles.('detectionCheckbox'), 'Value') && ~isempty(handles.detection{k})
%             detection = handles.detection{k}(f);
%         else
%             detection = [];
%         end
%         if ~isempty(handles.tracks{c})
%             idx = [handles.tracks{c}.start]<=handles.f & handles.f<=[handles.tracks{c}.end];
%             plotFrame(handles.data, handles.tracks{c}(idx), handles.f, c,...
%                 'iRange', handles.dRange,...
%                 'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
%                 'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%                 'ShowGaps', get(handles.gapCheckbox, 'Value')==1, 'Detection', detection,...
%                 'Colormap', handles.colorMap{c}(idx,:), 'Print', 'on', 'Visible', 'off');
%         else
%             plotFrame(handles.data, [], handles.f, c,...
%                 'iRange', handles.dRange,...
%                 'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
%                 'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%                 'ShowGaps', get(handles.gapCheckbox, 'Value')==1, 'Detection', detection,...
%                 'Print', 'on', 'Visible', 'off');
%         end
%     end
% end
% 
% fprintf([' done. Figures saved in ' getShortPath(handles.data) filesep 'Figures.\n']);
% 
% 
% 
% function movieButton_Callback(~, ~, hfig)
% 
% handles = getappdata(hfig, 'handles');
% 
% if get(handles.('detectionCheckbox'), 'Value') && ~isempty(handles.detection{k})
%     detection = handles.detection{handles.mCh};
% else
%     detection = [];
% end
% 
% makeMovieCME(handles.data, handles.tracks{handles.mCh}, 'Mode', handles.displayType,...
%     'Detection', detection,...
%     'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%     'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
%     'Displaytype', handles.trackMode, 'Colormap', handles.colorMap{handles.mCh});
% 
% 
% function keyListener(src, evnt)
% 
% handles = getappdata(src, 'handles');
% 
% selMask = ~isnan(handles.selectedTrack);
% itrack = handles.selectedTrack(selMask);
% 
% trackSelect = false;
% switch evnt.Key
%     case 'uparrow'
%         if itrack < numel(handles.tracks{1})
%             itrack = itrack + 1;
%         end
%         trackSelect = true;
%     case 'downarrow'
%         if itrack > 1
%             itrack = itrack - 1;
%         end
%         trackSelect = true;
%     case 'leftarrow'
%         if handles.f>1
%             handles.f = handles.f-1;
%         end
%     case 'rightarrow'
%         if handles.f<handles.data.movieLength
%             handles.f = handles.f+1;
%         end
% end
% 
% if trackSelect
%     handles.selectedTrack(selMask) = itrack;
%     set(handles.trackSlider, 'Value', itrack);
%     set(handles.trackLabel, 'String', ['Track ' num2str(itrack)]);
%     % if track not visible, jump to first frame
%     t = handles.tracks{1}(itrack);
%     if handles.f < t.start || handles.f > t.end
%         handles.f = t.start;
%     end
% end
% 
% % set frame number
% set(handles.frameLabel, 'String', ['Frame ' num2str(handles.f)]);
% % set frame slider
% set(handles.frameSlider, 'Value', handles.f);
% 
% setappdata(src, 'handles', handles);
% 
% refreshFrameDisplay(src);
% refreshTrackDisplay(src);

function ha = setupStackViewer(hf, dims)
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
    ha(1) = axes('Position', [(1-f0)/2+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h], 'Parent', hf);
    ha(2) = axes('Position', [(1-f0)/2 0 f0*f*nz/w f*ny/h], 'Parent', hf);
    ha(3) = axes('Position', [(1-f0)/2+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h], 'Parent', hf);
else
    f0 = h/w * rxy;
    ha(1) = axes('Position', [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h], 'Parent', hf);
    ha(2) = axes('Position', [0 1-f0 f*nz/w f0*f*ny/h], 'Parent', hf);
    ha(3) = axes('Position', [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h], 'Parent', hf);
end

set(hf, 'ResizeFcn', @pResize);

    function pResize(~,~)
        ipos = get(hf, 'Position');
        rxy = ipos(3)/ipos(4);
        dx = spc/ipos(3);
        dy = spc/ipos(4);
        if rxy > w/h % figure is too wide
            f0 = w/h / rxy;
            set(ha(1), 'Position', [(1-f0)/2+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h]);
            set(ha(2), 'Position', [(1-f0)/2 0 f0*f*nz/w f*ny/h]);
            set(ha(3), 'Position', [(1-f0)/2+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h]);
        else
            f0 = h/w * rxy;
            set(ha(1), 'Position', [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h]);
            set(ha(2), 'Position', [0 1-f0 f*nz/w f0*f*ny/h]);
            set(ha(3), 'Position', [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h]);
        end
    end
end
