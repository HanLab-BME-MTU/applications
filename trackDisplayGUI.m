
% Handles/settings are stored in 'appdata' of the figure handle

function hfig = trackDisplayGUI(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('tracks', cell(1,length(data.channels)), @(x) isstruct(x) || (iscell(x) && numel(x)==numel(data.channels)));
ip.parse(data, varargin{:});
handles.data = data;
tracks = ip.Results.tracks;

% detect number of channels (up to 4)
nCh = length(data.channels);
handles.nCh = nCh;
% exclude master from list of channels
handles.mCh = find(strcmp(data.source, data.channels));
nt = length(tracks);
handles.colorMap = hsv2rgb([rand(nt,1) ones(nt,2)]);

if nCh>4
    error('Only data with up to 4 channels are supported.');
end
    
if isstruct(tracks)
    handles.tracks = cell(1,nCh);
    handles.tracks{1} = tracks;
else
    handles.tracks = tracks;
end

if ~isempty(handles.tracks{handles.mCh})
    handles.maxLifetime_f = max([handles.tracks{handles.mCh}.end]-[handles.tracks{handles.mCh}.start]+1);
else
    handles.maxLifetime_f = [];
end

handles.displayType = 'raw';
if ~all(cellfun(@(x) isempty(x), handles.tracks))
    handles.selectedTrack = ones(1,handles.nCh);
    handles.f = handles.tracks{handles.mCh}(1).start;
else
    handles.selectedTrack = [];
    handles.f = 2; % valid tracks start in frame 2 at the earliest
end


hfig = figure('Units', 'normalized', 'Position', [0.1 0.2 0.8 0.7],...
    'Toolbar', 'figure', 'ResizeFcn', @figResize,...
    'Color', get(0,'defaultUicontrolBackgroundColor'));

set(hfig, 'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels');
pos = get(hfig, 'Position');


%---------------------
% Frames
%---------------------

handles.frameLabel = uicontrol('Style', 'text', 'String', ['Frame ' num2str(handles.f)], ...
    'Position', [20 pos(4)-40, 100 20], 'HorizontalAlignment', 'left');


% Slider
handles.frameSlider = uicontrol('Style', 'slider',...
    'Value', handles.f, 'SliderStep', [1/(data.movieLength-1) 0.05], 'Min', 1, 'Max', data.movieLength,...
    'Position', [20 60 0.6*pos(3) 20], 'Callback', {@frameSlider_Callback, hfig});

uicontrol('Style', 'text', 'String', 'Display: ',...
    'Position', [20 20, 80 20], 'HorizontalAlignment', 'left');

handles.frameChoice = uicontrol('Style', 'popup',...
    'String', {'Raw frames', 'Detection', 'RGB'},...
    'Position', [90 20 120 20], 'Callback', {@frameChoice_Callback, hfig});

% Checkboxes
handles.detectionCheckbox = uicontrol('Style', 'checkbox', 'String', 'Positions',...
    'Position', [250 35 140 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.labelCheckbox = uicontrol('Style', 'checkbox', 'String', 'Channel labels',...
    'Position', [250 10 140 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.trackCheckbox = uicontrol('Style', 'checkbox', 'String', 'Tracks', 'Value', true,...
    'Position', [390 35 100 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.gapCheckbox = uicontrol('Style', 'checkbox', 'String', 'Gaps',...
    'Position', [390 10 60 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});
handles.trackEventCheckbox = uicontrol('Style', 'checkbox', 'String', 'Births/Deaths',...
    'Position', [450 10 100 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});

handles.eapCheckbox = uicontrol('Style', 'checkbox', 'String', 'EAP status',...
    'Position', [560 10 100 20], 'HorizontalAlignment', 'left',...
    'Callback', {@refresh_Callback, hfig});

handles.trackChoice = uicontrol('Style', 'popup',...
    'String', {'Lifetime', 'Category', 'All'},...
    'Position', [460 35 100 20], 'Callback', {@trackChoice_Callback, hfig});

handles.trackButton = uicontrol('Style', 'pushbutton', 'String', 'Select track',...
    'Position', [20+0.6*pos(3)-100 30, 100 28], 'HorizontalAlignment', 'left',...
    'Callback', {@trackButton_Callback, hfig});


%---------------------
% Tracks
%---------------------

handles.trackLabel = uicontrol('Style', 'text', 'String', 'Track 1',...
    'Position', [40+0.6*pos(3) pos(4)-40, 100 20], 'HorizontalAlignment', 'left');

handles.trackSlider = uicontrol('Style', 'slider',...
    'Value', 1, 'SliderStep', [1 1], 'Min', 1, 'Max', 100,...
    'Position', [pos(3)-35 60 20 pos(4)-80],...
    'Callback', {@trackSlider_Callback, hfig});


% Output panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Output', 'Position', [pos(3)-180 5 140 70]);

handles.printButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Print figures',...
    'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.45],...
    'Callback', {@printButton_Callback, hfig});

handles.movieButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Make movie',...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.45],...
    'Callback', {@movieButton_Callback, hfig});
handles.outputPanel = ph;

% Montage panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Montage', 'Position', [pos(3)-390 5 180 70]);
handles.montageButton = uicontrol(ph,'Style','pushbutton','String','Generate',...
    'Units','normalized', 'Position',[.1 .55 .6 .4],...
    'Callback', {@montageButton_Callback, hfig});     
handles.montageText = uicontrol(ph, 'Style', 'text', 'String', 'Align to: ',...
    'Units', 'normalized', 'Position', [0.1 0.1 0.35 0.4], 'HorizontalAlignment', 'left');
handles.montageOptions = uicontrol(ph, 'Style', 'popup',...
    'String', {'Track', 'Frame'},...
    'Units', 'normalized', 'Position', [0.45 0.1 0.5 0.4]);
handles.montageCheckbox = uicontrol('Style', 'checkbox', 'String', 'Show track',...
    'Position', [pos(3)-120 10, 100 20], 'HorizontalAlignment', 'left', 'Visible', 'off');
handles.montagePanel = ph;


setappdata(hfig, 'handles', handles);

%================================


handles.fAspectRatio = handles.data.imagesize(1) / handles.data.imagesize(2);


handles.detection = cell(1,nCh);
handles.dRange = cell(1,nCh);

detectionFile = [data.source 'Detection' filesep 'detection_v2.mat'];
if exist(detectionFile, 'file')==2
    load(detectionFile);
    handles.detection{handles.mCh} = frameInfo;
    if isfield(frameInfo, 'dRange')
        for c = 1:nCh
            M = arrayfun(@(x) x.dRange{c}, frameInfo, 'UniformOutput', false);
            M = vertcat(M{:});
            handles.dRange{c} = [min(M(1,:)) max(M(2,:))];
        end
    end
end

for c = 1:nCh
    if isempty(handles.dRange{c})        
        % determine dynamic range
        firstFrame = double(imread(data.framePaths{c}{1}));
        lastFrame = double(imread(data.framePaths{c}{data.movieLength}));
        handles.dRange{c} = [min(min(firstFrame(:)),min(lastFrame(:))) max(max(firstFrame(:)),max(lastFrame(:)))];
    end
end


% initialize handles
handles.trackMode = 'Lifetime';
handles.hues = getFluorophoreHues(data.markers);
handles.rgbColors = arrayfun(@(x) hsv2rgb([x 1 1]), handles.hues, 'UniformOutput', false);

settings.zoom = 1;
setappdata(hfig, 'settings', settings);


%=================================================
% Set initial values for sliders and checkboxes
%=================================================
if ~isempty([handles.tracks{:}]) && length(handles.tracks{handles.mCh}) > 1
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




%=================================================
% Generate axes
%=================================================
% hFig = findall(0, '-regexp', 'Name', 'trackDisplayGUI')

handles = setupFrameAxes(handles);
dx = 1/23; % unit
dy = 1/12;
switch nCh
    case 1
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 6*dy 7*dx 5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 2
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 7*dy 7*dx 4*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 2*dy 7*dx 4*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 3
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 8.5*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 5.25*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{3} = axes('Parent', gcf, 'Position', [15*dx 2*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 4        
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 9*dy 7*dx 2*dy]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 6.5*dy 7*dx 2*dy]);
        handles.tAxes{3} = axes('Parent', gcf, 'Position', [15*dx 4*dy 7*dx 2*dy]);
        handles.tAxes{4} = axes('Parent', gcf, 'Position', [15*dx 1.5*dy 7*dx 2*dy]);
end
xlabel('Time (s)');


%===========================
% initialize figures/plots
%===========================
for c = 1:nCh
    set(handles.fAxes{c}, 'XLim', [0.5 data.imagesize(2)+0.5], 'YLim', [0.5 data.imagesize(1)+0.5]);
end
linkaxes([handles.tAxes{:}], 'x');
axis([handles.fAxes{:}], 'image');

% save XLim diff. for zoom reference
handles.refXLimDiff = data.imagesize(2)-1;

setappdata(hfig, 'handles', handles);
refreshFrameDisplay(hfig);
refreshTrackDisplay(hfig);

set(zoom, 'ActionPostCallback', {@zoompostcallback, hfig});
% UIWAIT makes trackDisplayGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%===================================
% Automatic actions after zoom
%===================================
function zoompostcallback(~, eventdata, hfig)

XLim = get(eventdata.Axes, 'XLim');

settings = getappdata(hfig, 'settings');
handles = getappdata(hfig, 'handles');

settings.zoom = handles.refXLimDiff / diff(XLim);

for c = 1:length(settings.selectedTrackMarkerID)
    id = settings.selectedTrackMarkerID(c);
    if ~isnan(id)
        set(id, 'MarkerSize', 10*settings.zoom);
    end
end

setappdata(hfig, 'settings', settings);







function figResize(src,~)
pos = get(src, 'Position');
handles = getappdata(src, 'handles');

% frames
set(handles.frameLabel, 'Position', [20 pos(4)-40, 100 20]);
set(handles.frameSlider, 'Position', [20 60 0.6*pos(3) 20]);

% tracks
set(handles.trackLabel, 'Position', [40+0.6*pos(3) pos(4)-40, 100 20]);
set(handles.trackButton, 'Position', [20+0.6*pos(3)-100 30, 100 30]);
set(handles.trackSlider, 'Position', [pos(3)-35 60 20 pos(4)-80]);
set(handles.outputPanel, 'Position', [pos(3)-180 5 140 70]);
set(handles.montagePanel, 'Position', [pos(3)-390 5 180 70]);





function handles = setupFrameAxes(handles, N)

if nargin<2
    N = handles.nCh;
end

dx = 1/23; % unit
dy = 1/12;

if isfield(handles, 'fAxes') && ~isempty(handles.fAxes)
    cellfun(@(x) delete(x), handles.fAxes);
end
handles.fAxes = cell(1,N);
switch N
    case 1
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 9*dy]);
    case 2
        if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 9*dy]);
            handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 9*dy]);
        else
            handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 13*dx 4*dy]);
            handles.fAxes{2} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 4*dy]);
        end
    case 3
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.fAxes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
    case 4
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.fAxes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
        handles.fAxes{4} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 4*dy]);
end
linkaxes([handles.fAxes{:}]);



% % --- Outputs from this function are returned to the command line.
% function varargout = trackDisplayGUI_OutputFcn(~, ~, handles)
% % varargout  cell array for returning output args (see VARARGOUT);
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Get default command line output from handles structure
% varargout{1} = handles.output;



%===================================
% Plot frames with overlaid tracks
%===================================
function handles = refreshFrameDisplay(hfig)

handles = getappdata(hfig, 'handles');
settings = getappdata(hfig, 'settings');

% save zoom settings
XLim = get(handles.fAxes{1}, 'XLim');
YLim = get(handles.fAxes{1}, 'YLim');

% zoomFactor = handles.refXLimDiff / diff(XLim);

f = handles.f;

mc = handles.mCh;

isRGB = strcmpi(handles.displayType, 'RGB');

if isRGB
    if length(handles.fAxes)>1
        handles = setupFrameAxes(handles, 1);
    end
    cvec = mc;
    
else 
    if length(handles.fAxes)~=handles.nCh
        handles = setupFrameAxes(handles);
    end
    cvec = 1:handles.nCh;
end
nAxes = length(cvec);

markerHandles = NaN(1, nAxes);
textHandles = NaN(1, nAxes);


for k = 1:nAxes
    
    cla(handles.fAxes{k}); % clear axis content
    
    % channel index for RGB display
    if isRGB
        cidx = 1:min(handles.nCh,3);
    else
        cidx = cvec(k);
    end
    
    if ~isempty(handles.tracks{k})
        chIdx = k;
    else
        chIdx = mc;
    end
    
    if get(handles.('trackCheckbox'), 'Value') && ~isempty(handles.tracks{cvec(k)})
        plotFrame(handles.data, handles.tracks{cvec(k)}, f, cidx,...
            'Handle', handles.fAxes{cvec(k)}, 'iRange', handles.dRange,...
            'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
            'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
            'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
            'Colormap', handles.colorMap);
    else
        plotFrame(handles.data, [], f, cidx,...
            'Handle', handles.fAxes{cvec(k)}, 'iRange', handles.dRange,...
            'Mode', handles.displayType);
    end
    
    if k==1 && get(handles.('trackCheckbox'), 'Value') 
        setColorbar(hfig, handles.trackMode);
    end
    
    hold(handles.fAxes{k}, 'on');
    
    % plot selected track marker
    if ~isempty(handles.selectedTrack) && get(handles.('trackCheckbox'), 'Value') 
        t = handles.tracks{chIdx}(handles.selectedTrack(k));
        ci = f-t.start+1;
        if 1 <= ci && ci <= length(t.x)
            markerHandles(k) = plot(handles.fAxes{k}, t.x(chIdx,ci), t.y(chIdx,ci), 'ws', 'MarkerSize', 10*settings.zoom);
            textHandles(k) = text(t.x(chIdx,ci)+15, t.y(chIdx,ci)+10, num2str(handles.selectedTrack(k)), 'Color', 'w', 'Parent', handles.fAxes{k});            
        end
    end
    
    % show detection localization values
    if get(handles.('detectionCheckbox'), 'Value') && ~isempty(handles.detection{k})
        d = handles.detection{k}(f);
        if ~isempty(d.x)
            plot(handles.fAxes{k}, d.x(k,:), d.y(k,:), 'ro', 'MarkerSize', 8);
        end
    end
    
    if ~isRGB && get(handles.('labelCheckbox'), 'Value')
        % plot channel name
        %[getDirFromPath(handles.data.channels{k}) '-' handles.data.markers{k}],...
        dx = 0.03;
        text(1-dx*handles.fAspectRatio, dx,...
            handles.data.markers{k},...
            'Color', handles.rgbColors{k}, 'Units', 'normalized',...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
            'Parent', handles.fAxes{k});
    end
    
    % plot EAP status
    if ~isRGB && get(handles.('eapCheckbox'), 'Value') &&...
            isfield(handles.tracks{chIdx}, 'significantSignal') && cvec(k) ~= handles.mCh
        % all tracks
        tracks = handles.tracks{chIdx};
        % tracks visible in current frame
        idx = [tracks.start]<=f & f<=[tracks.end];
        tracks = tracks(idx);
        % EAP status
        eapIdx = [tracks.significantSignal];
        eapIdx = eapIdx(k,:);
        % relative position in track
        fIdx = f-[tracks.start]+1;
        x = arrayfun(@(i) tracks(i).x(k,fIdx(i)), 1:length(tracks));
        y = arrayfun(@(i) tracks(i).y(k,fIdx(i)), 1:length(tracks));
        
        plot(handles.fAxes{k}, x(eapIdx==1), y(eapIdx==1), 'go', 'MarkerSize', 8);
        plot(handles.fAxes{k}, x(eapIdx==0), y(eapIdx==0), 'ro', 'MarkerSize', 8);
    end
    
    hold(handles.fAxes{k}, 'off');
end 

settings.selectedTrackMarkerID = markerHandles;
settings.selectedTrackLabelID = textHandles;

% write zoom level
set(handles.fAxes{1}, 'XLim', XLim);
set(handles.fAxes{1}, 'YLim', YLim);

setappdata(hfig, 'settings', settings);
setappdata(hfig, 'handles', handles);




function setColorbar(hfig, mode)
handles = getappdata(hfig, 'handles');
colorbar('peer', handles.fAxes{1}, 'delete');

if ~isempty(handles.tracks{handles.mCh}) && any(strcmpi(mode, {'lifetime', 'category'}))
    hc = colorbar('peer', handles.fAxes{1}, 'Units', 'normalized', 'Location', 'South');
    
    % Position colorbar on top of figures
    cpos = get(hc, 'Position');
    fpos = get(hfig, 'Position');
    % Top
    cpos(1) = cpos(1)+cpos(3)-150/fpos(3);
    cpos(3) = 150/fpos(3);
    cpos(2) = 1-40/fpos(4);
    cpos(4) = 0.8*cpos(4);
    
    % Right
    %pos(1) = pos(1)+15;
    %pos(2) = pos(2)+pos(4)-100;
    %pos(3) = 0.66*pos(3);
    %pos(4) = 100;
    %ml = handles.data.movieLength*handles.data.framerate;
    
    switch mode
        case 'Lifetime'
            if handles.maxLifetime_f>120
                df = handles.maxLifetime_f-120;
                dc = 0.25/df;
                cmap = [jet(120); (0.5:-dc:0.25+dc)' zeros(df,2)];
            else
                cmap = jet(handles.maxLifetime_f);
            end
            colormap(handles.fAxes{1}, cmap);
            caxis(handles.fAxes{1}, [0 handles.maxLifetime_f])
            set(hc, 'Position', cpos);
            XTick = get(hc, 'XTick');
            if XTick(1)==0
                XTick(1) = handles.data.framerate;
            end
            if XTick(end) < 0.85*handles.maxLifetime_f
                XTick = [XTick handles.maxLifetime_f];
            end
            set(hc, 'XTick', XTick);
        case 'Category'
            colormap(handles.fAxes{1}, [0 1 0; 0 1 1; 1 1 0; 1 0 0]);
            caxis(handles.fAxes{1}, [0 4])
            set(hc, 'Position', cpos);
            set(hc, 'XTick', 0.5:1:3.5, 'TickLength', [0 0], 'XTickLabel', {'Valid', 'M/S', 'Pers.', 'Invalid'});
        %otherwise
            
    end
end



%=========================
% Plot tracks
%=========================
function refreshTrackDisplay(hfig)

handles = getappdata(hfig, 'handles');

if ~isempty(handles.selectedTrack)

    for ci = 1:handles.nCh
        
        h = handles.tAxes{ci};
        hold(h, 'off');

        if ~isempty(handles.tracks{ci})
            sTrack = handles.tracks{ci}(handles.selectedTrack(1));
        else
            sTrack = handles.tracks{handles.mCh}(handles.selectedTrack(1));
        end
        
        if isfield(sTrack, 'startBuffer') && ~isempty(sTrack.startBuffer)
            bStart = size(sTrack.startBuffer.A,2);
        else
            bStart = 0;
        end
        if isfield(sTrack, 'endBuffer') && ~isempty(sTrack.endBuffer)
            bEnd = size(sTrack.endBuffer.A,2);
        else
            bEnd = 0;
        end
        
        if size(sTrack.A, 1)==1
            cx = 1;
        else
            cx = ci;
        end
        
        plotTrack(handles.data, sTrack, handles.selectedTrack, cx, 'Handle', h, 'Legend', 'hide');
        box on;
        l = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
        set(l, 'FontSize', 7);
                     
        % plot current frame position
        ybounds = get(h, 'YLim');
        plot(h, ([handles.f handles.f]-1)*handles.data.framerate, ybounds, '--', 'Color', 0.7*[1 1 1], 'HandleVisibility', 'off');
        axis(handles.tAxes{ci}, [0 handles.data.movieLength ybounds]);
        
        
        % display result of classification, if available
        %if isfield(handles.tracks{1}, 'cStatus')
        %    cStatus = handles.tracks{1}(handles.selectedTrack(1)).cStatus(2);
        %    if cStatus == 1
        %        set(handles.statusLabel, 'String', 'Ch. 2: EAF+');
        %    else
        %        set(handles.statusLabel, 'String', 'Ch. 2: EAF-');
        %    end
        %end
        %pos = get(handles.
        %aspectRatio = 
        dx = 0.03;
        if isfield(handles.tracks{1}, 'significantSignal')
            s = handles.tracks{1}(handles.selectedTrack).significantSignal;
            if s(ci)==1
                slabel = 'yes';
                scolor = [0 0.8 0];
            else
                slabel = 'no';
                scolor = [0.8 0 0];
            end
            text(1-dx, 1-dx,...
                ['Significant: ' slabel],...
                'Color', scolor, 'Units', 'normalized',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
                'Parent', handles.tAxes{ci});
        end
        
        %     if ~isRGB && get(handles.('labelCheckbox'), 'Value')
        %         dx = 0.03;
        %         text(1-dx*handles.fAspectRatio, dx,...
        %             handles.data.markers{k},...
        %             'Color', handles.rgbColors{k}, 'Units', 'normalized',...
        %             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
        %             'Parent', handles.fAxes{k});
        %     end
        
        
    end
    
    % retain zoom level
    set(h, 'XLim', [max(sTrack.start-bStart-11,0) min(sTrack.end+bEnd+9,handles.data.movieLength-1)]*handles.data.framerate);
    xlabel(h, 'Time (s)');
end
setappdata(hfig, 'handles', handles);


%========================
% Callback functions
%========================

function refresh_Callback(~,~,hfig)
refreshFrameDisplay(hfig);



% function frameSlider_CreateFcn(hObject, ~, ~)
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


function frameSlider_Callback(hObject, ~, hfig)
f = round(get(hObject, 'value'));
set(hObject, 'Value', f);
handles = getappdata(hfig, 'handles');
set(handles.frameLabel, 'String', ['Frame ' num2str(f)]);
handles.f = f;
setappdata(hfig, 'handles', handles);
refreshFrameDisplay(hfig);
refreshTrackDisplay(hfig);



function trackButton_Callback(~, ~, hfig)

handles = getappdata(hfig, 'handles');

% set focus for next input
axes(handles.fAxes{1}); % linked to axes2, selection possible in both
[x,y] = ginput(1);

for c = 1:handles.nCh
    % mean position of visible tracks
    if ~isempty(handles.tracks{c})
        chIdx = c;
    else
        chIdx = handles.mCh;
    end   
    idx = find([handles.tracks{chIdx}.start] <= handles.f & handles.f <= [handles.tracks{chIdx}.end]);

    np = length(idx);
    mu_x = zeros(1,length(np));
    mu_y = zeros(1,length(np));
    for k = 1:np
        fi = 1:handles.f-handles.tracks{chIdx}(idx(k)).start+1;
        mu_x(k) = mean(handles.tracks{chIdx}(idx(k)).x(fi));
        mu_y(k) = mean(handles.tracks{chIdx}(idx(k)).y(fi));
    end
    % nearest point
    d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
    selectedIdx = idx(d==min(d));
    handles.selectedTrack(c) = selectedIdx;
end

set(handles.trackSlider, 'Value', handles.selectedTrack(1));
set(handles.trackLabel, 'String', ['Track ' num2str(handles.selectedTrack(1))]);


% % mean position of visible tracks
% idx = handles.visibleIdx{2};
% np = length(idx);
% mu_x = zeros(1,length(np));
% mu_y = zeros(1,length(np));
% for k = 1:np
%     fi = 1:handles.f-handles.tracks2(idx(k)).start+1;
%     mu_x(k) = mean(handles.tracks2(idx(k)).x(fi));
%     mu_y(k) = mean(handles.tracks2(idx(k)).y(fi));
% end
% % nearest point
% d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
% handles.selectedTrack(2) = idx(d==min(d));
% % handles.selectedTrack(2) = handles.trackLinks(selectedIdx);


setappdata(hfig, 'handles', handles);
% axis(handles.axes3, [0 handles.data.movieLength 0 1]);
refreshFrameDisplay(hfig);
refreshTrackDisplay(hfig);



% --- Executes on button press in montageButton.
function montageButton_Callback(~, ~, hfig)
handles = getappdata(hfig, 'handles');

% Creates a montage based on the master track
if ~isempty(handles.selectedTrack)
    fprintf('Generating montage...');
    options = get(handles.montageOptions, 'String');
    
    [stack, x, y] = getTrackStack(handles.data, handles.tracks{handles.mCh}(handles.selectedTrack(1)),...
        'WindowWidth', 6, 'Reference', options{get(handles.montageOptions, 'Value')});
    
    if get(handles.montageCheckbox, 'Value')
        plotTrackMontage(stack, 'Labels', handles.data.markers, 'Mode', 'gray', 'TrackCoords', {x,y});
    else
        plotTrackMontage(stack, 'Labels', handles.data.markers, 'Mode', 'gray');
    end
    fprintf(' done.\n');
else
    fprintf('Cannot create montage: no track selected.');
end



% --- Executes on selection change in popupmenu1.
function frameChoice_Callback(hObject, ~, hfig)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles = getappdata(hfig, 'handles');

contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'Raw frames'
        handles.displayType = 'raw';
    case 'RGB'
        handles.displayType = 'RGB';
    case 'Detection'
        handles.displayType = 'mask';
end
setappdata(hfig, 'handles', handles);
refreshFrameDisplay(hfig);



function trackChoice_Callback(hObject, ~, hfig)
handles = getappdata(hfig, 'handles');
contents = cellstr(get(hObject, 'String'));
handles.trackMode = contents{get(hObject,'Value')};
setColorbar(hfig, handles.trackMode);
setappdata(hfig, 'handles', handles);
refreshFrameDisplay(hfig);



function trackSlider_Callback(hObject, ~, hfig)
handles = getappdata(hfig, 'handles');

t = round(get(hObject, 'value'));
set(hObject, 'Value', t);
set(handles.trackLabel, 'String', ['Track ' num2str(t)]);

handles.selectedTrack = t * ones(1,handles.nCh);

guidata(hObject,handles);

% if track not visible, jump to first frame
t = handles.tracks{1}(t);
if handles.f < t.start || handles.f > t.end
    handles.f = t.start;
    % set frame number
    set(handles.frameLabel, 'String', ['Frame ' num2str(handles.f)]);
    % set frame slider
    set(handles.frameSlider, 'Value', handles.f);
end

setappdata(hfig, 'handles', handles);

refreshFrameDisplay(hfig);
refreshTrackDisplay(hfig);






% --- Executes on button press in printButton.
function printButton_Callback(~, ~, hfig)
% hObject    handle to printButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('Printing...');
handles = getappdata(hfig, 'handles');

for ch = 1:handles.nCh
    if ~isempty(handles.tracks{ch})
        tracks = handles.tracks{ch};
    else
        tracks = handles.tracks{handles.mCh};
    end
    plotTrack(handles.data, tracks, handles.selectedTrack(ch), ch,...
        'Print', 'on', 'Visible', 'off', 'Legend', 'hide');
end




% if get(handles.('trackCheckbox'), 'Value') && ~isempty(handles.tracks{cvec(k)})
%         plotFrame(handles.data, handles.tracks{cvec(k)}, f, cidx,...
%             'Handle', handles.fAxes{cvec(k)}, 'iRange', handles.dRange,...
%             'Mode', handles.displayType, 'DisplayType', handles.trackMode,...
%             'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
%             'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
%             'Colormap', handles.colorMap);
%     else
%         plotFrame(handles.data, [], f, cidx,...
%             'Handle', handles.fAxes{cvec(k)}, 'iRange', handles.dRange,...
%             'Mode', handles.displayType);
%     end


if strcmp(handles.displayType, 'RGB')
    plotFrame(handles.data, handles.tracks{handles.mCh}, handles.f, 1:min(handles.nCh,3),...
        'iRange', handles.dRange, 'Mode', handles.displayType,...
        'ShowDetection', get(handles.('detectionCheckbox'), 'Value')==1,...
        'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
        'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
        'Print', 'on', 'Visible', 'off');
else
    for c = 1:handles.nCh
        plotFrame(handles.data, handles.tracks{c}, handles.f, c,...
            'iRange', handles.dRange, 'Mode', handles.displayType,...
            'ShowDetection', get(handles.('detectionCheckbox'), 'Value')==1,...
            'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
            'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
            'Print', 'on', 'Visible', 'off');     
    end
end


h = handles.montageOptions;
options = get(h, 'String');
stack = getTrackStack(handles.data, handles.tracks{handles.mCh}(handles.selectedTrack(1)), 'WindowWidth', 5, 'Reference', options{get(h, 'Value')});
fpath = [handles.data.source 'Figures' filesep 'track_' num2str(handles.selectedTrack(1)) '_montage.eps'];
plotTrackMontage(stack, 'Labels', handles.data.markers, 'Visible', 'off', 'epsPath', fpath, 'Mode', 'gray');

fprintf(' done.\n');



function movieButton_Callback(~, ~, hfig)

handles = getappdata(hfig, 'handles');
makeMovieCME(handles.data, handles.tracks{handles.mCh}, 'Mode', handles.displayType,...
    'ShowDetection', get(handles.('detectionCheckbox'), 'Value')==1,...
    'ShowEvents', get(handles.trackEventCheckbox, 'Value')==1,...
    'ShowGaps', get(handles.gapCheckbox, 'Value')==1,...
    'Displaytype', handles.trackMode);
