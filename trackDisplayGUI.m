function varargout = trackDisplayGUI(varargin)
% TRACKDISPLAYGUI M-file for trackDisplayGUI.fig
%      TRACKDISPLAYGUI, by itself, creates a new TRACKDISPLAYGUI or raises the existing
%      singleton*.
%
%      H = TRACKDISPLAYGUI returns the handle to a new TRACKDISPLAYGUI or the handle to
%      the existing singleton*.
%
%      TRACKDISPLAYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKDISPLAYGUI.M with the given input arguments.
%
%      TRACKDISPLAYGUI('Property','Value',...) creates a new TRACKDISPLAYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackDisplayGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackDisplayGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackDisplayGUI

% Last Modified by GUIDE v2.5 17-Mar-2011 10:29:06

% Francois Aguet, September 2010

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackDisplayGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackDisplayGUI_OutputFcn, ...
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


% --- Executes just before trackDisplayGUI is made visible.
function trackDisplayGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackDisplayGUI (see VARARGIN)

data = varargin{1};
handles.data = data;
handles.fAspectRatio = handles.data.imagesize(1) / handles.data.imagesize(2);


% detect number of channels (up to 4)
nCh = length(data.channels);
% exclude master from list of channels
handles.masterChannel = find(strcmp(data.source, data.channels));
handles.slaveChannels = setdiff(1:nCh, handles.masterChannel);


for k = 1:length(varargin)-1
    handles.tracks{k} = varargin{k+1};
end
for k = length(varargin):nCh
    handles.tracks{k} = [];
end
% if numel(varargin)>3
%     handles.trackLinks = varargin{4};
% end



frameList = cell(1,nCh);
maskList = cell(1,nCh);
for c = 1:nCh
    frames = dir([data.channels{c} '*.tif']);
    frameList{c} = cellfun(@(x) [data.channels{c} x], {frames.name}, 'UniformOutput', false);
    maskPath = [data.channels{c} 'Detection' filesep 'Masks' filesep];
    masks = dir([maskPath '*.tif']);
    maskList{c} = cellfun(@(x) [maskPath x], {masks.name}, 'UniformOutput', false);
    
    detectionFile = [data.channels{c} 'Detection' filesep 'detectionResults.mat'];
    if exist(detectionFile, 'file')==2
        load(detectionFile);
        handles.detection{c} = frameInfo;
        %handles.dRange{c} = [min([frameInfo.minI]) max([frameInfo.maxI])];
    else
        handles.detection{c} = [];
    end
    % determine dynamic range
    firstFrame = double(imread(frameList{c}{1}));
    lastFrame = double(imread(frameList{c}{data.movieLength}));
    handles.dRange{c} = [min(min(firstFrame(:)),min(lastFrame(:))) max(max(firstFrame(:)),max(lastFrame(:)))];
end
handles.frameList = frameList;
handles.maskList = maskList;
handles.nCh = nCh;


% initialize handles
handles.f = 2; % valid tracks start in frame 2 at the earliest
set(handles.frameLabel, 'String', 'Frame 2');
handles.displayType = 'raw';
if ~all(cellfun(@(x) isempty(x), handles.tracks))
    handles.selectedTrack = ones(1,handles.nCh);
else
    handles.selectedTrack = [];
end

handles.hues = getHuesFromMarkers(data.markers);
handles.rgbColors = arrayfun(@(x) hsv2rgb([x 1 1]), handles.hues, 'UniformOutput', false);


settings.zoom = 1;
setappdata(handles.figure1, 'mydata', settings);


%=================================================
% Set initial values for sliders and checkboxes
%=================================================
h = handles.frameSlider;
set(h, 'Min', 1);
set(h, 'Max', data.movieLength);
set(h, 'Value', 2);
set(h, 'SliderStep', [1/(data.movieLength-1) 0.05]);

h = handles.trackSlider;
if ~isempty([handles.tracks{:}]) && length(handles.tracks{handles.masterChannel}) > 1
    set(h, 'Min', 1);
    nTracks = length(handles.tracks{handles.masterChannel});
    set(h, 'Max', nTracks);
    set(h, 'SliderStep', [1/(nTracks-1) 0.05]);
else
    set(h, 'Visible', 'off');
end
% Choose default command line output for trackDisplayGUI
handles.output = hObject;


if nCh > 2
    set(handles.('labelCheckbox'), 'Value', 1);
end




%=================================================
% Generate axes
%=================================================
% hFig = findall(0, '-regexp', 'Name', 'trackDisplayGUI')
% uicontrol(hFig(1), 'Style', 'slider');

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
    set(handles.fAxes{c}, 'XLim', [1 data.imagesize(2)], 'YLim', [1 data.imagesize(1)]);
end
colormap(gray(256));
linkaxes([handles.tAxes{:}], 'x');
axis([handles.fAxes{:}], 'image');

% save XLim diff. for zoom reference
handles.refXLimDiff = data.imagesize(2)-1;
handles = refreshFrameDisplay(hObject, handles);


% init. track display
set(handles.trackLabel, 'String', 'Track 1');
refreshTrackDisplay(handles);

guidata(hObject, handles);
% set(zoom, 'ActionPostCallback', {@zoompostcallback, handles, hObject});
set(zoom, 'ActionPostCallback', {@zoompostcallback, handles});
guidata(hObject, handles);
% UIWAIT makes trackDisplayGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



%===================================
% Automatic actions after zoom
%===================================
function zoompostcallback(~, eventdata, handles)

XLim = get(eventdata.Axes, 'XLim');

settings = getappdata(handles.figure1, 'mydata');
settings.zoom = handles.refXLimDiff / diff(XLim);

for c = 1:length(settings.selectedTrackMarkerID)
    id = settings.selectedTrackMarkerID(c);
    if ~isnan(id)
        set(id, 'MarkerSize', 10*settings.zoom);
    end
end

setappdata(handles.figure1, 'mydata', settings);




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



% --- Outputs from this function are returned to the command line.
function varargout = trackDisplayGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%===================================
% Plot frames with overlaid tracks
%===================================
function handles = refreshFrameDisplay(hObject, handles)

% save zoom settings
XLim = get(handles.fAxes{1}, 'XLim');
YLim = get(handles.fAxes{1}, 'YLim');

% zoomFactor = handles.refXLimDiff / diff(XLim);

f = handles.f;
settings = getappdata(handles.figure1, 'mydata');


if strcmp(handles.displayType, 'RGB')
    if length(handles.fAxes)>1
        handles = setupFrameAxes(handles, 1);
    end

    mc = handles.masterChannel;
    if get(handles.('trackCheckbox'), 'Value')
        plotFrame(handles.data, handles.tracks{mc}, f, 1:min(handles.nCh,3),...
            'Handle', handles.fAxes{1}, 'iRange', handles.dRange,...
            'Mode', handles.displayType);
    else
        plotFrame(handles.data, [], f, 1:min(handles.nCh,3),...
            'Handle', handles.fAxes{1}, 'iRange', handles.dRange,...
            'Mode', handles.displayType);
    end
    markerHandles = NaN;
    textHandles = NaN;
    % plot selected track marker
    if ~isempty(handles.selectedTrack) && get(handles.('trackCheckbox'), 'Value')
        hold(handles.fAxes{1}, 'on');
        t = handles.tracks{mc}(handles.selectedTrack(mc));
        ci = f-t.start+1;
        if 1 <= ci && ci <= length(t.x)
            markerHandles = plot(handles.fAxes{mc}, t.x(ci), t.y(ci), 'ws', 'MarkerSize', 10*settings.zoom);
            textHandles = text(t.x(ci)+15, t.y(ci)+10, num2str(handles.selectedTrack(mc)), 'Color', 'w', 'Parent', handles.fAxes{mc});
        end
        hold(handles.fAxes{1}, 'off');
    end
    
else
    if length(handles.fAxes)~=handles.nCh
        handles = setupFrameAxes(handles);
    end

    markerHandles = NaN(1, handles.nCh);
    textHandles = NaN(1, handles.nCh);
    for c = 1:handles.nCh
        if ~isempty(handles.tracks{c})
            chIdx = c;
        else
            chIdx = handles.masterChannel;
        end
        
        if get(handles.('trackCheckbox'), 'Value')
            plotFrame(handles.data, handles.tracks{c}, f, c,...
                'Handle', handles.fAxes{c}, 'iRange', handles.dRange,...
                'Mode', handles.displayType);
        else
            plotFrame(handles.data, [], f, c,...
                'Handle', handles.fAxes{c}, 'iRange', handles.dRange,...
                'Mode', handles.displayType);
        end
        
        hold(handles.fAxes{c}, 'on');
              
        % plot selected track marker
        if ~isempty(handles.selectedTrack) && get(handles.('trackCheckbox'), 'Value')
            t = handles.tracks{chIdx}(handles.selectedTrack(c));
            ci = f-t.start+1;
            if 1 <= ci && ci <= length(t.x)
                markerHandles(c) = plot(handles.fAxes{c}, t.x(ci), t.y(ci), 'ws', 'MarkerSize', 10*settings.zoom);
                textHandles(c) = text(t.x(ci)+15, t.y(ci)+10, num2str(handles.selectedTrack(c)), 'Color', 'w', 'Parent', handles.fAxes{c});
            end
        end
        
        % show detection COM values
        if get(handles.('detectionCheckbox'), 'Value') && ~isempty(handles.detection{c})
            d = handles.detection{c}(f);
            plot(handles.fAxes{c}, d.xcom, d.ycom, 'ro', 'MarkerSize', 8);
            if isfield(d, 'xloc') && ~isempty(d.xloc)
                plot(handles.fAxes{c}, d.xloc, d.yloc, 'gx', 'MarkerSize', 8);
            end
        end
        hold(handles.fAxes{c}, 'off');
        
        if get(handles.('labelCheckbox'), 'Value')
            % plot channel name
            %[getDirFromPath(handles.data.channels{c}) '-' handles.data.markers{c}],...
            dx = 0.03;
            text(1-dx*handles.fAspectRatio, dx,...
                handles.data.markers{c},...
                'Color', handles.rgbColors{c}, 'Units', 'normalized',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                'Parent', handles.fAxes{c});
        end
    end
end

settings.selectedTrackMarkerID = markerHandles;
settings.selectedTrackLabelID = textHandles;
setappdata(handles.figure1, 'mydata', settings);

% write zoom level
set(handles.fAxes{1}, 'XLim', XLim);
set(handles.fAxes{1}, 'YLim', YLim);
guidata(hObject, handles);







%=========================
% Plot tracks
%=========================
function refreshTrackDisplay(handles)

if ~isempty(handles.selectedTrack)

    for ci = 1:handles.nCh
        
        h = handles.tAxes{ci};
        hold(h, 'off');

        if ~isempty(handles.tracks{ci})
            sTrack = handles.tracks{ci}(handles.selectedTrack(1));
        else
            sTrack = handles.tracks{handles.masterChannel}(handles.selectedTrack(1));
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
        
        plotTrack(handles.data, sTrack, handles.selectedTrack, cx, 'Handle', h);
        box on;
        l = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
        set(l, 'FontSize', 7);
                     
        % plot current frame position
        ybounds = get(h, 'YLim');
        plot(h, ([handles.f handles.f]-1)*handles.data.framerate, ybounds, '--', 'Color', 0.7*[1 1 1], 'HandleVisibility', 'off');
        axis(handles.tAxes{ci}, [0 handles.data.movieLength ybounds]);
    end
  
    % display result of classification, if available
    if isfield(handles.tracks{1}, 'cStatus')
        cStatus = handles.tracks{1}(handles.selectedTrack(1)).cStatus(2);
        if cStatus == 1
            set(handles.statusLabel, 'String', 'Ch. 2: EAF+');
        else
            set(handles.statusLabel, 'String', 'Ch. 2: EAF-');
        end
    end
    % retain zoom level
    set(h, 'XLim', [max(sTrack.start-bStart-11,0) min(sTrack.end+bEnd+9,handles.data.movieLength-1)]*handles.data.framerate);
    xlabel(h, 'Time (s)');
end




% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, ~, ~)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function frameSlider_Callback(hObject, ~, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

f = round(get(hObject, 'value'));
set(hObject, 'Value', f);
set(handles.frameLabel, 'String', ['Frame ' num2str(f)]);
handles.f = f;
guidata(hObject,handles);
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles);



% --- Executes when figure1 is resized.
function figure1_ResizeFcn(~, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set focus for next input
axes(handles.fAxes{1}); % linked to axes2, selection possible in both
[x,y] = ginput(1);

for c = 1:handles.nCh
    % mean position of visible tracks
    if ~isempty(handles.tracks{c})
        chIdx = c;
    else
        chIdx = handles.masterChannel;
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

guidata(hObject,handles);
% axis(handles.axes3, [0 handles.data.movieLength 0 1]);
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles)



% --- Executes on button press in montageButton.
function montageButton_Callback(~, ~, handles)
% hObject    handle to montageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Creates a montage based on the master track
if ~isempty(handles.selectedTrack)
    stack = getTrackStack(handles.data, handles.tracks{handles.masterChannel}(handles.selectedTrack(1)));
    montagePlot(stack, 'Labels', handles.data.markers, 'Mode', 'RGB+gray');
else
    fprintf('Cannot create montage: no track selected.');
end



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, ~, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'Raw frames'
        handles.displayType = 'raw';
    case 'RGB'
        handles.displayType = 'RGB';
    case 'Detection'
        handles.displayType = 'mask';
end
guidata(hObject,handles);
refreshFrameDisplay(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function trackSlider_Callback(hObject, ~, handles)
% hObject    handle to trackSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

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

refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles);




% --- Executes during object creation, after setting all properties.
function trackSlider_CreateFcn(hObject, ~, ~)
% hObject    handle to trackSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in printButton.
function printButton_Callback(~, ~, handles)
% hObject    handle to printButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for ch = 1:handles.nCh
    if ~isempty(handles.tracks{ch})
        tracks = handles.tracks{ch};
    else
        tracks = handles.tracks{handles.masterChannel};
    end
    plotTrack(handles.data, tracks, handles.selectedTrack(ch), ch,...
        'Print', 'on', 'Visible', 'off');
end


if strcmp(handles.displayType, 'RGB')
    mc = handles.masterChannel;
    plotFrame(handles.data, handles.tracks{mc}, handles.f, 1:min(handles.nCh,3),...
        'iRange', handles.dRange, 'Mode', handles.displayType,...
            'Print', 'on', 'Visible', 'off');
else
    for c = 1:handles.nCh
        plotFrame(handles.data, handles.tracks{c}, handles.f, c,...
            'iRange', handles.dRange, 'Mode', handles.displayType,...
            'Print', 'on', 'Visible', 'off');
    end
end

stack = getTrackStack(handles.data, handles.tracks{handles.masterChannel}(handles.selectedTrack(1)));
fpath = [handles.data.source 'Figures' filesep 'track_' num2str(handles.f) '_montage.eps'];
montagePlot(stack, 'Labels', handles.data.markers, 'Visible', 'off', 'Print', fpath, 'Mode', 'RGB+gray');

fprintf('Printing done.\n');


% --- Executes on button press in detectionCheckbox.
function detectionCheckbox_Callback(hObject, ~, handles)
% hObject    handle to detectionCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshFrameDisplay(hObject, handles);


% --- Executes on button press in labelCheckbox.
function labelCheckbox_Callback(hObject, ~, handles)
% hObject    handle to labelCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshFrameDisplay(hObject, handles);


% --- Executes on button press in trackCheckbox.
function trackCheckbox_Callback(hObject, ~, handles)
% hObject    handle to trackCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshFrameDisplay(hObject, handles);
