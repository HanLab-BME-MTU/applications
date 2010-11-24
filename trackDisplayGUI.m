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

% Last Modified by GUIDE v2.5 12-Nov-2010 11:08:33

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



% detect number of channels (up to 4)
nChannels = length(data.channels);
% exclude master from list of channels
handles.masterChannel = find(strcmp(data.source, data.channels));
handles.slaveChannels = setdiff(1:nChannels, handles.masterChannel);


for k = 1:length(varargin)-1
    handles.tracks{k} = varargin{k+1};
end
for k = length(varargin):nChannels
    handles.tracks{k} = [];
end
% if numel(varargin)>3
%     handles.trackLinks = varargin{4};
% end



frameList = cell(1,nChannels);
maskList = cell(1,nChannels);
for c = 1:nChannels
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
handles.nCh = nChannels;


% initialize handles
handles.f = 1;
handles.displayType = 'raw';
handles.visibleIdx = [];
handles.selectedTrack = [];


% Set slider values
h = handles.frameSlider;
set(h, 'Min', 1);
set(h, 'Max', data.movieLength);
set(h, 'SliderStep', [1/(data.movieLength-1) 0.05]);

h = handles.trackSlider;
if ~isempty([handles.tracks{:}])
    set(h, 'Min', 1);
    nTracks = length(handles.tracks{handles.masterChannel});
    set(h, 'Max', nTracks);
    set(h, 'SliderStep', [1/(nTracks-1) 0.05]);
else
    set(h, 'Visible', 'off');
end
% Choose default command line output for trackDisplayGUI
handles.output = hObject;


%=====================
% Generate axes
%=====================
% hFig = findall(0, '-regexp', 'Name', 'trackDisplayGUI')
% uicontrol(hFig(1), 'Style', 'slider');

dx = 1/23; % unit
dy = 1/12;
switch nChannels
    case 1
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 9*dy]);
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 6*dy 7*dx 5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 2
        if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 9*dy]);
            handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 9*dy]);
        else
            handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 13*dx 4*dy]);
            handles.fAxes{2} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 4*dy]);
        end
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 7*dy 7*dx 4*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 2*dy 7*dx 4*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 3
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.fAxes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
        
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 8.5*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 5.25*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
        handles.tAxes{3} = axes('Parent', gcf, 'Position', [15*dx 2*dy 7*dx 2.5*dy], 'Box', 'on', 'XLim', [0 handles.data.movieLength]);
    case 4
        handles.fAxes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.fAxes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.fAxes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
        handles.fAxes{4} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 4*dy]);
        
        handles.tAxes{1} = axes('Parent', gcf, 'Position', [15*dx 9*dy 7*dx 2*dy]);
        handles.tAxes{2} = axes('Parent', gcf, 'Position', [15*dx 6.5*dy 7*dx 2*dy]);
        handles.tAxes{3} = axes('Parent', gcf, 'Position', [15*dx 4*dy 7*dx 2*dy]);
        handles.tAxes{4} = axes('Parent', gcf, 'Position', [15*dx 1.5*dy 7*dx 2*dy]);
end

% Update handles structure
guidata(hObject, handles);


% initialize figures/plots
for c = 1:nChannels
    set(handles.fAxes{c}, 'XLim', [1 data.imagesize(2)], 'YLim', [1 data.imagesize(1)]);
end
colormap(gray(256));
linkaxes([handles.fAxes{:}]);
linkaxes([handles.tAxes{:}], 'x');
axis([handles.fAxes{:}], 'image');

refreshFrameDisplay(hObject, handles);
set(zoom, 'ActionPostCallback', @zoompostcallback);
set(gcf, 'UserData', handles.fAxes);

% UIWAIT makes trackDisplayGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);




function zoompostcallback(~,evd)
h = get(get(evd.Axes, 'Parent'), 'UserData');
na = length(h);
if na>2
    XLim = get(evd.Axes, 'XLim');
    YLim = get(evd.Axes, 'YLim');
    dx = min(0.02*diff(XLim), 0.02*diff(YLim)); 
    for k = 1:na
        c = get(h{k}, 'Children');
        set(c(1), 'Position', [XLim(2)-dx, YLim(2)-dx 0]);
    end
end




% --- Outputs from this function are returned to the command line.
function varargout = trackDisplayGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function refreshFrameDisplay(hObject, handles)

cmap = jet(handles.data.movieLength);
f = handles.f;

% save zoom settings
XLim = get(handles.fAxes{1}, 'XLim');
YLim = get(handles.fAxes{1}, 'YLim');

if strcmp(handles.displayType, 'RGB')
    [~, idx] = sort(handles.data.markers);
    switch length(idx)
        case 1
            G = double(imread(handles.frameList{1}{f}));
            imagesc(ch2RGB([], G, []), 'Parent', handles.fAxes{1});
        case 2
            G = double(imread(handles.frameList{idx(1)}{f}));
            R = double(imread(handles.frameList{idx(2)}{f}));
            imagesc(ch2RGB(R, G, []), 'Parent', handles.fAxes{1});
        case 3
            B = double(imread(handles.frameList{idx(1)}{f}));
            G = double(imread(handles.frameList{idx(2)}{f}));
            R = double(imread(handles.frameList{idx(3)}{f}));
            imagesc(ch2RGB(R, G, B), 'Parent', handles.fAxes{1});
    end
    axis(handles.fAxes{1}, 'image');
else
    for c = 1:handles.nCh
        if ~isempty(handles.detection{c})
            switch handles.displayType
                case 'raw'
                    imagesc(imread(handles.frameList{c}{f}), 'Parent', handles.fAxes{c});%, handles.dRange1);
                case 'WT'
                    imagesc(imread(handles.maskList{c}{f}), 'Parent', handles.fAxes{c});%, handles.dRange1);
                case 'overlayWT'
                    overlayColor = [1 0 0];
                    mask = double(imread(handles.maskList{c}{f}));
                    frame = double(imread(handles.frameList{c}{f}));
                    
                    maskIdx = mask~=0;
                    chR = frame;
                    chR(maskIdx) = chR(maskIdx)*overlayColor(1);
                    chG = frame;
                    chG(maskIdx) = chG(maskIdx)*overlayColor(2);
                    chB = frame;
                    chB(maskIdx) = chB(maskIdx)*overlayColor(3);
                    imRGB = uint8(cat(3, scaleContrast(chR, handles.dRange{c}), scaleContrast(chG, handles.dRange{c}), scaleContrast(chB, handles.dRange{c})));
                    imagesc(imRGB, 'Parent', handles.fAxes{c});
            end
        else
            imagesc(imread(handles.frameList{c}{f}), 'Parent', handles.fAxes{c});%, handles.dRange1);
        end
        axis(handles.fAxes{c}, 'image');
        
        hold(handles.fAxes{c}, 'on');
        if ~isempty(handles.tracks{c})
            idx = find([handles.tracks{c}.start] <= f & f <= [handles.tracks{c}.end]);
            handles.visibleIdx{c} = idx;
            for k = idx
                fi = 1:f-handles.tracks{c}(k).start+1;
                %if f == handles.tracks1(k).start
                plot(handles.fAxes{c}, handles.tracks{c}(k).x(1), handles.tracks{c}(k).y(1), '*', 'Color', cmap(handles.tracks{c}(k).lifetime,:));
                %end
                plot(handles.fAxes{c}, handles.tracks{c}(k).x(fi), handles.tracks{c}(k).y(fi), '-', 'Color', cmap(handles.tracks{c}(k).lifetime,:));
                %if f == handles.tracks1(k).end
                %plot(handles.axes1, handles.tracks1(k).x(end), handles.tracks1(k).y(end), '+', 'Color', cmap(handles.tracks1(k).lifetime,:), 'MarkerSize', 8);
                %end
            end
        end
        
        % plot selected track marker
        if ~isempty(handles.selectedTrack)
            if ~isempty(handles.tracks{c})
                chIdx = c;
            else
                chIdx = handles.masterChannel;
            end
            
            t = handles.tracks{chIdx}(handles.selectedTrack(c));
            ci = f-t.start+1;
            if 1 <= ci && ci <= t.lifetime
                plot(handles.fAxes{c}, t.x(ci), t.y(ci), 'ro', 'MarkerSize', 15);
                text(t.x(ci), t.y(ci), num2str(handles.selectedTrack(c)), 'Color', [1 0 0], 'Parent', handles.fAxes{c});
            end
        end
        
        % show detection COM values
        if get(handles.('checkbox1'), 'Value')
            plot(handles.fAxes{c}, handles.detection{c}(f).xcom, handles.detection{c}(f).ycom, 'x', 'Color', hsv2rgb([0/360 0.5 0.5]));
        end
        hold(handles.fAxes{c}, 'off');
        
        if handles.nCh>2
        % plot channel name
        dx = min(0.02*diff(XLim), 0.02*diff(YLim));
        text(XLim(2)-dx, YLim(2)-dx, getDirFromPath(handles.data.channels{c}),...
            'Color', wavelength2rgb(name2wavelength(handles.data.markers{c})),...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Parent', handles.fAxes{c})
        end
    end
end

% write zoom level
set(handles.fAxes{1}, 'XLim', XLim);
set(handles.fAxes{1}, 'YLim', YLim);
guidata(hObject,handles);



function refreshTrackDisplay(handles)

if ~isempty(handles.selectedTrack)

    % Significance thresholds
    % sigmaT = icdf('normal', 1-alpha/2, 0, 1);
    sigmaL = icdf('normal', 0.95, 0, 1); % weaker, single-tailed
    sigmaH = icdf('normal', 0.99, 0, 1);
    
    lh = zeros(1,3);
%     h = handles.axes3;
%     XLim = get(h, 'XLim');
%     hold(h, 'off');
        
    for ci = 1:handles.nCh
        
        h = handles.tAxes{ci};
        %XLim = get(h, 'XLim');
        hold(h, 'off');
        
        if ~isempty(handles.tracks{ci})
            sTrack = handles.tracks{ci}(handles.selectedTrack(1));
        else
            sTrack = handles.tracks{handles.masterChannel}(handles.selectedTrack(1));
        end
        
        bStart = size(sTrack.startBuffer.A,2);
        bEnd = size(sTrack.endBuffer.A,2);

        if size(sTrack.A, 1)==1
            cx = 1;
        else
            cx = ci;
        end
        
        % colors
        trackColor = wavelength2rgb(name2wavelength(handles.data.markers{ci}));
        alpha5c = rgb2hsv(trackColor);
        alpha5c(2) = alpha5c(2)/4;
        alpha1c = alpha5c;
        alpha1c(2) = alpha1c(2)/4;
        
        alpha5cB = alpha5c;
        alpha5cB(3) = 0.9;
        alpha1cB = alpha1c;
        alpha1cB(3) = 0.9;%alpha1cB(
        
        alpha5c = hsv2rgb(alpha5c);
        alpha1c = hsv2rgb(alpha1c);
        alpha5cB = hsv2rgb(alpha5cB);
        alpha1cB = hsv2rgb(alpha1cB);
        
        
        % Plot track
        A = sTrack.A(cx,:);
        c = sTrack.c(cx,:);
        cStd = sTrack.cStd(cx,:);
        t = sTrack.start:sTrack.end;
        
        % alpha = 0.05 level
        lh(3) = fill([t t(end:-1:1)], [c c(end:-1:1)+sigmaL*cStd(end:-1:1)], alpha5c, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'on');
        hold(h, 'on');
        
        % alpha = 0.01 level
        fill([t t(end:-1:1)], [c+sigmaL*cStd c(end:-1:1)+sigmaH*cStd(end:-1:1)], alpha1c, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'off');

        lh(1) = plot(h, t, A+c, '.-', 'Color', trackColor, 'LineWidth', 1);
        lh(2) = plot(h, t, c, '-', 'Color', trackColor, 'HandleVisibility', 'on');

        % Plot left buffer
        if isfield(sTrack, 'startBuffer')
            A = [sTrack.startBuffer.A(cx,:) sTrack.A(cx,1)];
            c = [sTrack.startBuffer.c(cx,:) sTrack.c(cx,1)];
            cStd = [sTrack.startBuffer.cStd(cx,:) sTrack.cStd(cx,1)];
            t = sTrack.start-bStart:sTrack.start;
            
            fill([t t(end:-1:1)], [c c(end:-1:1)+sigmaL*cStd(end:-1:1)], alpha5cB, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'off');
            fill([t t(end:-1:1)], [c+sigmaL*cStd c(end:-1:1)+sigmaH*cStd(end:-1:1)], alpha1cB, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'off');
            plot(h, t, A+c, '.--', 'Color', trackColor, 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(h, t, c, '--', 'Color', trackColor, 'HandleVisibility', 'off');
        end
        
        % Plot right buffer
        if isfield(sTrack, 'endBuffer')
            A = [sTrack.A(cx,end) sTrack.endBuffer.A(cx,:)];
            c = [sTrack.c(cx,end) sTrack.endBuffer.c(cx,:)];
            cStd = [sTrack.cStd(cx,end) sTrack.endBuffer.cStd(cx,:)];
            t = sTrack.end:sTrack.end+bEnd;
            
            fill([t t(end:-1:1)], [c c(end:-1:1)+sigmaL*cStd(end:-1:1)], alpha5cB, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'off');
            fill([t t(end:-1:1)], [c+sigmaL*cStd c(end:-1:1)+sigmaH*cStd(end:-1:1)], alpha1cB, 'EdgeColor', 'none', 'Parent', h, 'HandleVisibility', 'off');
            plot(h, t, A+c, '.--', 'Color', trackColor, 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(h, t, c, '--', 'Color', trackColor, 'HandleVisibility', 'off');
        end
        
        ybounds = get(h, 'YLim');

        % plot current frame position
        plot(h, [handles.f handles.f], ybounds, '--', 'Color', 0.7*[1 1 1], 'HandleVisibility', 'off');

        axis(handles.tAxes{ci}, [0 handles.data.movieLength ybounds]);
        l = legend(lh, ['Amplitude ch. ' num2str(ci)], ['Background ch. ' num2str(ci)], '\alpha = 0.95 level');
        set(l, 'FontSize', 7);
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
    %set(h, 'XLim', XLim);
    set(h, 'XLim', [sTrack.start-bStart-10 sTrack.end+bEnd+10]);
end




% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, eventdata, handles)
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
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, ~, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

refreshFrameDisplay(hObject, handles);



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
    idx = handles.visibleIdx{chIdx};
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
    
    sigma = getGaussianPSFsigma(handles.data.NA, handles.data.M, handles.data.pixelSize, name2wavelength(handles.data.markers{handles.masterChannel}));
    w = ceil(4*sigma);
    
    t = handles.tracks{handles.masterChannel}(handles.selectedTrack(1));
    % buffer with 5 frames before and after
    buffer = 5;
    bStart = t.start - max(1, t.start-buffer);
    bEnd = min(handles.data.movieLength, t.end+buffer) - t.end;
    idx = t.start-bStart:t.end+bEnd;
    nf = length(idx);
    window = cell(length(handles.selectedTrack),nf);
    
    xi = round(t.x);
    yi = round(t.y);
    xi = [xi(1)*ones(1,bStart) xi xi(end)*ones(1,bEnd)];
    yi = [yi(1)*ones(1,bStart) yi yi(end)*ones(1,bEnd)];
    % load all visible frames of this track and store
    for c = [handles.masterChannel handles.slaveChannels]
        tifFiles = dir([handles.data.channels{c} '*.tif*']);
        tifFiles = tifFiles(idx);
        for k = 1:nf
            frame = imread([handles.data.channels{c} tifFiles(k).name]);
            window{c,k} = frame(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
        end
    end
    montagePlot(window);
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
    case 'RGB merge'
        handles.displayType = 'RGB';
    case 'Detection'
        handles.displayType = 'WT';
    case 'Detection overlay'
        handles.displayType = 'overlayWT';
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
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles)




% --- Executes during object creation, after setting all properties.
function trackSlider_CreateFcn(hObject, ~, ~)
% hObject    handle to trackSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
