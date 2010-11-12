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

% Last Modified by GUIDE v2.5 01-Oct-2010 10:04:05

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
masterChannel = regexp(data.source, data.channels);
handles.masterChannel = find([masterChannel{:}]);
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
        handles.dRange{c} = [min([frameInfo.minI]) max([frameInfo.maxI])];
    else
        handles.detection{c} = [];
        % determine dynamic range
        firstFrame = imread(frameList{c}{1});
        lastFrame = imread(frameList{c}{data.movieLength});
        handles.dRange{c} = [min(min(firstFrame(:)),min(lastFrame(:))) max(max(firstFrame(:)),max(lastFrame(:)))];
    end
end
handles.frameList = frameList;
handles.maskList = maskList;
handles.nCh = nChannels;


% initialize handles
handles.f = 1;
handles.displayType = 'raw';
handles.visibleIdx = [];
handles.selectedTrack = [];

h = handles.('slider1');
set(h, 'Min', 1);
set(h, 'Max', data.movieLength);
set(h, 'SliderStep', [1/(data.movieLength-1) 0.05]);

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
        handles.axes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 9*dy]);
    case 2
        if handles.data.imagesize(1) > handles.data.imagesize(2) % horiz.
            handles.axes{1} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 9*dy]);
            handles.axes{2} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 9*dy]);
        else
            handles.axes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 13*dx 4*dy]);
            handles.axes{2} = axes('Parent', gcf, 'Position', [dx 2*dy 13*dx 4*dy]);
        end
    case 3
        handles.axes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.axes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.axes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
    case 4
        handles.axes{1} = axes('Parent', gcf, 'Position', [dx 7*dy 6*dx 4*dy]);
        handles.axes{2} = axes('Parent', gcf, 'Position', [8*dx 7*dy 6*dx 4*dy]);
        handles.axes{3} = axes('Parent', gcf, 'Position', [dx 2*dy 6*dx 4*dy]);
        handles.axes{4} = axes('Parent', gcf, 'Position', [8*dx 2*dy 6*dx 4*dy]);
end

% Update handles structure
guidata(hObject, handles);


% initialize figures/plots
for c = 1:nChannels
    imagesc(imread(frameList{c}{1}), 'Parent', handles.axes{c});%, handles.dRange{c});
end
colormap(gray(256));
linkaxes([handles.axes{:}]);
axis([handles.axes{:}], 'image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename these
axis(handles.axes3, [0 handles.data.movieLength 0 1]);
box(handles.axes3, 'on');


% UIWAIT makes trackDisplayGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = trackDisplayGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

f = round(get(hObject, 'value'));
set(hObject, 'Value', f);
set(handles.('text1'), 'String', ['Frame ' num2str(f)]);
handles.f = f;
guidata(hObject,handles);
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles);



function refreshFrameDisplay(hObject, handles)

cmap = jet(handles.data.movieLength);
f = handles.f;

% save zoom settings
XLim = get(handles.axes{1}, 'XLim');
YLim = get(handles.axes{1}, 'YLim');

for c = 1:handles.nCh
    if ~isempty(handles.detection{c})
        switch handles.displayType
            case 'raw'
                imagesc(imread(handles.frameList{c}{f}), 'Parent', handles.axes{c});%, handles.dRange1);
            case 'WT'
                imagesc(imread(handles.maskList{c}{f}), 'Parent', handles.axes{c});%, handles.dRange1);
            case 'merge'
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
                imagesc(imRGB, 'Parent', handles.axes{c});
        end
    else
        imagesc(imread(handles.frameList{c}{f}), 'Parent', handles.axes{c});%, handles.dRange1);
    end
    axis(handles.axes{c}, 'image');

    hold(handles.axes{c}, 'on');
    if ~isempty(handles.tracks{c})
        idx = find([handles.tracks{c}.start] <= f & f <= [handles.tracks{c}.end]);
        handles.visibleIdx{c} = idx;
        for k = idx
            fi = 1:f-handles.tracks{c}(k).start+1;
            %if f == handles.tracks1(k).start
            plot(handles.axes{c}, handles.tracks{c}(k).x(1), handles.tracks{c}(k).y(1), '*', 'Color', cmap(handles.tracks{c}(k).lifetime,:));
            %end
            plot(handles.axes{c}, handles.tracks{c}(k).x(fi), handles.tracks{c}(k).y(fi), '-', 'Color', cmap(handles.tracks{c}(k).lifetime,:));
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
            plot(handles.axes{c}, t.x(ci), t.y(ci), 'ro', 'MarkerSize', 15);
            text(t.x(ci), t.y(ci), num2str(handles.selectedTrack(c)), 'Color', [1 0 0], 'Parent', handles.axes{c});
        end
    end
    
    % show detection COM values
    if get(handles.('checkbox1'), 'Value')
        plot(handles.axes{c}, handles.detection{c}(f).xcom, handles.detection{c}(f).ycom, 'x', 'Color', hsv2rgb([0/360 0.5 0.5]));
    end
    hold(handles.axes{c}, 'off');
end

% write zoom level
set(handles.axes{c}, 'XLim', XLim);
set(handles.axes{c}, 'YLim', YLim);
guidata(hObject,handles);



function refreshTrackDisplay(handles)

if ~isempty(handles.selectedTrack)

    h = handles.axes3;
    XLim = get(h, 'XLim');
    hold(h, 'off');
        
    for ci = 1:handles.nCh
        if ~isempty(handles.tracks{ci})
            sTrack = handles.tracks{ci}(handles.selectedTrack(1));
        else
            sTrack = handles.tracks{handles.masterChannel}(handles.selectedTrack(1));
        end
        
        bStart = length(sTrack.startBuffer);
        bEnd = length(sTrack.endBuffer);
        t = sTrack.start-bStart:sTrack.end+bEnd;

        %if ~isempty(handles.tracks{ci})
        if size(sTrack.A, 1)==1
            A = [sTrack.startBuffer sTrack.A sTrack.endBuffer];
            c = [sTrack.c(1)*ones(1,bStart) sTrack.c sTrack.c(end)*ones(1,bEnd)];
        else
            A = [sTrack.startBuffer sTrack.A(ci,:) sTrack.endBuffer];
            c = [sTrack.c(ci,1)*ones(1,bStart) sTrack.c(ci,:) sTrack.c(ci,end)*ones(1,bEnd)];
        end
        
        colorA = wavelength2rgb(name2wavelength(handles.data.markers{ci}));
        plot(h, t, A+c, '-', 'Color', colorA);
        hold(h, 'on');
        plot(h, t, c, '--', 'Color', colorA, 'HandleVisibility', 'off');
        %plot(h, t, c+3*cStd, 'k--');
    end
    
    % plot current frame position
    ybounds = get(h, 'YLim');
    plot(h, [handles.f handles.f], ybounds, '--', 'Color', 0.7*[1 1 1], 'HandleVisibility', 'off');
    
    chList = arrayfun(@(x) ['Ch. ' num2str(x) ': ' handles.data.markers{x}], 1:handles.nCh, 'UniformOutput', false);
    legend(h, chList, 'Location', 'SouthEast');
    %legend(h, 'Amplitude ch. 1', 'Background ch. 1', 'Amplitude ch. 2', 'Background ch. 2');

    % retain zoom level
    set(h, 'XLim', XLim);
end




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



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
axes(handles.axes{1}); % linked to axes2, selection possible in both
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
axis(handles.axes3, [0 handles.data.movieLength 0 1]);
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles)



% --- Executes on button press in montageButton.
function montageButton_Callback(~, ~, handles)
% hObject    handle to montageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Creates a montage based on the master track
if ~isempty(handles.selectedTrack)
    
    sigma = getGaussianPSFsigma(handles.data.NA, handles.data.M, handles.data.pixelSize, name2wavelength(handles.data.(['channel' num2str(handles.masterChannel) 'marker'])));
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
        tifFiles = dir([handles.data.(['channel' num2str(c)]) '*.tif*']);
        tifFiles = tifFiles(idx);
        for k = 1:nf
            frame = imread([handles.data.(['channel' num2str(c)]) tifFiles(k).name]);
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
    case 'Detection'
        handles.displayType = 'WT';
    case 'Overlay'
        handles.displayType = 'merge';
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
