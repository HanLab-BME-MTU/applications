function varargout = trackDisplay_GUI(varargin)
% TRACKDISPLAY_GUI M-file for trackDisplay_GUI.fig
%      TRACKDISPLAY_GUI, by itself, creates a new TRACKDISPLAY_GUI or raises the existing
%      singleton*.
%
%      H = TRACKDISPLAY_GUI returns the handle to a new TRACKDISPLAY_GUI or the handle to
%      the existing singleton*.
%
%      TRACKDISPLAY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKDISPLAY_GUI.M with the given input arguments.
%
%      TRACKDISPLAY_GUI('Property','Value',...) creates a new TRACKDISPLAY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackDisplay_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackDisplay_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackDisplay_GUI

% Last Modified by GUIDE v2.5 24-Sep-2010 00:46:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackDisplay_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackDisplay_GUI_OutputFcn, ...
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


% --- Executes just before trackDisplay_GUI is made visible.
function trackDisplay_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackDisplay_GUI (see VARARGIN)

data = varargin{1};

if isfield(data, 'channel1')
    framesCh1 = dir([data.channel1 '*.tif']);
    frameListCh1 = cellfun(@(x,y) [x y], repmat({data.channel1}, [1 data.movieLength]), {framesCh1.name}, 'UniformOutput', false);
    framesCh2 = dir([data.channel2 '*.tif']);
    frameListCh2 = cellfun(@(x,y) [x y], repmat({data.channel2}, [1 data.movieLength]), {framesCh2.name}, 'UniformOutput', false);
    handles.nCh = 2;
else
    framesCh1 = dir([data.source '*.tif']);
    frameListCh1 = cellfun(@(x,y) [x y], repmat({data.source}, [1 data.movieLength]), {framesCh1.name}, 'UniformOutput', false);

    mpath = [data.source 'Detection' filesep 'Masks' filesep];
    framesCh2 = dir([mpath '*.tif']);
    frameListCh2 = cellfun(@(x,y) [x y], repmat({mpath}, [1 data.movieLength]), {framesCh2.name}, 'UniformOutput', false);
    handles.nCh = 1;
end

frameCh1 = imread(frameListCh1{1});
frameCh2 = imread(frameListCh2{1});

% Choose default command line output for trackDisplay_GUI
handles.output = hObject;

% initialize handles
handles.frameListCh1 = frameListCh1;
handles.frameListCh2 = frameListCh2;
handles.f = 1;

handles.data = data;
handles.tracks1 = varargin{2};
if numel(varargin)>2
    handles.tracks2 = varargin{3};
end

if isfield(data, 'channel1')
    load([data.channel1 'Detection' filesep 'detectionResults.mat']);
    handles.detection1 = frameInfo;
    load([data.channel2 'Detection' filesep 'detectionResults.mat']);
    handles.detection2 = frameInfo;
else
    load([data.source 'Detection' filesep 'detectionResults.mat']);
    handles.detection = frameInfo;
end

handles.visibleIdx = [];
handles.selectedTrack = [];



% Update handles structure
guidata(hObject, handles);

% set(handles.frameList,'String',frameList);
% load([data.source 'TrackInfoMatrices' filesep 'trackinfo.mat']);

h = handles.('slider1');
set(h, 'Min', 1);
set(h, 'Max', data.movieLength);
set(h, 'SliderStep', [1/(data.movieLength-1) 0.05]);


imagesc(frameCh1, 'Parent', handles.('axes1'));
imagesc(frameCh2, 'Parent', handles.('axes2'));
colormap(gray(256));
linkaxes([handles.('axes1') handles.('axes2')]);
axis([handles.('axes1') handles.('axes2')],'image');

axis(handles.('axes3'), [0 handles.data.movieLength 0 1]);


% UIWAIT makes trackDisplay_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = trackDisplay_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% guidata(gcbo) % same as 'handles'
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

% save zoom properties
haxes = handles.('axes1');
XLim = get(haxes, 'XLim');
YLim = get(haxes, 'YLim');

frameCh1 = imread(handles.frameListCh1{f});
frameCh2 = imread(handles.frameListCh2{f});

idx1 = find([handles.tracks1.start] <= f & f <= [handles.tracks1.end]);
handles.visibleIdx1 = idx1;
idx2 = find([handles.tracks2.start] <= f & f <= [handles.tracks2.end]);
handles.visibleIdx2 = idx2;


imagesc(frameCh1, 'Parent', handles.('axes1'));
axis(handles.('axes1'),'image');
% overlay tracks only if 2-channel input
if isfield(handles, 'tracks2')
    hold(handles.('axes1'), 'on');
    for k = idx1
        fi = 1:f-handles.tracks1(k).start+1;
        %plot(handles.tracks1(k).x(fi), handles.tracks1(k).y(fi), '-', 'Color', cmap(handles.tracks1(k).lifetime,:));
        plot(handles.axes1, handles.tracks1(k).x(1), handles.tracks1(k).y(1), '*', 'Color', 'r');
        plot(handles.axes1, handles.tracks1(k).x(fi), handles.tracks1(k).y(fi), '-', 'Color', 'r');
    end
%     hold(handles.('axes1'), 'off');
end


imagesc(frameCh2, 'Parent', handles.('axes2'));
axis(handles.('axes2'),'image');
hold(handles.('axes2'), 'on');
% for k = idx
%     fi = 1:f-tracks(k).start+1;
%     plot(handles.('axes2'), handles.tracks1(k).x(1), handles.tracks1(k).y(1), '*', 'Color', cmap(tracks(k).lifetime,:));
%     plot(handles.('axes2'), handles.tracks1(k).x(fi), handles.tracks1(k).y(fi), '-', 'Color', cmap(tracks(k).lifetime,:));
% end
if isfield(handles, 'tracks2')
    tracks2 = handles.('tracks2');
    for k = idx2
        fi = 1:f-tracks2(k).start+1;
        plot(handles.axes2, handles.tracks2(k).x(1), handles.tracks2(k).y(1), '*', 'Color', 'b');
        plot(handles.axes2, handles.tracks2(k).x(fi), handles.tracks2(k).y(fi), '-', 'Color', 'b');
    end
end

% plot selected track marker
if length(handles.selectedTrack)==2
    t = handles.tracks1(handles.selectedTrack(1));
    ci = f-t.start+1;
    if 1 <= ci && ci <= t.lifetime
        plot(handles.('axes1'), t.x(ci), t.y(ci), 'ro', 'MarkerSize', 15);
        text(t.x(ci), t.y(ci), num2str(handles.selectedTrack(1)), 'Color', [1 0 0], 'Parent', handles.('axes1'));
    end
    t = handles.tracks2(handles.selectedTrack(2));
    ci = f-t.start+1;
    if 1 <= ci && ci <= t.lifetime
        plot(handles.('axes2'), t.x(ci), t.y(ci), 'ro', 'MarkerSize', 15);
        text(t.x(ci), t.y(ci), num2str(handles.selectedTrack(2)), 'Color', [1 0 0], 'Parent', handles.('axes2'));
    end
end

% show detection COM values
if get(handles.('checkbox1'), 'Value')
    if handles.nCh == 1
        plot(handles.('axes2'), handles.detection(f).xcom, handles.detection(f).ycom, 'x', 'Color', hsv2rgb([120/360 0.5 0.5]));
    else
        plot(handles.('axes1'), handles.detection1(f).xcom, handles.detection1(f).ycom, 'x', 'Color', hsv2rgb([0/360 0.5 0.5]));
        plot(handles.('axes2'), handles.detection2(f).xcom, handles.detection2(f).ycom, 'x', 'Color', hsv2rgb([0/360 0.5 0.5]));
    end
end
hold(handles.('axes2'), 'off');

% write zoom level
set(haxes, 'XLim', XLim);
set(haxes, 'YLim', YLim);
guidata(hObject,handles);



function refreshTrackDisplay(handles)

if ~isempty(handles.selectedTrack)

    h = handles.('axes3');
    
    sTrack = handles.tracks1(handles.selectedTrack(1));
    hold(h, 'off');
    plot(h, sTrack.t, sTrack.A + sTrack.c, 'r');
    hold(h, 'on');
    plot(h, sTrack.t, sTrack.c, 'k');
    plot(h, sTrack.t, sTrack.c + 3*sTrack.cStd, 'k--');
    
    sTrack = handles.tracks2(handles.selectedTrack(2));
    plot(h, sTrack.t, sTrack.A + sTrack.c, 'b');
    
    
    ybounds = get(h, 'YLim');
    plot(h, [handles.f handles.f], ybounds, '--', 'Color', 0.7*[1 1 1]);
    
    xlim(h, [0 handles.data.movieLength]);
    legend(h, 'Amplitude', 'Background');
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
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

refreshFrameDisplay(hObject, handles);



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set focus for next input
axes(handles.('axes1')); % linked to axes2, selection possible in both
[x,y] = ginput(1);

% mean position of visible tracks
idx = handles.visibleIdx1;
np = length(idx);
mu_x = zeros(1,length(np));
mu_y = zeros(1,length(np));
for k = 1:np
    fi = 1:handles.f-handles.tracks1(idx(k)).start+1;
    mu_x(k) = mean(handles.tracks1(idx(k)).x(fi));
    mu_y(k) = mean(handles.tracks1(idx(k)).y(fi));
end
% nearest point
d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
handles.selectedTrack(1) = idx(d==min(d));

% mean position of visible tracks
idx = handles.visibleIdx2;
np = length(idx);
mu_x = zeros(1,length(np));
mu_y = zeros(1,length(np));
for k = 1:np
    fi = 1:handles.f-handles.tracks1(idx(k)).start+1;
    mu_x(k) = mean(handles.tracks2(idx(k)).x(fi));
    mu_y(k) = mean(handles.tracks2(idx(k)).y(fi));
end
% nearest point
d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
handles.selectedTrack(2) = idx(d==min(d));

guidata(hObject,handles);
refreshFrameDisplay(hObject, handles);
refreshTrackDisplay(handles)



% --- Executes on button press in montageButton.
function montageButton_Callback(hObject, eventdata, handles)
% hObject    handle to montageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.selectedTrack)
    t = handles.tracks1(handles.selectedTrack(1));
    % load all visible frames of this track and store
    
    tifFiles = dir([handles.data.source '*.tif*']);
    
    % buffer with 5 frames before and after
    buffer = 5;
    bStart = t.start - max(1, t.start-buffer);
    bEnd = min(handles.data.movieLength, t.end+buffer) - t.end;
    
    xi = round(t.x);
    yi = round(t.y);
    xi = [xi(1)*ones(1,bStart) xi xi(end)*ones(1,bEnd)];
    yi = [yi(1)*ones(1,bStart) yi yi(end)*ones(1,bEnd)];
    
    tifFiles = tifFiles(t.start-bStart:t.end+bEnd);
    nf = length(tifFiles);
    sigma = 1.628;
    w = ceil(4*sigma);
    window = cell(1,nf);
    for k = 1:nf
        frame = imread([handles.data.source tifFiles(k).name]);
        window{k} = frame(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    end
    montagePlot(window, 12)
end



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end