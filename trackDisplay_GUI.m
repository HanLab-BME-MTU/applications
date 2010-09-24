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
else
    framesCh1 = dir([data.source '*.tif']);
    frameListCh1 = cellfun(@(x,y) [x y], repmat({data.source}, [1 data.movieLength]), {framesCh1.name}, 'UniformOutput', false);

    mpath = [data.source 'Detection' filesep 'Masks' filesep];
    framesCh2 = dir([mpath '*.tif']);
    frameListCh2 = cellfun(@(x,y) [x y], repmat({mpath}, [1 data.movieLength]), {framesCh2.name}, 'UniformOutput', false);
end

frameCh1 = imread(frameListCh1{1});
frameCh2 = imread(frameListCh2{1});


% Choose default command line output for trackDisplay_GUI
handles.output = hObject;


% initialize handles
handles.frameListCh1 = frameListCh1;
handles.frameListCh2 = frameListCh2;


handles.data = data;
handles.tracks = varargin{2};
if numel(varargin)>2
    handles.tracks2 = varargin{3};
end

if isfield(data, 'channel1')

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
refreshDisplay(hObject, handles);



function refreshDisplay(hObject, handles)

f = get(handles.('slider1'), 'value');

% save zoom properties
haxes = handles.('axes1');
XLim = get(haxes, 'XLim');
YLim = get(haxes, 'YLim');

fr = handles.('frameListCh1');
frameCh1 = imread(fr{f});
fr = handles.('frameListCh2');
frameCh2 = imread(fr{f});
tracks = handles.('tracks');

startV = [tracks.start];
endV = [tracks.end];

idx = find(startV <= f & f <= endV);
handles.visibleIdx = idx;

% nActive = sum(idx);
% maxLifetime = 299; % change to read from input

% lft = handles.('lifetimes');
% 
% activeIdx = ~isnan(x(:,f));
% lifetimes = lft(activeIdx, 1:f);
% maxLifetime = max(lifetimes, [], 2);
% xs = x(activeIdx,1:f);
% ys = y(activeIdx,1:f);
cmap = jet(301);
% 
% 


imagesc(frameCh1, 'Parent', handles.('axes1'));
axis(handles.('axes1'),'image');
 

% hold(handles.('axes1'), 'on');
% for k = idx
%     fi = 1:f-tracks(k).start+1;
%     plot(tracks(k).x(fi), tracks(k).y(fi), '-', 'Color', cmap(tracks(k).lifetime,:));
% end
% hold(handles.('axes1'), 'off');


imagesc(frameCh2, 'Parent', handles.('axes2'));
axis(handles.('axes2'),'image');

hold(handles.('axes2'), 'on');
for k = idx
    fi = 1:f-tracks(k).start+1;
    plot(handles.('axes2'), tracks(k).x(1), tracks(k).y(1), '*', 'Color', cmap(tracks(k).lifetime,:));
    plot(handles.('axes2'), tracks(k).x(fi), tracks(k).y(fi), '-', 'Color', cmap(tracks(k).lifetime,:));
end
if isfield(handles, 'tracks2')
    tracks2 = handles.('tracks2');
    for k = idx
        fi = 1:f-tracks2(k).start+1;
        plot(tracks2(k).x(1), tracks2(k).y(1), 'o', 'Color', 'r');
        plot(tracks2(k).x(fi), tracks2(k).y(fi), '--', 'Color', 'r');
    end
end

% plot selected track marker
if ~isempty(handles.selectedTrack)
    t = tracks(handles.selectedTrack);
    ci = f-t.start+1;
    if 1 <= ci && ci <= t.lifetime
        plot(handles.('axes2'), t.x(ci), t.y(ci), 'ro', 'MarkerSize', 15);
    end
end

% show detection COM values
if get(handles.('checkbox1'), 'Value')
    detection = handles.('detection');
    plot(handles.('axes2'), detection(f).xcom, detection(f).ycom, 'x', 'Color', hsv2rgb([120/360 0.5 0.5]));
end
hold(handles.('axes2'), 'off');

% frame indicator
% ybounds = get(handles.('axes3'), 'YLim');
% axis(handles.('axes3'), [0 handles.data.movieLength ybounds(1) ybounds(2)]);
% plot(handles.('axes3'), [f ybounds(1)], [f ybounds(2)], 'k--');

% write zoom level
set(haxes, 'XLim', XLim);
set(haxes, 'YLim', YLim);
guidata(hObject,handles);



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

refreshDisplay(hObject, handles);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set focus for next input
axes(handles.('axes1')); % linked to axes2, selection possible in both
[x,y] = ginput(1);



% mean position of visible tracks
idx = handles.visibleIdx;

% idx = idx(1:10);

np = length(idx);

mu_x = zeros(1,length(np));
mu_y = zeros(1,length(np));

f = get(handles.('slider1'), 'value');

tracks = handles.('tracks');
for k = 1:np
    fi = 1:f-tracks(idx(k)).start+1; % was tracks(k)
    mu_x(k) = mean(tracks(idx(k)).x(fi));
    mu_y(k) = mean(tracks(idx(k)).y(fi));
end

% mu_x = tracks(idx(k)).x(fi);
% mu_y = tracks(idx(k)).y(fi);

% nearest point
d = sqrt((x-mu_x).^2 + (y-mu_y).^2);
minIdx = d==min(d);

% mu_x(idx)
% mu_y(idx)

handles.selectedTrack = idx(minIdx);
guidata(hObject,handles);
refreshDisplay(hObject, handles);


% % plot selected tracks
% h = handles.('axes2');
% hold(h, 'on');
% % plot(h, mu_x, mu_y, 'yo');
% plot(h, tracks(idx(minIdx)).x(fi), tracks(idx(minIdx)).y(fi), 'ro', 'MarkerSize', 15);
% hold(h, 'off');

sigma = 1.628;
h = handles.('axes3');
hold(h, 'off');
plot(h, tracks(idx(minIdx)).t, tracks(idx(minIdx)).A + tracks(idx(minIdx)).c, 'r');
hold(h, 'on');
plot(h, tracks(idx(minIdx)).t, tracks(idx(minIdx)).I/(2*pi*sigma^2) + tracks(idx(minIdx)).c, 'b--');
plot(h, tracks(idx(minIdx)).t, tracks(idx(minIdx)).c, 'k');
plot(h, tracks(idx(minIdx)).t, tracks(idx(minIdx)).c + 3*tracks(idx(minIdx)).cStd, 'k--');


% tracks(idx).I
% tracks(idx).t

xlim(handles.('axes3'), [0 handles.data.movieLength]);
legend(handles.('axes3'), 'Amplitude', 'Background');


% --- Executes on button press in montageButton.
function montageButton_Callback(hObject, eventdata, handles)
% hObject    handle to montageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.selectedTrack)
    tracks = handles.('tracks');
    t = tracks(handles.selectedTrack);
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
