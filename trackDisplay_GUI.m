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

% Last Modified by GUIDE v2.5 19-Sep-2010 19:03:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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

frames = dir([varargin{1}.source '*.tif']);
frameList = cellfun(@(x,y) [x y], repmat({varargin{1}.source}, [1 varargin{1}.movieLength]), {frames.name}, 'UniformOutput', false);

mpath = [varargin{1}.source 'Detection' filesep 'Masks' filesep];
masks = dir([mpath '*.tif']);
maskList = cellfun(@(x,y) [x y], repmat({mpath}, [1 varargin{1}.movieLength]), {masks.name}, 'UniformOutput', false);


% Choose default command line output for trackDisplay_GUI
handles.output = hObject;
handles.frameList = frameList;
handles.maskList = maskList;
handles.tracks = varargin{2};


% Update handles structure
guidata(hObject, handles);

% set(handles.frameList,'String',frameList);
% load([varargin{1}.source 'TrackInfoMatrices' filesep 'trackinfo.mat']);

h = handles.('slider1');
set(h, 'Min', 1);
set(h, 'Max', varargin{1}.movieLength);
set(h, 'SliderStep', [1/(varargin{1}.movieLength-1) 0.05]);

frame = imread(frameList{1});
mask = imread(maskList{1});

imagesc(frame, 'Parent', handles.('axes1'));
imagesc(mask, 'Parent', handles.('axes2'));
colormap(gray(256));
linkaxes([handles.('axes1') handles.('axes2')]);
axis([handles.('axes1') handles.('axes2')],'image');

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

set(hObject, 'Value', round(get(hObject, 'value')));
% guidata(gcbo) % same as 'handles'

% save zoom properties
haxes = handles.('axes1');
XLim = get(haxes, 'XLim');
YLim = get(haxes, 'YLim');


f = get(hObject, 'Value');
fr = handles.('frameList');
frame = imread(fr{f});
fr = handles.('maskList');
mask = imread(fr{f});

tracks = handles.('tracks');

startV = [tracks.start];
endV = [tracks.end];

idx = find(startV <= f & f <= endV);
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

% axes(handles.('axes1'));
imagesc(frame, 'Parent', handles.('axes1'));
axis(handles.('axes1'),'image');

hold on;
for k = idx
    fi = 1:f-tracks(k).start+1;
    plot(tracks(k).x(fi), tracks(k).y(fi), '-', 'Color', cmap(tracks(k).lifetime,:));
end
hold off;


% axes(handles.('axes2'));
imagesc(mask, 'Parent', handles.('axes2'));
axis(handles.('axes2'),'image');
hold on;
for k = idx
    fi = 1:f-tracks(k).start+1;
    plot(tracks(k).x(fi), tracks(k).y(fi), '-', 'Color', cmap(tracks(k).lifetime,:));
end
hold off;

% write zoom level
set(haxes, 'XLim', XLim);
set(haxes, 'YLim', YLim);




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
