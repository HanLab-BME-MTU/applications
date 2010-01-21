function varargout = endocytosisGUI(varargin)
% ENDOCYTOSISGUI2PANNEL M-file for endocytosisGUI2Pannel.fig
%      ENDOCYTOSISGUI2PANNEL, by itself, creates a new ENDOCYTOSISGUI2PANNEL or raises the existing
%      singleton*.
%
%      H = ENDOCYTOSISGUI2PANNEL returns the handle to a new ENDOCYTOSISGUI2PANNEL or the handle to
%      the existing singleton*.
%
%      ENDOCYTOSISGUI2PANNEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENDOCYTOSISGUI2PANNEL.M with the given input arguments.
%
%      ENDOCYTOSISGUI2PANNEL('Property','Value',...) creates a new ENDOCYTOSISGUI2PANNEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before endocytosisGUI2Pannel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to endocytosisGUI2Pannel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help endocytosisGUI2Pannel

% Last Modified by GUIDE v2.5 15-Jan-2010 12:48:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @endocytosisGUI2Pannel_OpeningFcn, ...
    'gui_OutputFcn',  @endocytosisGUI2Pannel_OutputFcn, ...
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

% --- Executes just before endocytosisGUI2Pannel is made visible.
function endocytosisGUI2Pannel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to endocytosisGUI2Pannel (see VARARGIN)

% Choose default command line output for endocytosisGUI2Pannel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set initial parameters
handles.iframe = 1;
handles.dragTailLength = 0;
handles.lftMin = 0;
handles.lftMax = 0;
handles.markFilter = 0;
handles.markInitiation = 0;
handles.markIternalization = 0;
handles.markGap = 0;
handles.markDragtail = 0;
handles.markPersistent = 0;
handles.markCutOff = 0;
handles.showHotSpots = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes endocytosisGUI2Pannel wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = endocytosisGUI2Pannel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% FRAME NUMBER SLIDER
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% get frame number out of slider value
handles.iframe = round(get(hObject,'Value'));
%set text box next to dragtail input to the input value
set(handles.text1,'String', ['frame ' num2str(handles.iframe) ' out of ' num2str(size(handles.lftInfo.Mat_xcoord,2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% LOAD CHANNEL 1 PUSHBUTTON
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ask user to pick folder containing tifs for first channel
experiment.source = uigetdir;
%this doesn't work so well...idea is to display a busy text
set(handles.text2,'String','loading movie......please wait');
% load lifetime info
%it is loaded here to get the number of frames
load([experiment.source filesep 'LifetimeInfo' filesep 'lftInfo.mat']);
%set lifetime maximum to length of movie; done here for easy storage of
%movie length
handles.lftMax = size(lftInfo.Mat_xcoord,2);
%set text box next on top of frame slider to read 'frame 1 out of
%movielength
set(handles.text1,'String', ['frame 1 out of ' num2str(handles.lftMax)]);
%update this value of edit box for lifetime max
set(handles.edit3,'String',num2str(handles.lftMax));
%set slider max to last frame
set(handles.slider1,'Max',handles.lftMax);
%set min
set(handles.slider1,'Min',1);
set(handles.slider1,'Value',1);
set(handles.slider1,'SliderStep',[1/(handles.lftMax-1) 10/(handles.lftMax-1)]);
handles.iframe = 1;
% get all tif files
filenames = dir([experiment.source filesep '*.tif']);
% load images for channel one
handles.imagesChannel1 = [];
for iIm = 1:length(filenames)
    handles.imagesChannel1(:,:,iIm) = imread([experiment.source filesep filenames(iIm).name]);
end
%set slider values for channel one Intensity
set(handles.slider2,'Max',max(handles.imagesChannel1(:)));
set(handles.slider2,'Min',min(handles.imagesChannel1(:)));
set(handles.slider2,'Value',max(handles.imagesChannel1(:)));
handles.maxIntChannel1 = max(handles.imagesChannel1(:));
handles.minIntChannel1 = min(handles.imagesChannel1(:));
set(handles.slider2,'SliderStep',[1/(double(handles.maxIntChannel1)-double(handles.minIntChannel1)) ...
    10/(double(handles.maxIntChannel1)-double(handles.minIntChannel1))]);
%set slider values for channel two Intensity
set(handles.slider3,'Max',max(handles.imagesChannel1(:)));
set(handles.slider3,'Min',min(handles.imagesChannel1(:)));
set(handles.slider3,'Value',min(handles.imagesChannel1(:)));
set(handles.slider3,'SliderStep',[1/(double(handles.maxIntChannel1)-double(handles.minIntChannel1)) ...
    10/(double(handles.maxIntChannel1)-double(handles.minIntChannel1))]);
files = dir([experiment.source filesep 'ClusterData' filesep 'clusterResults*']);
%Cluster Results For Hot Spots
%load cluster results
if exist([experiment.source filesep 'ClusterData']) == 7
    load([experiment.source filesep 'ClusterData' filesep files(1).name]);
    %save cluster results into handles structure
    handles.clusterResults = clusterResults;
    handles.pitID = retievePitIDFromClusterResults(lftInfo,clusterResults);
    handles.hotOrNot = zeros(size(lftInfo.Mat_xcoord,1),1);
    handles.hotOrNot(handles.pitID(clusterResults.clusterResults(:,3)~=0)) = 1;
    handles.hotOrNot(handles.pitID(clusterResults.clusterResults(:,3)==0)) = -1;
    %make points in a circle around each
    %make nump in a circle around hot spot centroids
    nump = 8;
    %equally space them
    theta = linspace(0,2*pi,nump);
    rho = ones(1,nump)*clusterResults.hotSpotRadius;
    [x,y] = pol2cart(theta,rho);
    %add these points to each hotspot centroid to denote hotspot
    for ihotspot = 1:size(clusterResults.clusterCentroids,1)
        X(ihotspot,:) = x + repmat(clusterResults.clusterCentroids(ihotspot,1),1,nump);
        Y(ihotspot,:) = y + repmat(clusterResults.clusterCentroids(ihotspot,2),1,nump);
    end
    %save points to handles structure
    handles.hotSpotsX = X ;
    handles.hotSpotsY = Y ;
end %of if cluster results exist
%make zeros into nans so that the dragtails can be plotted
%sparse matrices take too long to change into nans
Mat_xcoord = full(lftInfo.Mat_xcoord);
Mat_xcoord(Mat_xcoord==0) = nan;
handles.Mat_xcoord = Mat_xcoord;
Mat_ycoord = full(lftInfo.Mat_ycoord);
Mat_ycoord(Mat_ycoord==0) = nan;
handles.Mat_ycoord = Mat_ycoord;
handles.lftInfo = lftInfo;
%plot image
plotImage(handles);
%write current movie path in text box next to load movie button
set(handles.text2,'String', experiment.source);
% Update handles structure
guidata(hObject, handles);

% INITIATION MARKER CHECKBOX
% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles.markInitiation = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

% INTERNALIZATION MARKER CHECKBOX
% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox2
handles.markIternalization = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

% GAP MARKER CHECKBOX
% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox3
handles.markGap = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

% PERSISTENT MARKER CHECKBOX
% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox5
handles.markPersistent = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

% CUT-OFF MARKER CHECKBOX
% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox6
handles.markCutOff = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

%DRAGTAIL LENGHT EDIT BOX
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a
%        double
handles.dragTailLength = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%MINIMUM LIFETIME EDIT BOX
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.lftMin = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles)
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%MAXIMUM LIFETIME EDIT BOX
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.lftMax = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles)
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% MAKE MOVIE PUSHBUTTON
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% CREATE CODE TO MAKE AND SAVE .dvi
for i = 1:size(handles.lftInfo.Mat_ycoord,2)
    handles.iframe = i;
    % Update handles structure
    guidata(hObject, handles);
    %plot image
    plotImage(handles);
end

% CHANNEL ONE AXES
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%TEXT ABOVE FRAME SLIDER
% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% PLOT IMAGE
function plotImage(handles)
%replacing makes the plotting go much faster
set([handles.axes1 handles.axes2],'NextPlot','replacechildren')
%plot image .... second channel
%if loaded
if isfield(handles,'imagesChannel2')
    imshow(handles.imagesChannel2(:,:,handles.iframe),[handles.minIntChannel2 handles.maxIntChannel2],'Parent',handles.axes2);
end
%plot image ..... first channel
imshow(handles.imagesChannel1(:,:,handles.iframe),[handles.minIntChannel1 handles.maxIntChannel1],'Parent',handles.axes1);
%put hold on on channel one axes to plot markers on top
hold(handles.axes1,'on')
hold(handles.axes2,'on')
% turn zoom on for axes of channel one
zoom(handles.axes1,'on')
%plot markers
%if markers are to be filtered according to lifetime
lftMax = str2num(get(handles.edit3,'String'));
lftMin = str2num(get(handles.edit2,'String'));

if handles.markInitiation
    plot(handles.axes1,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==1 & ...
        handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
        & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe), ...
        handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==1 & ...
        handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
        & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe),'g.');
end
if handles.markIternalization
    plot(handles.axes1,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==-1 & ...
        handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
        & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe), ...
        handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==-1 & ...
        handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
        & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe),'r.');
end
if handles.markGap
    plot(handles.axes1,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==4 | ...
        handles.lftInfo.Mat_status(:,handles.iframe)==5, handles.iframe), ...
        handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==4 | ...
        handles.lftInfo.Mat_status(:,handles.iframe)==5, ...
        handles.iframe),'m.');
end
if handles.dragTailLength ~= 0
    %had to do this because when X and Y are vectors instead of matrices
    %all the points are connected
    if handles.iframe == 1
        plot(handles.axes1,handles.Mat_xcoord(...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe)', ...
            handles.Mat_ycoord(handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, handles.iframe)','b.');
    elseif ~handles.showHotSpots
        plot(handles.axes1,handles.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
            handles.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)','b-');
    else
                plot(handles.axes1,handles.Mat_xcoord(handles.hotOrNot==1 & ...
                    handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
            handles.Mat_ycoord(handles.hotOrNot==1 & ...
            handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)','m-');
        plot(handles.axes1,handles.Mat_xcoord(handles.hotOrNot==-1 & ...
                    handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
            handles.Mat_ycoord(handles.hotOrNot==-1 & ...
            handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)','c-');
         plot(handles.axes2,handles.Mat_xcoord(handles.hotOrNot==1 & ...
                    handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
            handles.Mat_ycoord(handles.hotOrNot==1 & ...
            handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)','m-');
        plot(handles.axes2,handles.Mat_xcoord(handles.hotOrNot==-1 & ...
                    handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
            handles.Mat_ycoord(handles.hotOrNot==-1 & ...
            handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
            max(handles.iframe-handles.dragTailLength,1):handles.iframe)','c-');
    end %of if first frame
end
if handles.markPersistent
    plot(handles.axes1,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==3, ...
        handles.iframe), ...
        handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==3, ...
        handles.iframe),'c.');
end
if handles.markCutOff
    plot(handles.axes1,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==2, ...
        handles.iframe), ...
        handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==2, ...
        handles.iframe),'y.');
end
if handles.showHotSpots
    plot(handles.axes1,handles.hotSpotsX(:),handles.hotSpotsY(:),'y.')
    plot(handles.axes2,handles.hotSpotsX(:),handles.hotSpotsY(:),'y.')
end
% link the zoom factor for axes of channel one and channel two
linkaxes([handles.axes1 handles.axes2])

% CHANNEL ONE MAX INTENSITY SLIDER
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.maxIntChannel1 = get(hObject,'Value');
%set text box next to dragtail input to the input value
%set(handles.text1,'String', ['frame ' num2str(handles.iframe) ' out of ' num2str(size(handles.lftInfo.Mat_xcoord,2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% CHANNEL ONE MIN INTENSITY SLIDER
% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.minIntChannel1 = get(hObject,'Value');
%set text box next to dragtail input to the input value
%set(handles.text1,'String', ['frame ' num2str(handles.iframe) ' out of ' num2str(size(handles.lftInfo.Mat_xcoord,2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% CHANNEL 2 AXES
% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes2

% CHANNEL 2 MAX INTENSITY SLIDER
% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.maxIntChannel2 = get(hObject,'Value');
%set text box next to dragtail input to the input value
%set(handles.text1,'String', ['frame ' num2str(handles.iframe) ' out of ' num2str(size(handles.lftInfo.Mat_xcoord,2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% CHANNEL 2 MIN INTENSITY SLIDER
% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.minIntChannel2 = get(hObject,'Value');
%set text box next to dragtail input to the input value
%set(handles.text1,'String', ['frame ' num2str(handles.iframe) ' out of ' num2str(size(handles.lftInfo.Mat_xcoord,2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% LOAD CHANNEL 2 PUSHBUTTON
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
experiment.source = uigetdir;
% get all tif files
cd(experiment.source)
filenames = dir(['*.tif']);
% load image
handles.imagesChannel2 = [];
for iIm = 1:length(filenames)
    handles.imagesChannel2(:,:,iIm) = imread([experiment.source filesep filenames(iIm).name]);
end
%set slider values for max Intensity
set(handles.slider4,'Max',max(handles.imagesChannel2(:)));
set(handles.slider4,'Min',min(handles.imagesChannel2(:)));
set(handles.slider4,'Value',max(handles.imagesChannel2(:)));
handles.maxIntChannel2 = max(handles.imagesChannel2(:));
handles.minIntChannel2 = min(handles.imagesChannel2(:));
set(handles.slider4,'SliderStep',[1/(double(handles.maxIntChannel2)-double(handles.minIntChannel2)) ...
    10/(double(handles.maxIntChannel2)-double(handles.minIntChannel2))]);
%set slider values for min Intnsity
set(handles.slider5,'Max',max(handles.imagesChannel2(:)));
set(handles.slider5,'Min',min(handles.imagesChannel2(:)));
set(handles.slider5,'Value',min(handles.imagesChannel2(:)));
set(handles.slider5,'SliderStep',[1/(double(handles.maxIntChannel2)-double(handles.minIntChannel2)) ...
    10/(double(handles.maxIntChannel2)-double(handles.minIntChannel2))]);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);
%write current movie path in text box next to load movie button
set(handles.text14,'String', experiment.source);


% --- Executes during object creation, after setting all properties.
function text14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% HOT SPOT MARKER CHECKBOX
% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showHotSpots = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of checkbox8
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);