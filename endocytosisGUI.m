function varargout = endocytosisGUI(varargin)
% ENDOCYTOSISGUI M-file for endocytosisGUI.fig
%      ENDOCYTOSISGUI, by itself, creates a new ENDOCYTOSISGUI or raises the existing
%      singleton*.
%
%      H = ENDOCYTOSISGUI returns the handle to a new ENDOCYTOSISGUI or the handle to
%      the existing singleton*.
%
%      ENDOCYTOSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENDOCYTOSISGUI.M with the given input arguments.
%
%      ENDOCYTOSISGUI('Property','Value',...) creates a new ENDOCYTOSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before endocytosisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to endocytosisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help endocytosisGUI

% Last Modified by GUIDE v2.5 26-Jan-2010 17:57:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @endocytosisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @endocytosisGUI_OutputFcn, ...
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

% --- Executes just before endocytosisGUI is made visible.
function endocytosisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to endocytosisGUI (see VARARGIN)
% Choose default command line output for endocytosisGUI
handles.output = hObject;
% Set initial parameters (these parameters are added here for the first
% time to the handles structure)
handles.iframe = 1;
handles.dragTailLength = 0;
handles.lftMin = 0;
handles.markFilter = 0;
handles.markInitiation = 0;
handles.markIternalization = 0;
handles.markGap = 0;
handles.markDragtail = 0;
handles.markPersistent = 0;
handles.markCutOff = 0;
handles.showHotSpots = 0;
handles.dragTail2Channel = 0;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes endocytosisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = endocytosisGUI_OutputFcn(hObject, eventdata, handles)
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
experiment.source = uigetdir([],'choose directory under which the images can be found');
%save directory for saving movie frames under a child directory of this
%direcotry
handles.source = experiment.source;
% load lifetime info
%it is loaded here to get the number of frames
load([experiment.source filesep 'LifetimeInfo' filesep 'lftInfo.mat']);
%set lifetime maximum to length of movie; done here for easy storage of
%movie length
handles.lftMax = size(lftInfo.Mat_xcoord,2);
%set text box on top of frame slider to read 'frame 1 out of
%movielength
set(handles.text1,'String', ['frame 1 out of ' num2str(handles.lftMax)]);
%update this value of edit box for lifetime max
set(handles.edit3,'String',num2str(handles.lftMax));
%FRAME SLIDER VALUES
%set frame slider max to last frame
set(handles.slider1,'Max',handles.lftMax);
%set frame slider min
set(handles.slider1,'Min',1);
%set the frame slider value
set(handles.slider1,'Value',1);
%set frame slider
set(handles.slider1,'SliderStep',[1/(handles.lftMax-1) 10/(handles.lftMax-1)]);
handles.iframe = 1;
%LOAD IMAGES
% get all tif files
filenames = dir([experiment.source filesep '*.tif']);
% load images for channel one
handles.imagesChannel1 = [];
for iIm = 1:length(filenames)
    handles.imagesChannel1(:,:,iIm) = imread([experiment.source filesep filenames(iIm).name]);
end
%set slider values for channel one max Intensity
set(handles.slider2,'Max',max(handles.imagesChannel1(:)));
set(handles.slider2,'Min',min(handles.imagesChannel1(:)));
set(handles.slider2,'Value',max(handles.imagesChannel1(:)));
handles.maxIntChannel1 = max(handles.imagesChannel1(:));
handles.minIntChannel1 = min(handles.imagesChannel1(:));
set(handles.slider2,'SliderStep',[1/(double(handles.maxIntChannel1)-double(handles.minIntChannel1)) ...
    10/(double(handles.maxIntChannel1)-double(handles.minIntChannel1))]);
%set slider values for channel one min Intensity
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
set(handles.text7,'String', experiment.source);
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

% DRAGTAIL FOR SECOND CHANNEL CHECKBOX
% --- Executes on button press in checkbox6.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox6
handles.dragTail2Channel = get(hObject,'Value');
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
% CREATE CODE TO MAKE AND SAVE .avi
%make avi for channel1
aviobj1 = avifile([handles.source filesep handles.movieDir '_channel1']);
for i = 1:size(handles.lftInfo.Mat_ycoord,2)
    handles.iframe = i;
    % Update handles structure
    guidata(hObject, handles);
    %plot image
    plotImage(handles);
    %make new invisible figure
    figureHandle = figure('Visible','off');
    %copy axes unto new figure
    newHandle = copyobj(handles.axes1,figureHandle);
    %get new axes handles
    axHandles = get(figureHandle,'CurrentAxes');
    %make greyscale
    colormap(axHandles,'gray');
    %print new figure
    %print(figureHandle,'-dtiff',[handles.source filesep handles.movieDir filesep 'channel1' filesep 'frame' num2str(i)])
    %make matlab movie frame
    aviobj1 = addframe(aviobj1,figureHandle);
    %close figure
    close(figureHandle)
    %if second channel is plotted
    if isfield(handles,'imagesChannel2')
        %make new invisible figure
        figureHandle = figure('Visible','off');
        %copy axes unto new figure
        newHandle = copyobj(handles.axes2,figureHandle);
        %get new axes handles
        axHandles = get(figureHandle,'CurrentAxes');
        %make greyscale
        colormap(axHandles,'gray');
        %print new figure
        print(figureHandle,'-dtiff',[handles.source filesep handles.movieDir filesep 'channel2' filesep 'frame' num2str(i)])
        %close figure
        close(figureHandle)
    end 
end
%close movie when done
aviobj1 = close(aviobj1);

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
%if user wants markers plotted on both channels
if handles.dragTail2Channel
    numChannels = 2;
else
    numChannels = 1;
end

for ichannel = 1:numChannels
    eval(['axHandle = handles.axes' num2str(ichannel) ';']);
    if handles.markInitiation
        plot(axHandle,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==1 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe), ...
            handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==1 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe),'g.');
    end
    if handles.markIternalization
        plot(axHandle,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==-1 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe), ...
            handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_disapp(:,handles.iframe)==-1 & ...
            handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
            & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe),'r.');
    end
    if handles.markGap
        plot(axHandle,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==4 | ...
            handles.lftInfo.Mat_status(:,handles.iframe)==5, handles.iframe), ...
            handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==4 | ...
            handles.lftInfo.Mat_status(:,handles.iframe)==5, ...
            handles.iframe),'m.');
    end
    if handles.dragTailLength ~= 0
        %had to do this because when X and Y are vectors instead of matrices
        %all the points are connected
        if handles.iframe == 1
            plot(axHandle,handles.Mat_xcoord(...
                handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax,handles.iframe)', ...
                handles.Mat_ycoord(handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, handles.iframe)','b.');
        elseif ~handles.showHotSpots
            plot(axHandle,handles.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
                handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
                max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
                handles.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
                handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
                max(handles.iframe-handles.dragTailLength,1):handles.iframe)','b-');
        else
            plot(axHandle,handles.Mat_xcoord(handles.hotOrNot==1 & ...
                handles.lftInfo.Mat_status(:,handles.iframe) ~=0 & ...
                handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
                max(handles.iframe-handles.dragTailLength,1):handles.iframe)', ...
                handles.Mat_ycoord(handles.hotOrNot==1 & ...
                handles.lftInfo.Mat_status(:,handles.iframe) ~=0 &...
                handles.lftInfo.Mat_lifetime(:,handles.iframe) >= lftMin ...
                & handles.lftInfo.Mat_lifetime(:,handles.iframe) <= lftMax, ...
                max(handles.iframe-handles.dragTailLength,1):handles.iframe)','m-');
            plot(axHandle,handles.Mat_xcoord(handles.hotOrNot==-1 & ...
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
        plot(axHandle,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==3, ...
            handles.iframe), ...
            handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==3, ...
            handles.iframe),'c.');
    end
    if handles.markCutOff
        plot(axHandle,handles.lftInfo.Mat_xcoord(handles.lftInfo.Mat_status(:,handles.iframe)==2, ...
            handles.iframe), ...
            handles.lftInfo.Mat_ycoord(handles.lftInfo.Mat_status(:,handles.iframe)==2, ...
            handles.iframe),'y.');
    end
    if handles.showHotSpots
        plot(axHandle,handles.hotSpotsX(:),handles.hotSpotsY(:),'y.')
        plot(axHandle,handles.hotSpotsX(:),handles.hotSpotsY(:),'y.')
    end
end %of for each channel being plotted
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
source = uigetdir([],'choose directory under which the images can be found');
% get all tif files
filenames = dir([source filesep '*.tif']);
% load image
handles.imagesChannel2 = [];
for iIm = 1:length(filenames)
    handles.imagesChannel2(:,:,iIm) = imread([source filesep filenames(iIm).name]);
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
%write current movie path in text box next to load movie button
set(handles.text8,'String', source);
% Update handles structure
guidata(hObject, handles);
%plot image
plotImage(handles);

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

% --- Executes during object creation, after setting all properties.
function text7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%MAXIMUM LIFETIME EDIT BOX
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.movieDir = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end