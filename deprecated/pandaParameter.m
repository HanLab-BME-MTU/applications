function varargout = pandaParameter(varargin)
% PANDAPARAMETER M-file for pandaParameter.fig
%      PANDAPARAMETER, by itself, creates a new PANDAPARAMETER or raises the existing
%      singleton*.
%
%      H = PANDAPARAMETER returns the handle to a new PANDAPARAMETER or the handle to
%      the existing singleton*.
%
%      PANDAPARAMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PANDAPARAMETER.M with the given input arguments.
%
%      PANDAPARAMETER('Property','Value',...) creates a new PANDAPARAMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pandaParameter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pandaParameter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pandaParameter

% Last Modified by GUIDE v2.5 02-Jun-2009 14:47:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pandaParameter_OpeningFcn, ...
    'gui_OutputFcn',  @pandaParameter_OutputFcn, ...
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


% --- Executes just before pandaParameter is made visible.
function pandaParameter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pandaParameter (see VARARGIN)

% Choose default command line output for pandaParameter
handles.output = hObject;

if nargin > 2
    if strcmp(varargin{1},'myhandles')
        handles.parenthandles= varargin{2};
        %set the fields with the given parameters
        update_mask_axes(handles);
    end
end

% Update handles structure
guidata(hObject, handles);
set(handles.figure2,'Visible','on');
waitfor(handles.figure2,'Visible');
% UIWAIT makes pandaParameter wait for user response (see UIRESUME)
% uiwait(handles.figure2);


function parameters = update_mask_axes(handles)
thismask = imread([handles.parenthandles.directory_name filesep handles.parenthandles.maskinfo.maskdir(1).name]);

parameters.h100 = figure(100);
imshow(thismask, []);
title('Start Position Determination (red) and Spline Fitting (yellow)');
[rowSize, colSize] = size(thismask);
if isfield(handles,'v')
    axis(handles.v);
end

if get(handles.extrema, 'Value');
    if get(handles.orient_left,       'Value');
        parameters.ORIENT_CELL  = 0;
    elseif get(handles.orient_right,  'Value');
        parameters.ORIENT_CELL  = 1;
    elseif get(handles.orient_up,     'Value');
        parameters.ORIENT_CELL  = 2;
    elseif get(handles.orient_low,    'Value');
        parameters.ORIENT_CELL  = 3;
    end
elseif get(handles.crosshair, 'Value');
    if get(handles.orient_left,       'Value');
        parameters.ORIENT_CELL  = 4;
    elseif get(handles.orient_right,  'Value');
        parameters.ORIENT_CELL  = 5;
    elseif get(handles.orient_low,     'Value');
        parameters.ORIENT_CELL  = 6;
    elseif get(handles.orient_up,    'Value');
        parameters.ORIENT_CELL  = 7;
    end
else
    parameters.ORIENT_CELL  = 4;
end

if get(handles.optimized,     'Value');
    parameters.BESTSTART  = 1;
    parameters.ORIENT_CELL  = 4;    
elseif get(handles.follow,    'Value');
    parameters.BESTSTART  = 0;
elseif get(handles.fixed,     'Value');
    parameters.BESTSTART  = -1;
else
    parameters.BESTSTART  = 0;
end

tolerance_idx = get(handles.spline_tolerance,'Value');
parameters.TOLERANCE = (tolerance_idx-1)*10;

%%%%%  find the starting position, codes copied from protrusionAnalysis.m
L = bwlabel(thismask);
s  = regionprops(L, 'Area','Centroid');
Allarea=zeros(length(s),1);
for j=1:length(s)
    Allarea(j) = s(j).Area;
end
[tmp1, tmp2] = max(Allarea);
thismask = (L == tmp2);

% test if cell is isolated or touches the image border
bi1=find(thismask(1:rowSize,1)>0);
bi2=find(thismask(1:rowSize,colSize)>0);
bi3=find(thismask(1,1:colSize)>0);
bi4=find(thismask(rowSize,1:colSize)>0);

if isempty(bi1) & isempty(bi2) & isempty(bi3) & isempty(bi4)
    ISCLOSE = 1;
else
    ISCLOSE = 0;    
    set(handles.orient_left,  'Enable', 'off');
    set(handles.orient_right, 'Enable', 'off');
    set(handles.orient_up,    'Enable', 'off');
    set(handles.orient_low,   'Enable', 'off');
    set(handles.follow,       'Enable', 'off');
    set(handles.fixed,        'Enable', 'off');
    set(handles.extrema,       'Enable', 'off');
    set(handles.crosshair,     'Enable', 'off');
    set(handles.optimized,     'Enable', 'off');
end

if ISCLOSE
    pixel_list = prOrientEdge(thismask, parameters.ORIENT_CELL, 1);
else
    c = contourc(double(thismask), [0 0]);
    pixel_list = [c(1,2:end)' c(2,2:end)'];
end

dupPixel = find(sum(diff(pixel_list).^2,2) <0.0001);
pixel_list(dupPixel,:) = [];

if ~ISCLOSE
    % close contour
    c = bwboundaries(thismask);
    if isCurveClockwise([c{1}(:,2), c{1}(:,1)])
        pixel_list = pixel_list(end:-1:1,:);
    end
else
    % open contour
    if isCurveClockwise(pixel_list)
        % make the pixel Clockwise in the image coordinate
        pixel_list = pixel_list(end:-1:1,:);
    end
end

%%%%%  fit a spline, codes copied from prSamProtrusion.m
l2 = size(pixel_list,1);
if ISCLOSE == 1
    % For close contour, rotate the pixels so that the head of the splines is at the intersection of two splines
    pixel_list_extended = [pixel_list; pixel_list; pixel_list];
    sp_t = spaps([1:3*l2],pixel_list_extended',parameters.TOLERANCE);
    spline_pixel_list=fnval(sp_t,[l2+1:2*l2])';
else
    sp_t = spaps([1:l2],pixel_list',parameters.TOLERANCE);
    spline_pixel_list=fnval(sp_t,[1:l2])';
end
hold on;  plot(spline_pixel_list(:,1), spline_pixel_list(:,2),'y-o','MarkerSize',1.5,'MarkerFaceColor','y');  hold off;
% Display the starting pixel position
if parameters.ORIENT_CELL > -1
    hold on;  plot(spline_pixel_list(1,1), spline_pixel_list(1,2),'ro','MarkerFaceColor','r','MarkerSize',5);  hold off;
end

% --- Outputs from this function are returned to the command line.
function varargout = pandaParameter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

parameters = update_mask_axes(handles);
h100 = figure(100);
close(h100);
varargout{1} = parameters;
delete(hObject);


% % --- Executes on button press in Ok.
% function Ok_Callback(hObject, eventdata, handles)
% % hObject    handle to Ok (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % parameters = update_mask_axes(handles);
% % close(parameters.h100);
% 
% h100 = figure(100);
% close(h100);
% set(handles.figure2,'Visible','off');
% keyboard;

function pandaParameter_CloseRequestFcn(hObject, eventdata, handles)
% hide the figure here, it will be deleted in the OutputFcn
set(handles.figure2,'Visible','off');


function Ok_Callback(hObject, eventdata, handles)
pandaParameter_CloseRequestFcn(hObject, eventdata, handles);


% --- Executes on button press in orient_left.
function orient_left_Callback(hObject, eventdata, handles)
% hObject    handle to orient_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_left
set(handles.orient_left,      'Value', 1);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);

update_mask_axes(handles);

% --- Executes on button press in orient_right.
function orient_right_Callback(hObject, eventdata, handles)
% hObject    handle to orient_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_right
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 1);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);

update_mask_axes(handles);

% --- Executes on button press in orient_up.
function orient_up_Callback(hObject, eventdata, handles)
% hObject    handle to orient_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_up
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 1);
set(handles.orient_low,       'Value', 0);

update_mask_axes(handles);

% --- Executes on button press in orient_low.
function orient_low_Callback(hObject, eventdata, handles)
% hObject    handle to orient_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_low
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 1);

update_mask_axes(handles);

% --- Executes on button press in extrema.
function extrema_Callback(hObject, eventdata, handles)
% hObject    handle to extrema (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of extrema
set(handles.orient_left,  'Enable', 'on');
set(handles.orient_right, 'Enable', 'on');
set(handles.orient_up,    'Enable', 'on');
set(handles.orient_low,   'Enable', 'on');
set(handles.follow,    'Enable', 'on');
set(handles.fixed,   'Enable', 'on');

set(handles.crosshair,      'Value', 0);
set(handles.optimized,      'Value', 0);
update_mask_axes(handles);


% --- Executes on button press in crosshair.
function crosshair_Callback(hObject, eventdata, handles)
% hObject    handle to crosshair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of crosshair

set(handles.orient_left,  'Enable', 'on');
set(handles.orient_right, 'Enable', 'on');
set(handles.orient_up,    'Enable', 'on');
set(handles.orient_low,   'Enable', 'on');
set(handles.follow,    'Enable', 'on');
set(handles.fixed,   'Enable', 'on');

set(handles.extrema,       'Value', 0);
set(handles.optimized,     'Value', 0);
update_mask_axes(handles);


% --- Executes on button press in optimized.
function optimized_Callback(hObject, eventdata, handles)
% hObject    handle to optimized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optimized
set(handles.orient_left,  'Enable', 'off');
set(handles.orient_right, 'Enable', 'off');
set(handles.orient_up,    'Enable', 'off');
set(handles.orient_low,   'Enable', 'off');
set(handles.follow,    'Enable', 'off');
set(handles.fixed,   'Enable', 'off');

set(handles.extrema,       'Value', 0);
set(handles.crosshair,     'Value', 0);



% --- Executes on selection change in spline_tolerance.
function spline_tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to spline_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns spline_tolerance contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spline_tolerance

h100 = figure(100);
a = axis;
if a(1) == 0 & a(2) == 1 & a(3) == 0 & a(4) == 1
    close(h100);
    if isfield(handles,'v')
        handles = rmfield(handles,'v');
    end
else
    handles.v = axis;
end
update_mask_axes(handles);


% --- Executes during object creation, after setting all properties.
function spline_tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spline_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in follow.
function follow_Callback(hObject, eventdata, handles)
% hObject    handle to follow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of follow
set(handles.fixed,       'Value', 0);
set(handles.follow,       'Value', 1);

% --- Executes on button press in fixed.
function fixed_Callback(hObject, eventdata, handles)
% hObject    handle to fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixed
set(handles.fixed,       'Value', 1);
set(handles.follow,       'Value', 0);


% --- Executes during object creation, after setting all properties.
function figure2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


