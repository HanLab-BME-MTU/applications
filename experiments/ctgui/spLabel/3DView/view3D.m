function varargout  =  view3D(varargin)
% VIEW3D M-file for view3D.fig
%      VIEW3D, by itself, creates a new VIEW3D or raises the existing
%      singleton*.
%
%      H = VIEW3D returns the handle to a new VIEW3D or the handle to
%      the existing singleton*.
%
%      VIEW3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW3D.M with the given input arguments.
%
%      VIEW3D('Property','Value',...) creates a new VIEW3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view3D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view3D

% Last Modified by GUIDE v2.5 20-Mar-2003 08:57:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view3D_OpeningFcn, ...
                   'gui_OutputFcn',  @view3D_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before view3D is made visible.
function view3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view3D (see VARARGIN)

% Choose default command line output for view3D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');

%immediately save figure handle
SetUserData(labelPanelH,hObject,1,'view3DGUIH');

%if there is no labelpanel: exit
if isempty(labelPanelH)
    close(hObject);
    error('no movie loaded (start labelgui)')
    return
end
labelguiH = findall(0,'Tag','labelgui');
label_handles = guihandles(labelguiH); %handle structure of labelgui

%initialize timesliders
idlist = GetUserData(labelPanelH,'idlist');
maxTime = size(idlist,2);
set(handles.view3D_tstart_slider,'Max',maxTime);
set(handles.view3D_tend_slider,'Max',maxTime);
set(handles.view3D_tstart_slider,'Min',1);
set(handles.view3D_tend_slider,'Min',1);

%set timesliders to labelgui time - independently of setTime
labelguiTime = get(label_handles.slider3,'Value');
%if nothing in labelgui: look for good time to plot. (t+1,t-1,t+2,t-2 etc) 
step = 1;
done = 0;
while ~done
    if labelguiTime==0|labelguiTime==maxTime+1
        done = 1;
    elseif ~isempty(idlist(labelguiTime).linklist)
        done = 1;
    else
        labelguiTime = labelguiTime+step;
        step = -(step+sign(step));
    end
end
%if still empty: start looking in only one direction (into which step is already pointing)
if labelguiTime==maxTime+1|labelguiTime==0
    labelguiTime = labelguiTime+step;
    while labelguiTime<maxTime&labelguiTime>1&isempty(idlist(labelguiTime).linklist)
        labelguiTime = labelguiTime+sign(step);
    end
end
%if still empty: bad idlist
if isempty(idlist(labelguiTime).linklist)
    error('bad idlist - all empty')
    close(hObject)
end
set(label_handles.slider3,'Value',labelguiTime);
labelgui('slider3_Callback',label_handles.slider3,[],label_handles);
set(handles.view3D_tstart_slider,'Value',labelguiTime);
set(handles.view3D_tend_slider,'Value',labelguiTime);
set(handles.view3D_tstart_txt,'String',num2str(labelguiTime));
set(handles.view3D_tend_txt,'String',num2str(labelguiTime));

%remember current status of PD-menus (some might not be needed)
handles.lastAlign = get(handles.view3D_align_PD,'Value');
handles.lastCenter = get(handles.view3D_center_PD,'Value');
handles.lastLinks = get(handles.view3D_link_PD,'Value');

%predefine values
handles.lastLabelguiTime = labelguiTime;
handles.lastSetTime = get(handles.view3D_setTime_PD,'Value');
handles.centerTagNumber = 1; %to prevent errors if default is set to center tag
handles.axisDefined = 0;
handles.lastAxis = 0;
handles.newAxis = []; %should not be necessary, as this is only checked if an axis has been defined before, but you never know...
handles.rotMatrix = [];

%save data
guidata(hObject,handles);

%generate handle structure; generateHandles then calls refresh.
view3D_generateHandles;



% --- Outputs from this function are returned to the command line.
function varargout = view3D_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function view3D_tstart_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_tstart_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function view3D_tstart_slider_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_tstart_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%make sure only integer times are set 
tstart = round(get(gcbo,'Value'));
set(gcbo,'Value',tstart);
%set time_text
set(handles.view3D_tstart_txt,'String',num2str(tstart));

%adjust tend, if necessary
tend = get(handles.view3D_tend_slider,'Value');
if tstart>tend
    set(handles.view3D_tend_slider,'Value',tstart);
    set(handles.view3D_tend_txt,'String',num2str(tstart));
end

view3D_refresh;

% --- Executes during object creation, after setting all properties.
function view3D_tend_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_tend_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function view3D_tend_slider_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_tend_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%make sure only integer times are set 
tend = round(get(gcbo,'Value'));
set(gcbo,'Value',tend);

%set time_text
set(handles.view3D_tend_txt,'String',num2str(tend));

%adjust tend, if necessary
tstart = get(handles.view3D_tstart_slider,'Value');
if tstart>tend
    set(handles.view3D_tstart_slider,'Value',tend);
    set(handles.view3D_tstart_txt,'String',num2str(tend));
end

view3D_refresh;


% --- Executes during object creation, after setting all properties.
function view3D_link_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_link_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in view3D_link_PD.
function view3D_link_PD_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_link_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%just update plot
view3D_refresh;

% --- Executes during object creation, after setting all properties.
function view3D_center_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_center_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in view3D_center_PD.
function view3D_center_PD_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_center_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currValue = get(hObject,'Value');
handles = guidata(hObject); %get complete handle structure

switch currValue
    case 1 %no center
        set(hObject,'ForegroundColor','k');
        handles.lastCenter = 1;
    case 2 %centroid
        set(hObject,'ForegroundColor','k');
        handles.lastCenter = 2;
    case 3 %tag
        lastTag = handles.centerTagNumber;
        h = selectTagGUI;
        sHandles = guidata(h);
        sHandles.num2select = 1;
        guidata(h,sHandles);
        uiwait(h);
        %load latest handles
        handles = guidata(hObject);
        if handles.centerTagNumber==0 %user canceled: reset
            set(hObject,'Value',handles.lastCenter);
            handles.centerTagNumber = lastTag;
            guidata(hObject,handles);
            return %end evaluation here
        else
            %find colormap
            cMap = handles.cMap;
            cMapFact = handles.cMapFact;
            %calculate color
            fgCol = cMap(2^(handles.centerTagNumber-1)*cMapFact,:);
            %set foregroundColor
            set(hObject,'ForegroundColor',fgCol);
            handles.lastCenter = 3;
        end
end

%save data and update plot
guidata(hObject,handles);
view3D_generateHandles;
        

% --- Executes during object creation, after setting all properties.
function view3D_align_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_align_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in view3D_align_PD.
function view3D_align_PD_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_align_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
currValue = get(hObject,'Value');
handles = guidata(hObject); %get complete handle structure
 
switch currValue
    case 1 %no align
        if handles.lastAlign==1;
            %don't change anything
        else
            handles.lastAlign = 1;
            guidata(hObject,handles);
            view3D_generateHandles;
        end
        
    case 2 %align axis
        %check if axis has already been defined; else define axis
        if handles.axisDefined==0
            view3D_defineAxis_PB_Callback(handles.view3D_defineAxis_PB, eventdata, handles);
            handles = guidata(hObject);
            %if user cancelled: handles.axisTagNumber==0 -> stop evaluation
            if handles.axisTagNumber==0 %user canceled: reset
                handles.axisTagNumber = handles.lastAxis;
                set(hObject,'Value',handles.lastAlign);
                guidata(hObject,handles);
                return %end evaluation here
            end
        end
        
        %store data
        handles.lastAlign = 2;
        guidata(hObject,handles);
        
        view3D_calcAxis(0);
        %calculate axis vectors and rotation matrices if necessary
        if handles.newAxis
            
            view3D_calcRotMat;
        end
        
        %generate new handle structure
        view3D_generateHandles;
        
        %warn user if not correctly centered
        if ~any(handles.axisTagNumber==handles.centerTagNumber)|get(handles.view3D_center_PD,'Value')~=3
            h = warndlg('The plot is not centered on one of the end of the axis; it might therefore look misaligned.','No bug!');
            uiwait(h);
        end
         
    case 3 %align metaphase coordinates
        
        %check if axis has already been defined; else define axis
        if handles.axisDefined==0
            view3D_defineAxis_PB_Callback(handles.view3D_defineAxis_PB, eventdata, handles);
            handles = guidata(hObject);
            %if user cancelled: handles.axisTagNumber==0 -> stop evaluation
            if handles.axisTagNumber==0 %user canceled: reset
                handles.axisTagNumber = handles.lastAxis;
                set(hObject,'Value',handles.lastAlign);
                guidata(hObject,handles);
                return %end evaluation here
            end
        end
        
        %store data
        handles.lastAlign = 3;
        guidata(hObject,handles);
        
        view3D_calcAxis(0);
        %calculate align axis and rotation matrices if neccessary
        if handles.newAxis
            
            view3D_calcRotMat;
        end
        
        %calculate metaphase-align-rotation matrices
        view3D_calcMetaphaseAlign;
        
        %generate new handle structure
        view3D_generateHandles;
        
        %warn user if not correctly centered
        if ~any(handles.axisTagNumber==handles.centerTagNumber)|get(handles.view3D_center_PD,'Value')~=3
            h = warndlg('The plot is not centered on one of the end of the axis; it might therefore look misaligned.','No bug!');
            uiwait(h);
        end
end

% --- Executes during object creation, after setting all properties.
function view3D_check_drawAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3D_check_drawAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in view3D_check_drawAxis.
function view3D_check_drawAxis_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_check_drawAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

view3D_refresh;


% --- Executes during object creation, after setting all properties.
function view3D_setTime_PD_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in view3D_setTime_PD.
function view3D_setTime_PD_Callback(hObject, eventdata, handles)

%just update plot
view3D_refresh;


 


% --- Executes on button press in view3D_defineAxis_PB.
function view3D_defineAxis_PB_Callback(hObject, eventdata, handles)
% hObject    handle to view3D_defineAxis_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

%check if more than one tag (necessary to define axis!)
if handles.nTags<2
    h = warndlg('Can''t define axis: not enough tags!','bad data');
    uiwait(h);
    return
    %end eval here
end

%remember old state (in case of cancel)
lastAxis = handles.lastAxis;

%launch GUI to select two tags to define axis
h = selectTagGUI;
sHandles = guidata(h);
sHandles.num2select = 2;
guidata(h,sHandles);
uiwait(h);
%load latest handles
handles = guidata(hObject);
if handles.axisTagNumber==0 %user canceled: reset
    handles.axisTagNumber = lastAxis;
    guidata(hObject,handles);
    return %end evaluation here
else %show axis, turn on option to draw axis
    %show colors
    %find colormap
    cMap = handles.cMap;
    cMapFact = handles.cMapFact;
    %calculate color
    frameCol = cMap(2.^(handles.axisTagNumber-1)*cMapFact,:);
    %show colors
    set(handles.view3D_alignFrame1,'BackgroundColor',frameCol(1,:),'Visible','on');
    set(handles.view3D_alignFrame2,'BackgroundColor',frameCol(2,:),'Visible','on');
    %turn on checkbox
    set(handles.view3D_check_drawAxis,'Enable','on');
    
    handles.lastAxis = handles.axisTagNumber;
    %remember that a new axis has been defined
    handles.newAxis = 1;
    
    guidata(hObject,handles);
end

%calculate axis
if gcbo==handles.view3D_align_PD
    %don't calculate axis -> it will be done in the other callback
else
    view3D_calcAxis(1);
end


handles = guidata(hObject);


if (get(handles.view3D_check_drawAxis,'Value')==1)
    view3D_refresh;
end